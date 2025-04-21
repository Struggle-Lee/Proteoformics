##version1: OLR+amino acid
#Prepare packages
library(data.table)
library(stringr)
library(intervals)
library(mclust)
library(nlme)
library(ggplot2)
library(openxlsx)
library(parallel)
library(doParallel)
library(dplyr)
library(extrafont)
library(extrafontdb)

massdata<-fread("/public/home/yangjh/psm.txt")

#111
massdata[, AminoAcid_Mod := str_extract(Modified.Peptide,"[A-Z]{1}(?=\\[)")]
massdata.pre<-massdata
massdata[,PeptidePOS:=str_extract(massdata$Modifications,"[0-9]{1,2}(?=\\()")]
massdata[,PeptideL:=str_length(Sequence)]
massdata<-massdata[PeptidePOS==1 | PeptidePOS==PeptideL]
#222
massdata[,PeptidePOS:=str_extract(massdata$Modifications,"[0-9]{1,2}(?=\\()")]
massdata[,PeptideL:=str_length(Sequence)]
massdata<-massdata[PeptidePOS!=1 & PeptidePOS!=PeptideL]


#Make output dictionary
directoryName <- "C:/Users/202005/Desktop/11"
# output file directory
output_directory <- str_c(directoryName, "/ncaasis_multi_output_data")
if (dir.exists(output_directory))
  file.remove(str_c(output_directory, "/", list.files(str_c(output_directory, "/"))))
if (!dir.exists(output_directory))
  dir.create(output_directory)

# log output directory
log_directory <- str_c(directoryName, "/log")
if (!dir.exists(log_directory))
  dir.create(log_directory)

# log file
log_file <- file(str_c(log_directory, "/", "ncaasis_multi.log"))
# log function

# log function
Log <- function(text, ...) {
  msg <- sprintf(paste(as.character(Sys.time()), ": ", text, "\n"), ..., collapse= ',')
  cat(msg, file = str_c(log_directory, "/", "ncaasis_multi.log"), append = T)
}

print("配置 completed")

minDeltaMass <- round(min(massdata$DeltaMass))
maxDeltaMass <- round(max(massdata$DeltaMass))
print(minDeltaMass )
print(maxDeltaMass )

## Normal fitting function definition
f <- function(x, meanValue, sdValue){
  return ((1 / (sqrt(2 * pi) * sdValue)) * exp(-((x-meanValue)*(x-meanValue)) / (2*sdValue*sdValue)))}

## nls iteration
nlsfunc<-function(intvMods,h){
nlsResult <- tryCatch(
  {
    # parameter estimate
    nls(y~f(x,m,s),
        data = data.frame(x = h$mids,y = h$density), 
        control = list(maxiter = 50, tol = 1e-05, minFactor = 1/8096,
                       printEval = FALSE, warnOnly = FALSE), 
        start = list(m = mean(h$mids), 
                     s = sd(h$mids)))
  },
  error = function(err) {
    Log("first nls error:")
    Log(err)
  },
  warning = function(war) {
    Log("first nls warning:")
    Log(war)
  },
  finally = {
  })  

if (is.null(nlsResult) | is.numeric(nlsResult)) 
{
  nlsResult <- NULL
  nlsResult <- tryCatch(
    {
      nls(y~f(x,m,s),
          data = data.frame(x = h$mids,y = h$density), 
          control = list(maxiter = 100, tol = 1e-05, minFactor = 1/8096,
                         printEval = FALSE, warnOnly = FALSE), 
          start = list(m = mean(h$mids), 
                       s = sd(h$mids) / 2))
    },
    error = function(err) {
      Log("second nls error:")
      Log(err)
    },
    warning = function(war) {
      Log("second nls warning:")
      Log(war)
    },
    finally = {
    })
  if (is.null(nlsResult) | is.numeric(nlsResult)) 
  {
    nlsResult <- NULL
    nlsResult <- tryCatch(
      {
        nls(y~f(x,m,s),
            data = data.frame(x = h$mids,y = h$density), 
            control = list(maxiter = 200, tol = 1e-05, minFactor = 1/8096,
                           printEval = FALSE, warnOnly = FALSE), 
            start = list(m = mean(h$mids), 
                         s = sd(h$mids) / 4))
      },
      error = function(err) {
        Log("third nls error:")
        Log(err)
      },
      warning = function(war) {
        Log("third nls warning:")
        Log(war)
      },
      finally = {
      })    
  }else {
    Log("second NLS success.")
  }
}else {
  Log("first NLS success.")
}
return(nlsResult)
}

### definition of OLR calculating function
G<-function(x,miu,sigma){return(exp(-1/2*((x-miu)/sigma)^2)*(1/(sqrt(2*pi)*sigma)))}
pdf<-function(x,miu1,miu2,sigma1,sigma2,p1,p2){return(p1*G(x,miu1,sigma1)+p2*G(x,miu2,sigma2))}
OLR<-function(miu1,miu2,sigma1,sigma2,p1,p2,lidu){
  index<-1
  step<-round(abs((miu2-miu1))/lidu,0)
  if(miu2>miu1){
  pdf.pre<-pdf(miu2,miu1,miu2,sigma1,sigma2,p1,p2)
  pdf.next<-NULL
  miu<-miu2
  pdfc<-pdf.pre
  while (miu>miu1) {miu<-miu-lidu
  pdf.next<-(pdf(miu,miu1,miu2,sigma1,sigma2,p1,p2))
  pdfc<-c(pdfc,pdf.next)
  }    
  Overlaprate<-min(pdfc)/min(pdf(miu1,miu1,miu2,sigma1,sigma2,p1,p2),pdf(miu2,miu1,miu2,sigma1,sigma2,p1,p2))
  }else if(miu2=miu1){
    Overlaprate<-1
  }else{
    pdf.pre<-pdf(miu1,miu2,miu1,sigma2,sigma1,p2,p1)
    pdf.next<-NULL
    miu<-miu1
    pdfc<-pdf.pre
    while (miu>miu2) {miu<-miu-lidu
    pdf.next<-(pdf(miu,miu1,miu2,sigma1,sigma2,p1,p2))
    pdfc<-c(pdfc,pdf.next)
    }    
    Overlaprate<-min(pdfc)/min(pdf(miu1,miu1,miu2,sigma1,sigma2,p1,p2),pdf(miu2,miu1,miu2,sigma1,sigma2,p1,p2))
  }
  return(Overlaprate)
}

### ploting parameters
ggtheme=theme_light()+theme(legend.title = element_blank(),
                            axis.text.y = element_text(colour="black", size=14,face = "bold",family = "serif"),
                            axis.text.x = element_text(colour="black", size=14,face = "bold",family = "serif",angle = 60,hjust = 1),
                            axis.title = element_text(colour="black", size=14,face = "bold",family = "serif"),
                            legend.position = "right",title =element_text(size=14,face = "bold",family = "serif"),
                            legend.text = element_text(size=14,face = "bold",family = "serif"))

### definition of ploting
plot.clust<-function(intervalmods,mclustResult,binwidth,nd,output_directory,intervalname,intv,label){
  p<-ggplot()+ggtheme+ylab("PSM count")+geom_histogram(intervalmods,mapping = aes(x=DeltaMass,y=..count..),binwidth = binwidth,stat="bin",fill="gray50",alpha=0.9)
  a<-ggplot_build(p)$data[[1]]$count%>%max()
  p<-p+stat_density(intervalmods,mapping = aes(x=DeltaMass,y=after_stat(scaled*a)),colour="black",fill="transparent",n=nd)+ggtitle(label = label)
  U <- sort(unique(mclustResult[["classification"]]))
  L <- length(U)
  bd<-seq(from=0-a/L,to=a,length.out=3*L+4)
  df<-data.frame(DM=intervalmods$DeltaMass,cls=mclustResult$classification)
  p<-p+geom_segment(mapping = aes(x=DM,xend=DM,y=bd[3*cls+2],yend=bd[3*cls+3],color=as.factor(cls)),data=df,show.legend = FALSE)+scale_color_manual(values=mclust.options("classPlotColors")[1:L])
  p<-p+geom_segment(mapping = aes(x=DM,xend=DM,y=bd[2],yend=bd[3]),data=df,colour="black",show.legend = FALSE)+geom_line(data = df,aes(x=min(DM),y=0,color=as.factor(cls)),size=5)
  int<-as.data.frame(intv)
  if(nrow(int)>0){
    p<-p+geom_point(intv%>%as.vector()%>%as.data.frame(),mapping = aes(.,y=a/2),fill="black",size=2)+geom_segment(int,mapping = aes(x=V1,xend=V2,y=a/2,yend=a/2),size=1)
    png(filename=str_c(output_directory, "/", intervalname, ".png"),width = 700,height = 500)
    print(p)
    dev.off()
  }else{
    png(filename=str_c(output_directory, "/", intervalname, ".png"))
    print(p)
    dev.off()
  }
  print(paste0(intervalname,"plot completed"))
}

### definition of interval union
union_intv<-function(mclustResult,intervalmods,lidu,criteria,overlap){
  validVector <- c()
  intv<-Intervals()
  if(mclustResult$G==1){
    intv<-Intervals(rbind(intv,range(intervalmods[mclustResult$classification == 1]$DeltaMass)))
  }else{
  for (i in 1:mclustResult$G){
    validVector <- c(validVector, 
                     nrow(intervalmods[mclustResult$classification == names(mclustResult$parameters$mean[i])]) > 0)            ##########groups which nrow!=0########
  }
  
  if(mclustResult$modelName=="E"){
    meanVector <- mclustResult$parameters$mean[validVector]
    sdVector<-rep(sqrt(mclustResult$parameters$variance$sigmasq),length(meanVector))
    }else{
      abnormal<-quantile(sqrt(mclustResult$parameters$variance$sigmasq))[3]+1.5*IQR(sqrt(mclustResult$parameters$variance$sigmasq))
      validVector <- validVector & sqrt(mclustResult$parameters$variance$sigmasq) < criteria & sqrt(mclustResult$parameters$variance$sigmasq) <= abnormal                                        
      meanVector <- mclustResult$parameters$mean[validVector]
      sdVector <- sqrt(mclustResult$parameters$variance$sigmasq)[validVector]
    }
  order<-order(meanVector)
  if (length(meanVector) > 0){
    if(length(meanVector) == 1){
      intv<-Intervals(rbind(intv,range(intervalmods[mclustResult$classification == names(meanVector[1])]$DeltaMass)))
      }else{
    intv<-Intervals()
    z<-1
    miu1<-meanVector[order[z]]
    sigma1<-sdVector[order[z]]
    wtf1<-names(meanVector[order[z]])
    a1<-nrow(intervalmods[mclustResult$classification == wtf1])
    ordervector<-c(intervalmods[mclustResult$classification == wtf1]$DeltaMass)
    while (z<length(meanVector)){
      miu2<-meanVector[order[z+1]]
      sigma2<-sdVector[order[z+1]]
      wtf2<-names(meanVector[order[z+1]])
      a2<-nrow(intervalmods[mclustResult$classification == wtf2])         
      p1<-a1/(a1+a2)
      p2<-a2/(a1+a2)
      if(miu2>miu1){
        OLR.r<-OLR(miu1,miu2,sigma1,sigma2,p1,p2,lidu)
      }else{
        OLR.r<-OLR(miu2,miu1,sigma2,sigma1,p2,p1,lidu)
      }
      if(OLR.r>overlap){
        miu.1<-(miu1*a1+miu2*a2)/(a1+a2)
        sigma1<-sqrt((miu1^2*a1+miu2^2*a2+(a1-1)*sigma1^2+(a2-1)*sigma2^2-((a1+a2)*miu.1^2))/(a1+a2-1))
        a1<-a1+a2
        miu1<-miu.1
        ordervector<-c(ordervector,intervalmods[mclustResult$classification == wtf2]$DeltaMass)
        z<-z+1
      }else{
        intv<-Intervals(rbind(intv,range(ordervector)))
        ordervector<-c(intervalmods[mclustResult$classification == wtf2]$DeltaMass) 
        miu1<-meanVector[order[z+1]]
        sigma1<-sdVector[order[z+1]]
        a1<-a2
        z<-z+1 
      }
    }
    intv<-Intervals(rbind(intv,range(ordervector)))
      }
  }
  }
  intvUnion<-intv
  return(intvUnion)
}

# clustering method
print("function1 completed")
DeltaMass.Func <- function(DM) {

  rValue <- data.table()
  
  lowerLimit <- DM - 0.5
  upperLimit <- DM + 0.5
  lidu<-0.0001
  criteria<-0.01483
  overlap<-0.8
  # name every interval like "-200.5_-199.5"
  intervalname <- paste(lowerLimit, "_", upperLimit)
  #print(paste("data interval:", intervalname))
  Log(paste("Data segment:", intervalname))
  
  # get interval data
  intervalmods <- massdata[DeltaMass %between% c(lowerLimit, upperLimit)]     ##########lowerLimit<=DeltaMass<=upperLimit##########
  
  #print(paste(intervalname, "segment contains", nrow(intervalmods), "datasets"))
  Log(paste(intervalname, "segment contains", nrow(intervalmods), "datasets"))
  
  # numbers of dataset in each segment should greater than 10, drop those less than 10, because mclust does not behave well with few input
  if (nrow(intervalmods) > 10)
  {  
    # process each segment with mclust
    print(str_c(intervalname, " ", nrow(intervalmods)))
    mclustResult <- Mclust(intervalmods$DeltaMass)
    print("Mclust")
    
    # save clustering parameter
    intvUnion<-union_intv(mclustResult,intervalmods,lidu,criteria,overlap)
    variancecount<-nrow(intvUnion)
    # plot histgram
     plot.clust(intervalmods,mclustResult,binwidth = lidu,nd=1/lidu,output_directory,intervalname,intvUnion,"") 
      
      #iteration of nls
      print(paste0("cluster number =",variancecount))
      Log(paste0("cluster number =",variancecount))
      if(!is.null(variancecount)){
      if (variancecount > 0) {
        for (j in 1:variancecount)
        {
          intvStartValue <- as.numeric(intvUnion[j,1])
          intvStopValue <- as.numeric(intvUnion[j,2])
          intvMods <- massdata[DeltaMass %between% c(intvStartValue, intvStopValue)][order(DeltaMass)]
          
          if (nrow(intvMods) > 10)
          {
            # normalize data
            print(str_c("breaks:", (max(intvMods$DeltaMass) - min(intvMods$DeltaMass))/lidu))
            h <-  hist(intvMods$DeltaMass,freq = TRUE,breaks = round((max(intvMods$DeltaMass) - min(intvMods$DeltaMass))/lidu))
            nlsResult<-nlsfunc(intvMods,h)
            
            if (!(is.null(nlsResult) | is.numeric(nlsResult))) {
              if (nlsResult$convInfo$isConv)
              {
                nlsResult_summary <- summary(nlsResult)
                meanValue <- round(as.numeric(nlsResult_summary$coefficients[1]), 5)
                sdValue <- as.numeric(nlsResult_summary$coefficients[2])
                countValue <- nrow(intvMods)
                x <- h$mids
                y <- dnorm(x, meanValue, sdValue)
                rSquared <-cor(h$density,predict(nlsResult))
                
                sequenceNumber <- str_c("Mean", meanValue, 
                                        "_SD", round(sdValue, 5), 
                                        "_Count", countValue,
                                        "_RS", round(rSquared, 2))
                
                # save the result
                rValue <- rbind(rValue, data.table(mean = meanValue,
                                                   variance = sdValue,
                                                   frenquency = countValue, 
                                                   RSquared = rSquared,
                                                   SequenceNumber = sequenceNumber))
                
                png(str_c(output_directory, "/", sequenceNumber, ".png"),width = 700,height = 500,bg="transparent")
                p1<-ggplot(intvMods,mapping = aes(x=DeltaMass))+ylab("PSM count")+geom_histogram(aes(y=..count..),binwidth = lidu,stat="bin",fill="steelblue",alpha=0.9)
                a<-ggplot_build(p1)$data[[1]]$count%>%max()
                p1<-p1+stat_density(aes(x=DeltaMass,y=after_stat(scaled*a)),colour="black",fill="transparent",n=1/lidu)+geom_line(aes(x=DeltaMass,y = dnorm(DeltaMass, meanValue, sdValue)*mean(h$counts/h$density),na.rm=TRUE), color = 'firebrick')+ggtitle(paste0("Mean=",meanValue," SD=",round(sdValue,5)," R-Square=",round(rSquared*100,3),"%"))+ggtheme+theme(title = element_text(family="serif",size=14))
                print(p1)
                dev.off()
                
                Log(paste("mean =",meanValue, 
                          "variance =", sdValue, 
                          "frenquence =", countValue))
                
                massdata[DeltaMass %between% c(intvStartValue, intvStopValue), Cluster_Index := sequenceNumber]
                write.table(intvMods, str_c(output_directory, "/", sequenceNumber, '.csv'), sep = ',', row.names = FALSE)
                
                #creat output files
                output_sub<-str_c(output_directory,"/",sequenceNumber)
                dir.create(output_sub)
                
                #amino acid cluster
                group_aa <- c("D|E","I|L","S|T","Y","H|K","N|Q","G","A","C","R","W","V","M","P")
                color_aa<-mclust.options("classPlotColors")[1:14]
                seq_1<-intvMods$Sequence
                seq_1<-str_remove(seq_1,"^[A-Z]{1}")
                seq_1<-str_remove(seq_1,"[A-Z]{1}$")
                aminoname<-str_replace(group_aa,"\\|"," or ")
                seq_s<-str_split(seq_1,"")%>% unlist() %>% table() %>% as.data.frame()
                colnames(seq_s)<-c("AA","Freq")
                seq_s1<-as.data.frame(table(str_extract(intvMods$Modifications,"[A-Z]{1}")))
                colnames(seq_s1)<-c("AA","Freq_location")
                seq_s<-merge(seq_s,seq_s1,by="AA",all.x=TRUE)
                seq_s$P_all<-seq_s$Freq/sum(seq_s$Freq)
                seq_s$P_location<-seq_s$Freq_location/seq_s$Freq
                seq_s$Rela<-seq_s$P_location-seq_s$P_all
                write.csv(seq_s,paste0(output_sub,"/AminoAcid Proportion.csv"),row.names = FALSE)
                intvMods_sub<-list()
                for (aa in 1:length(group_aa)) {
                  intvMods_sub[[aa]]<-list(total=setDT(setDF(intvMods)[grep(group_aa[aa],intvMods$Modifications),]))
                }
                subcluster_result<-NULL
                color_vector<-c()
                aa.plot<-data.frame()
                for (aa in 1:length(group_aa)) {
                  if(nrow(intvMods_sub[[aa]]$total)>10){
                  mclustResult.sub<-Mclust(intvMods_sub[[aa]]$total$DeltaMass)
                  intvUnion.sub<-union_intv(mclustResult.sub,intvMods_sub[[aa]]$total,lidu,criteria,overlap)
                  nd.sub<-round((range(intvMods_sub[[aa]]$total$DeltaMass)[2]-range(intvMods_sub[[aa]]$total$DeltaMass)[1])/lidu,0)
                  
                  plot.clust(intvMods_sub[[aa]]$total,mclustResult.sub,lidu,nd.sub,output_sub,aminoname[aa],intvUnion.sub,group_aa[aa])
                  intvcount<-nrow(intvUnion.sub)
                    if (!is.null(intvcount)) {
                      if(intvcount > 0 ){
                      for (int in 1:nrow(intvUnion.sub)) 
                      {
                        intvStartValue_sub <- as.numeric(intvUnion.sub[int,1])
                        intvStopValue_sub <- as.numeric(intvUnion.sub[int,2])
                        intvMods_sub[[aa]][[int+1]]<-intvMods_sub[[aa]]$total[DeltaMass %between% c(intvStartValue_sub, intvStopValue_sub)][order(DeltaMass)]
                       
                        if (nrow(intvMods_sub[[aa]][[int+1]]) > 10)
                        {
                          # normalize data
                          print(str_c(group_aa[aa],":",int,"breaks:", (max(intvMods_sub[[aa]][[int+1]]$DeltaMass) - min(intvMods_sub[[aa]][[int+1]]$DeltaMass))/lidu))
                          h1 <-  hist(intvMods_sub[[aa]][[int+1]]$DeltaMass,freq = TRUE,breaks = round((max(intvMods_sub[[aa]][[int+1]]$DeltaMass) - min(intvMods_sub[[aa]][[int+1]]$DeltaMass))/lidu)+1)
                          nlsResult.sub<-nlsfunc(intvMods_sub[[aa]][[int+1]],h1)
                          if (!(is.null(nlsResult.sub) | is.numeric(nlsResult.sub))) {
                            if (nlsResult.sub$convInfo$isConv)
                            {
                              color_vector<-c(color_vector,color_aa[aa])
                              nlsResult_summary.sub <- summary(nlsResult.sub)
                              meanValue.sub <- round(as.numeric(nlsResult_summary.sub$coefficients[1]), 5)
                              sdValue.sub <- as.numeric(nlsResult_summary.sub$coefficients[2])
                              countValue.sub <- nrow(intvMods_sub[[aa]][[int+1]])
                              x.sub <- h1$mids
                              y.sub <- dnorm(x.sub, meanValue.sub, sdValue.sub)
                              rSquared.sub <- cor(h1$density,predict(nlsResult.sub))
                              vector_output<-c(meanValue.sub,sdValue.sub,countValue.sub,rSquared.sub,aminoname[aa],int,intvUnion.sub[int,1],intvUnion.sub[int,2])
                              subcluster_result<-rbind(subcluster_result,vector_output)
                              sub<-data.frame(DeltaMass=seq(from=meanValue.sub-3*sdValue.sub,to=meanValue.sub+3*sdValue.sub,by=0.00001),mean=meanValue.sub,sd=sdValue.sub,bei=sum(h1$counts)/100000,AminoAcid=group_aa[aa],Color=paste0(group_aa[aa],"_",int))
                              sub$y<-dnorm(sub$DeltaMass,sub$mean,sub$sd)*(sum(h1$counts)/sum(h1$density))
                              aa.plot<-rbind(aa.plot,sub)
                              }
                           }else{Log(paste0("sub",int,"[",intvStartValue_sub,",",intvStopValue_sub,"]","third times nls fails."))}
                          }
                         }
                        }
                    }
                  }
                }
                if(nrow(aa.plot)>0){
                png(filename=str_c(output_sub, "/AminoAcid_Group_all.png"),width = 700,height = 500)
                aa.plot$Color<-factor(aa.plot$Color,levels=unique(aa.plot$Color))
                pp<-ggplot(aa.plot)+geom_point(aes(x=DeltaMass,y = y,color=Color),show.legend = FALSE,size=.2)+geom_line(aes(x=min(DeltaMass),y=0,color=Color),size=5)+scale_color_manual(values = color_vector)+ggtheme+ylab("PSM Count")+ guides(guide_legend(title = "AminoAcid Group"))
                colnames(subcluster_result)<-c("Mean","Sd","Count","R-Square","AA-Group","sub-AA","Lower","Upper")
                print(pp)
                dev.off()
                write.csv(subcluster_result,file = paste0(output_sub,"/","subcluster.csv"),row.names = FALSE)}
                     } 
                    }else {
                Log("third times nls fails.")
              }
            }
          }
  }
      }
  return(rValue)
      }
  }
Start<-Sys.time()
sepplementalTable3<-data.table()
sepplementalTable3 <- do.call("rbind", apply(as.array(276:500), 1, DeltaMass.Func))

End<-Sys.time()
