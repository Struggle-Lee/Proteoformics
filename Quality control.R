Sepsis_a11<-fread("Z:/ICU/-a1/pFindTask5/result/pFind-Filtered.spectra",fill=TRUE)
Sepsis_a21<-fread("Z:/ICU/-a2/A2/result/pFind-Filtered.spectra",fill=TRUE)

################################################
###################结果重聚类###################
################################################
allmodification<-unique(unlist(str_split(c(Sepsis_a11$Modification,Sepsis_a21$Modification),";")))  
allmodification<-unique(str_remove(allmodification,"^[0-9]{1,},"))
#allmodification<-unique(str_remove_all(allmodification,"(?<=\\[)[a-zA-Z]{3,7}[N|C]{1}-term[A-Z]{0,1}(?=\\])|(?<=\\[)[A-Z]{1}(?=\\])"))
#allmodification<-unique(str_remove_all(allmodification,"\\[[a-zA-Z]{3,7}[N|C]{1}-term[A-Z]{0,1}\\]$|\\[[A-Z]{1}\\]$"))

allmodification<-unique(str_remove_all(allmodification[-1],"\\(.{2,}\\)"))  
allmodification<-data.frame(Mod=allmodification,Count=0)
allmodification$Mod<-str_remove(allmodification$Mod,"\\]$")
allmodification$Mod<-str_remove(allmodification$Mod,"\\[$")
allmodification$Mod.pattern<-str_replace_all(allmodification$Mod,"\\[","\\\\[")
allmodification$Mod.pattern<-str_replace_all(allmodification$Mod.pattern,"\\]","\\\\]")
allmodification$Mod.pattern<-str_replace_all(allmodification$Mod.pattern,"\\(","\\\\(")
allmodification$Mod.pattern<-str_replace_all(allmodification$Mod.pattern,"\\)","\\\\)")
allmodification$Mod.pattern<-str_replace_all(allmodification$Mod.pattern,"\\+","\\\\+")
nlsfunc2<-function(intvMods,h){
  f <- function(x, meanValue, sdValue){
    return ((1 / (sqrt(2 * pi) * sdValue)) * exp(-((x-meanValue)*(x-meanValue)) / (2*sdValue*sdValue)))}
  nlsResult <- tryCatch(
    {
      # parameter estimate
      nls(y~f(x,m,s),
          data = data.frame(x = h$mids,y = h$density), 
          control = list(maxiter = 50, tol = 1e-05, minFactor = 1/1024,
                         printEval = FALSE, warnOnly = FALSE), 
          start = list(m = mean(intvMods$ppm), 
                       s = sd(intvMods$ppm)))
    },
    error = function(err) {
      print("first nls error:")
      print(err)
    },
    warning = function(war) {
      print("first nls warning:")
      print(war)
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
            control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024,
                           printEval = FALSE, warnOnly = FALSE), 
            start = list(m = mean(intvMods$ppm), 
                         s = sd(intvMods$ppm)/2))
      },
      error = function(err) {
        print("second nls error:")
        print(err)
      },
      warning = function(war) {
        print("second nls warning:")
        print(war)
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
              control = list(maxiter = 200, tol = 1e-05, minFactor = 1/1024,
                             printEval = FALSE, warnOnly = FALSE), 
              start = list(m = mean(intvMods$ppm), 
                           s = sd(intvMods$ppm)/4))
        },
        error = function(err) {
          print("third nls error:")
          print(err)
        },
        warning = function(war) {
          print("third nls warning:")
          print(war)
        },
        finally = {
        })    
    }else {
      print("second NLS success.")
    }
  }else {
    print("first NLS success.")
  }
  return(nlsResult)
}

output_directory<-"D:/Clust0709/pFind_sepsis_deltamass"
for (i in 1:nrow(allmodification)) {
  allmodification$Count[i]<-length(which(str_detect(Sepsis_a11$Modification,paste0(",",allmodification$Mod.pattern[i]))))+length(which(str_detect(Sepsis_a21$Modification,paste0(",",allmodification$Mod.pattern[i]))))
  print(i)
}
rejudge.mass<-rejudge
rejudge<-data.frame()
for (i in 1:nrow(allmodification)) {
  allmodification$Count[i]<-length(which(str_detect(Sepsis_a11$Modification,paste0(",",allmodification$Mod.pattern[i]))))+length(which(str_detect(Sepsis_a21$Modification,paste0(",",allmodification$Mod.pattern[i]))))
  sub<-rbind(Sepsis_a11[str_detect(Sepsis_a11$Modification,paste0(",",allmodification$Mod.pattern[i])),],Sepsis_a21[str_detect(Sepsis_a21$Modification,paste0(",",allmodification$Mod.pattern[i])),])
  if(nrow(sub)>1){
    sub$ppm<-(sub$`Mass_Shift(Exp.-Calc.)`/sub$`Calc.MH+`)*1000000
    h <-  hist(sub$ppm,freq = TRUE,breaks = round((max(sub$ppm) - min(sub$ppm))/0.05))
    nlsResult<-nlsfunc2(sub,h)                      
    if (!(is.null(nlsResult) | is.numeric(nlsResult))) {
      if (nlsResult$convInfo$isConv)
      {
        nlsResult_summary <- summary(nlsResult)
        name<-str_replace_all(allmodification$Mod[i],"\\-\\>","\\_")
        name<-str_replace_all(name,"\\/","\\_")
        meanValue <- round(as.numeric(nlsResult_summary$coefficients[1]), 5)
        sdValue <- as.numeric(nlsResult_summary$coefficients[2])
        countValue <- nrow(sub)
        x <- h$mids
        y <- dnorm(x, meanValue, sdValue)
        rSquared <-cor(h$density,predict(nlsResult))
        sequenceNumber <- str_c("Name",name,
                                "Mean", meanValue, 
                                "_SD", round(sdValue, 5), 
                                "_Count", countValue,
                                "_RS", round(rSquared, 2))
        rvalue<-c(Name=name,Mean=meanValue,SD=round(sdValue, 5),Count=countValue,RS=round(rSquared, 2),lillie=lillie.test(sub$`Mass_Shift(Exp.-Calc.)`)$p.value)
        rejudge<-rbind(rejudge,rvalue)
        png(str_c(output_directory, "/", str_replace(sequenceNumber,"\\-\\>","\\_"), ".png"),width = 700,height = 500,bg="transparent")
        p1<-ggplot(sub,mapping = aes(x=ppm))+ylab("PSM count")+geom_histogram(aes(y=..count..),binwidth = 0.05,stat="bin",fill="steelblue",alpha=0.9)
        a<-ggplot_build(p1)$data[[1]]$count%>%max()
        ind<-sum(h$counts)/sum(h$density)
        llline<-data.frame(DeltaMass=seq(min(sub$ppm),max(sub$ppm),0.05),Ratio=dnorm(seq(min(sub$ppm),max(sub$ppm),0.05), meanValue, sdValue)*ind)
        p1<-p1+stat_density(aes(x=ppm,y=after_stat(scaled*a)),colour="black",fill="transparent",n=10000)+geom_line(data=llline,aes(x=DeltaMass,y =Ratio), color = 'firebrick',size=1.5)+ggtitle(paste0(name," Mean=",meanValue," SD=",round(sdValue,5)," R-Square=",round(rSquared*100,3),"%"))+geom_vline(xintercept = mean(sub$ppm),size=1.5,fill="red",colour="red")+ggtheme+theme(title = element_text(family="serif",size=14))
        print(p1)
        dev.off()
      }
    }
    write.csv(sub,paste0("D:/Clust0709/pFind_sepsis_deltamass/",str_replace(sequenceNumber,"\\-\\>","\\_"),".csv"))
  }
  print(i)
}
colnames(rejudge)<-c("Name","Mean","SD","Count","RSquare","P")  
rejudge<-list.files("D:/Clust0709/pFind_sepsis_deltamass/11/","\\.png")  
rejudge<-str_remove(rejudge,"\\.png")  
rejudge<-data.frame(Name=str_extract(rejudge,"(?<=Name).*(?=Mean)"),
                    Mean=as.numeric(str_extract(rejudge,"(?<=Mean).*(?=_SD)")),
                    SD=as.numeric(str_extract(rejudge,"(?<=_SD).*(?=_Count)")),
                    Count=as.numeric(str_extract(rejudge,"(?<=_Count).*(?=_RS)")),
                    RS=as.numeric(str_extract(rejudge,"(?<=_RS).*")))  

rejudge2<-list.files("D:/Clust0709/pFind_sepsis_deltamass/","\\.png")  
rejudge2<-str_remove(rejudge2,"\\.png")  
rejudge2<-data.frame(Name=str_extract(rejudge2,"(?<=Name).*(?=Mean)"),
                     Mean=as.numeric(str_extract(rejudge2,"(?<=Mean).*(?=_SD)")),
                     SD=as.numeric(str_extract(rejudge2,"(?<=_SD).*(?=_Count)")),
                     Count=as.numeric(str_extract(rejudge2,"(?<=_Count).*(?=_RS)")),
                     RS=as.numeric(str_extract(rejudge2,"(?<=_RS).*")))

ggplot(rejudge,mapping = aes(RS,Count))+geom_point()+scale_y_continuous(trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+ggtheme


#########################################
#################数据质控################
#########################################
########以pFind-filtered.spectra作为输入，总结各文件信息，生成主要的质控信息
Sepsis_all<-rbind(Sepsis_a11,Sepsis_a21)
Sepsis_all<-setDF(Sepsis_all)
Sepsis_summary<-data.frame()
Sepsis_all$File<-str_extract(Sepsis_all$File_Name,".*_[a-zA-Z0-9]{1,}(?=\\.)")
File<-unique(Sepsis_all$File)
originFile<-list.files("Z:/ICU/-A","\\.ms2")
for (i in 1:length(unique(Sepsis_all$File))) {
  Sepsis_sub<-Sepsis_all[which(Sepsis_all$File==File[i]),]
  ms2<-read_lines(paste0("Z:/ICU/-A/",originFile[which(str_detect(originFile,File[i]))]))
  row1<-c(Total=nrow(Sepsis_sub),
          ToalMS2=length(which(str_detect(ms2,"^Z\\\t"))),
          Unmod=length(which(Sepsis_sub$Modification=="")),
          Mod=length(which(Sepsis_sub$Modification!="")),
          Sequence=length(unique(Sepsis_sub$Sequence)),
          Peptide=length(unique(paste0(Sepsis_sub$Sequence,Sepsis_sub$Modification))),
          Precursor.Mean=mean(Sepsis_sub$`Mass_Shift(Exp.-Calc.)`/Sepsis_sub$`Calc.MH+`)*1000000,
          Precursor.SD=sd((Sepsis_sub$`Mass_Shift(Exp.-Calc.)`/Sepsis_sub$`Calc.MH+`)*1000000),
          Fragment.Mean=mean(Sepsis_sub$Avg.Frag.Mass.Shift),
          Fragment.SD=sd(Sepsis_sub$Avg.Frag.Mass.Shift),
          Precursor.Mean.Unmod=mean(Sepsis_sub$`Mass_Shift(Exp.-Calc.)`[which(Sepsis_sub$Modification=="")]/Sepsis_sub$`Calc.MH+`[which(Sepsis_sub$Modification=="")])*1000000,
          Precursor.SD.Unmod=sd((Sepsis_sub$`Mass_Shift(Exp.-Calc.)`[which(Sepsis_sub$Modification=="")]/Sepsis_sub$`Calc.MH+`[which(Sepsis_sub$Modification=="")])*1000000),
          Fragment.Mean.Unmod=mean(Sepsis_sub$Avg.Frag.Mass.Shift[which(Sepsis_sub$Modification=="")]),
          Fragment.SD.Unmod=sd(Sepsis_sub$Avg.Frag.Mass.Shift[which(Sepsis_sub$Modification!="")]),
          Precursor.Mean.Mod=mean(Sepsis_sub$`Mass_Shift(Exp.-Calc.)`[which(Sepsis_sub$Modification!="")]/Sepsis_sub$`Calc.MH+`[which(Sepsis_sub$Modification!="")])*1000000,
          Precursor.SD.Mod=sd((Sepsis_sub$`Mass_Shift(Exp.-Calc.)`[which(Sepsis_sub$Modification!="")]/Sepsis_sub$`Calc.MH+`[which(Sepsis_sub$Modification!="")])*1000000),
          Fragment.Mean.Mod=mean(Sepsis_sub$Avg.Frag.Mass.Shift[which(Sepsis_sub$Modification!="")]),
          Fragment.SD.Mod=sd(Sepsis_sub$Avg.Frag.Mass.Shift[which(Sepsis_sub$Modification!="")]),
          Spectra.0=length(which(str_detect(Sepsis_sub$File_Name,"\\.0\\.dta"))),
          Enzymatic=paste0(paste(c("Non-Specific","C-term Specific","N-term Specific","Specific"),as.vector(table(Sepsis_sub$Specificity)),sep="\t"),collapse = "\t"),
          Charge=paste0(paste(names(table(Sepsis_sub$Charge)),as.vector(table(Sepsis_sub$Charge)),sep="\t"),collapse = "\t"),
          Misscleavage=paste0(paste(names(table(Sepsis_sub$Miss.Clv.Sites)),as.vector(table(Sepsis_sub$Miss.Clv.Sites)),sep="\t"),collapse = "\t")
  )
  Sepsis_summary<-rbind(Sepsis_summary,row1)
  print(i)
}
remove(ms2)
colnames(Sepsis_summary)<-names(row1)
colnames(Sepsis_summary)[20:22]<-c("Enzymatic","Charge","Misscleavage")
for (i in 1:19) {
  Sepsis_summary[,i]<-as.numeric(Sepsis_summary[,i])
}
Sepsis_summary<-cbind(File=File,Sepsis_summary)

####
library(ggplot2)
library(scales)
library(ggthemes)

#mass error
masserror2<-data.frame(Mod=Sepsis_all$Modification,Error=(Sepsis_all$`Mass_Shift(Exp.-Calc.)`/Sepsis_all$`Calc.MH+`)*1000000)
masserror2$Mod[which(masserror2$Mod!="")]<-"Mod"
masserror2$Mod[which(masserror2$Mod=="")]<-"UnMod"
ggplot(data = masserror2,mapping = aes(x=Error,fill=Mod))+geom_density(alpha=0.7,adjust=.5,colour="transparent")+ggtheme+scale_fill_manual(values = c(cols_first[2],cols_first[3]))+ylab("Density")+xlab("Mass Error(PPM)")                       
cairo_pdf(file="D:/Clust0709/pFind_sepsis/MassError.pdf",width = 7,height = 5,bg="transparent",family="Serif") 
print(ggplot(data = masserror2,mapping = aes(x=Error,fill=Mod))+geom_density(alpha=0.7,adjust=.5,colour="transparent")+ggtheme+scale_fill_manual(values = c(cols_first[2],cols_first[3]))+ylab("Density")+xlab("Mass Error(PPM)"))                      
dev.off()

#Target Distribution                       
TDA<-rbind(Sepsis_a11[,c(9,11)],Sepsis_a21[,c(9,11)])                       
TDA$Mod[which(TDA$Modification!="")]<-"Mod"
TDA$Mod[which(TDA$Modification=="")]<-"UnMod"          
ggplot(data = TDA,mapping = aes(x=Raw_Score,fill=Mod))+geom_density()+ggtheme                       

TDA<-rbind(Sepsis_all[,c(9,11,16)])    
TDA$Mod[which(TDA$Modification!="")]<-"Mod"
TDA$Mod[which(TDA$Modification=="")]<-"UnMod"
TDA$Mod[which(TDA$`Target/Decoy`=="decoy")]<-"Decoy"
ggplot(data = TDA,mapping = aes(x=Raw_Score,fill=Mod))+geom_density(alpha=0.7,adjust=.5,colour="transparent")+ggtheme+scale_fill_manual(values = c(cols_first[2],cols_first[3]))+ylab("Density")
cairo_pdf(file="D:/Clust0709/pFind_sepsis/Raw_score.pdf",width = 7,height = 5,bg="transparent",family="Serif") 
print(ggplot(data = TDA,mapping = aes(x=Raw_Score,fill=Mod))+geom_density(alpha=0.7,adjust=.5,colour="transparent")+ggtheme+scale_fill_manual(values = c(cols_first[2],cols_first[3]))+ylab("Density"))                      
dev.off()                       

####Summary                        
#remove "SEP_00202_42L_a" "SEP_00176_42L_a" "SEP_00168_42L_a" "SEP_00174_42L_a" "HYC_00046_42L_a" "SEP_00226_42L_a" "SEP_00172_42L_a" "HYC_00049_42L_a" "SEP_00166_42L_a" "SEP_00100_40L_a" "SEP_00169_42L_a" "SEP_00175_42L_a" "SEP_00173_42L_a" "SEP_00167_42L_a" "HYC_00047_42L_a" "SEP_00177_42L_a" "SEP_00191_42L_a" "SEP_00240_42L_a" "SEP_00180_42L_a" "HYC_00045_42L_a" "SEP_00301_42L_a" "SEP_00131_40L_a"
Sepsis_summary<-Sepsis_summary[-which((Sepsis_summary$Total/Sepsis_summary$ToalMS2)<0.2),]
Sepsis_summary.sub<-Sepsis_summary
for (i in 1:nrow(Sepsis_summary.sub)) {
  Sepsis_sub<-Sepsis_all[which(Sepsis_all$File==Sepsis_summary$File[i]),]
  Sepsis_sub$Info<-paste0(Sepsis_sub$Sequence,Sepsis_sub$Modification)
  Sepsis_sub2<-data.frame(Info=unique(Sepsis_sub$Info),Enzymatic=Sepsis_sub$Specific[match(unique(Sepsis_sub$Info),Sepsis_sub$Info)],Charge=Sepsis_sub$Charge[match(unique(Sepsis_sub$Info),Sepsis_sub$Info,)],Misscleavage=Sepsis_sub$Miss.Clv.Sites[match(unique(Sepsis_sub$Info),Sepsis_sub$Info,)])
  Sepsis_summary.sub$Misscleavage_seq[i]<-paste0(paste(names(table(Sepsis_sub2$Misscleavage)),as.vector(table(Sepsis_sub2$Misscleavage)),sep="\t"),collapse = "\t")
  Sepsis_summary.sub$Enzymatic_seq[i]<-paste0(paste(names(table(Sepsis_sub2$Enzymatic)),as.vector(table(Sepsis_sub2$Enzymatic)),sep="\t"),collapse = "\t")
  Sepsis_summary.sub$Charge_seq[i]<-paste0(paste(names(table(Sepsis_sub2$Charge)),as.vector(table(Sepsis_sub2$Charge)),sep="\t"),collapse = "\t")
  print(i)
}

arrang<-function(x){
  a<-unlist(str_split(x,"\\\t"))
  tab<-as.data.frame(cbind(a[seq(1,(length(a)-1),2)],a[seq(2,length(a),2)]))
  tab$ratio<-(as.numeric(tab[,2])/sum(as.numeric(tab[,2])))*100
  return(tab)
}

##Enzymatic
Enzymatic<-data.frame()
#Enzymatic<-as.data.frame(arrang(Sepsis_summary.sub$Enzymatic_seq[1]))
#colnames(Enzymatic)<-c("Enzymatic","Spectra Count","Ratio")
#Enzymatic<-cbind(File=Sepsis_summary.sub$File[1])
for (i in 1:nrow(Sepsis_summary.sub2)) {
  sub<-data.frame(arrang(Sepsis_summary.sub2$Enzymatic_seq[i]))
  sub<-cbind(File=Sepsis_summary.sub2$File[i],sub)
  colnames(sub)<-c("File","Enzymatic","Spectra Count","Ratio")
  Enzymatic<-rbind(Enzymatic,sub)
  print(i)}
Enzymatic$Enzymatic<-str_replace(Enzymatic$Enzymatic,"0","Non Specific")
Enzymatic$Enzymatic<-str_replace(Enzymatic$Enzymatic,"1","C-term Specific")
Enzymatic$Enzymatic<-str_replace(Enzymatic$Enzymatic,"2","N-term Specific")
Enzymatic$Enzymatic<-str_replace(Enzymatic$Enzymatic,"3","Specific")
P.Enzymatic<-ggplot(Enzymatic,mapping = aes(x=File,y=Ratio,fill=Enzymatic))+
  geom_bar(stat="identity",position = "stack",width=.9)+
  theme_light()+theme(legend.title = element_text(size=14,face = "bold",family = "serif"),
                      axis.text.y = element_text(colour="black", size=14,face = "bold",family = "serif"),
                      axis.text.x = element_blank(),axis.title.x = element_blank(),
                      axis.title = element_text(colour="black", size=14,face = "bold",family = "serif"),
                      legend.position = "right",title =element_text(size=14,face = "bold",family = "serif"),
                      legend.text = element_text(size=14,face = "bold",family = "serif"))+
  scale_fill_manual(values = c("#F28E2B","#FFBE7D","#A0CBE8","#8CD17D"))+
  ylab("PSM Count (%)")+scale_y_continuous(expand = c(0,0))
P.Enzymatic%>%
  insert_bottom(P.anno,height = 0.05) 
cairo_pdf(file="D:/Clust0709/pFind_sepsis/Enzymatic.pdf",width = 7,height = 5,bg="transparent",family="Serif")
print(P.Enzymatic%>%
        insert_bottom(P.anno,height = 0.05))
dev.off()

##Misscleavge
Misscleavage<-data.frame()
for (i in 1:nrow(Sepsis_summary.sub2)) {
  sub<-data.frame(arrang(Sepsis_summary.sub2$Misscleavage_seq[i]))
  sub<-cbind(File=Sepsis_summary.sub2$File[i],sub)
  colnames(sub)<-c("File","Misscleavage","Spectra Count","Ratio")
  Misscleavage<-rbind(Misscleavage,sub)
  print(i)}
P.Misscleavage<-ggplot(Misscleavage,mapping = aes(x=File,y=Ratio,fill=factor(Misscleavage,levels = c(8,7,6,5,4,3,2,1,0))))+
  geom_bar(stat="identity",position = "stack",width=.9)+
  theme_light()+theme(legend.title = element_blank(),
                      axis.text.y = element_text(colour="black", size=14,face = "bold",family = "serif"),
                      axis.text.x = element_blank(),axis.title.x = element_blank(),
                      axis.title = element_text(colour="black", size=14,face = "bold",family = "serif"),
                      legend.position = "right",title =element_text(size=14,face = "bold",family = "serif"),
                      legend.text = element_text(size=14,face = "bold",family = "serif"))+
  scale_fill_manual(values = c("gray10","gray30","gray50","gray70","gray90",rev(cols_first[1:4])))+
  ylab("PSM Count (%)")+scale_y_continuous(expand = c(0,0))
P.Misscleavage%>%
  insert_bottom(P.anno,height = 0.05)                    
cairo_pdf(file="D:/Clust0709/pFind_sepsis/Misscleavage.pdf",width = 7,height = 5,bg="transparent",family="Serif")
print(P.Misscleavage%>%
        insert_bottom(P.anno,height = 0.05))
dev.off()

##Charge
Charge<-data.frame()
for (i in 1:nrow(Sepsis_summary.sub2)) {
  sub<-data.frame(arrang(Sepsis_summary.sub2$Charge_seq[i]))
  sub<-cbind(File=Sepsis_summary.sub2$File[i],sub)
  colnames(sub)<-c("File","Charge","Spectra Count","Ratio")
  Charge<-rbind(Charge,sub)
  print(i)}
P.Charge<-ggplot(Charge,mapping = aes(x=File,y=Ratio,fill=factor(Charge,levels = c(7,6,5,4,3,2))))+
  geom_bar(stat="identity",position = "stack",width=.9)+
  theme_light()+theme(legend.title = element_blank(),
                      axis.text.y = element_text(colour="black", size=14,face = "bold",family = "serif"),
                      axis.text.x = element_blank(),axis.title.x = element_blank(),
                      axis.title = element_text(colour="black", size=14,face = "bold",family = "serif"),
                      legend.position = "right",title =element_text(size=14,face = "bold",family = "serif"),
                      legend.text = element_text(size=14,face = "bold",family = "serif"))+
  scale_fill_manual(values = c("gray10","gray30","gray50","gray70",cols_first[3],cols_first[2]))+
  ylab("PSM Count (%)")+scale_y_continuous(expand = c(0,0))
P.Charge%>%
  insert_bottom(P.anno,height = 0.05)                    
cairo_pdf(file="D:/Clust0709/pFind_sepsis/Charge.pdf",width = 7,height = 5,bg="transparent",family="Serif")
print(P.Charge %>%
        insert_bottom(P.anno,height = 0.05))
dev.off()

##ModUnMod 
Spectra<-data.frame(File=Sepsis_summary.sub2$File,Mod=Sepsis_summary.sub2$Mod,UnMod=Sepsis_summary.sub2$Unmod)
Spectra<-melt(Spectra,id.vars="File")
P.ModUnMod<-ggplot(Spectra,mapping = aes(x=File,y=`PSM Count`,fill=Type))+
  geom_bar(stat="identity",position = "fill",width=.9)+
  theme_light()+
  theme(legend.title = element_text(colour="black", size=14,face = "bold",family = "serif"),
        axis.text.y = element_text(colour="black", size=14,face = "bold",family = "serif"),
        axis.text.x = element_blank(),
        axis.title = element_text(colour="black", size=14,face = "bold",family = "serif"),
        legend.position = "right",axis.title.x = element_blank(),
        title =element_text(size=14,face = "bold",family = "serif"),
        legend.text = element_text(size=14,face = "bold",family = "serif"))+
  scale_fill_manual(values = c("#FC8D62","#8DA0CB"))+scale_y_continuous(expand = c(0,0))+ylab("PSM Count (%)")
P.anno<-ggplot(Sepsis_summary.sub2,aes(x=File,y=1,fill=Group))+
  geom_tile()+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_text(size=14,face = "bold",family = "serif"),
        legend.text = element_text(size=14,face = "bold",family = "serif"))+
  scale_fill_manual(values = c(cols_first[5],cols_first[7]))
P.ModUnMod%>%
  insert_bottom(P.anno,height = 0.05) 
cairo_pdf(file="D:/Clust0709/pFind_sepsis/Mod_Unmod.pdf",width = 7,height = 5,bg="transparent",family="Serif")
print(P.ModUnMod%>%
        insert_bottom(P.anno,height = 0.05))
dev.off()


########################################################
#########从总PSM表Sepsis_all里整理修饰信息##############
########################################################
vv<-rep(1:nrow(Sepsis_all),str_count(Sepsis_all$Modification,";"))
Sepsis_Mod<-data.frame(Proteins=Sepsis_all$Proteins[vv],File=Sepsis_all$File[vv],Group=Sepsis_all$Group[vv],Mod=Mod,Seq=Sepsis_all$Sequence[vv],Positions=Sepsis_all$Positions[vv])
Sepsis_Mod$Modification<-str_remove(Sepsis_Mod$Mod,"^[0-9]{1,},")
Sepsis_Mod$ID<-apply(Sepsis_Mod,1,function(x){a<-unlist(str_split(x[1],"\\/"))
a<-a[1]
return(a)})
Sepsis_Mod$Uniprot<-str_extract(Sepsis_Mod$ID,"(?<=\\|).*(?=\\|)")
Sepsis_Mod$Name<-str_extract(Sepsis_Mod$ID,"(?<=\\|[A-Z0-9]{4,15}\\|).*(?=_HUMAN)")
Sepsis_Mod$Position_start<-as.numeric(str_extract(Sepsis_Mod$Positions,"^[0-9]{1,}"))+1
Sepsis_Mod$Position_end<-
  Sepsis_Mod$Mod_Pos<-apply(Sepsis_Mod,1,function(x){a<-as.numeric(str_extract(x[4],"[0-9]{1,}(?=,)"))
  if(a==0){a<-1
  }else if(a==(nchar(x[5])+1)){
    a<-as.numeric(a)-1
  }
  return(a)
  })
Sepsis_Mod$Pos<-Sepsis_Mod$Position_start+Sepsis_Mod$Mod_Pos-1

Sepsis_Mod$Mod_Pos[which(Sepsis_Mod$Mod_Pos==0)]<-1
Sepsis_Mod$Mod_Pos[which(Sepsis_Mod$Mod_Pos==(nchar(Sepsis_Mod$Seq)+1))]<-Sepsis_Mod$Mod_Pos[which(Sepsis_Mod$Mod_Pos==(nchar(Sepsis_Mod$Seq)+1))]-1
Sepsis_Mod$Mod1<-str_remove(Sepsis_Mod$Modification,"\\[[a-zA-Z]{3,7}[N|C]{1}-term[A-Z]{0,1}\\]$|\\[[A-Z]{1}\\]$")
Sepsis_Mod$AA<-str_extract(Sepsis_Mod$Modification,"(?<=\\[)[a-zA-Z]{3,7}[N|C]{1}-term[A-Z]{0,1}(?=\\][\\)]{0,1}$)|(?<=\\[)[A-Z]{1}(?=\\][\\)]{0,1}$)")

Sepsis_Mod$Info<-paste0(Sepsis_Mod$Mod1,"@",Sepsis_Mod$Name,"_",Sepsis_Mod$Pos,Sepsis_Mod$AA)






