

pFindSummary<-function(path,outputpath){
  #path<-"Y:/LH/LH/20220201_TNF"
  #outputpath<-"D:/Clust0709/LH_control/TNF/"
  waitread<-paste0(path,"/",list.files(path),"/pFind.summary")
  waitread<-waitread[-length(waitread)]
  all_summary<-data.frame()
  width<-700+10*length(waitread)
  for (i in 1:length(waitread)) {
    all_summary<-rbind(all_summary,file_summary(waitread[i]))
  }
  colnames(all_summary)<-c("Spectra","Scans","Peptides","Sequence","Proteins","Protein Groups","Avg_Seq_Cov","Decoy spectra","T/D Distribution","Decoy Peptides","Decoy Proteins","Decoy Protein Groups","Specific","C-term Specific","N-term Specific","Non-Specific","Modification","Charge","Misscleavage","Mixture spectra","Precursor Mean","Precursor SD","Fragment Mean","Fragment SD")
  all_summary<-cbind(File=str_extract(waitread,paste0("(?<=",path,"/).*(?=_pFindResults/pFind.summary)")),all_summary)
  all_summary$Group<-"Batch1"
  png(paste0(outputpath,"Peptide_Scans.png"),width=width,height = 450)
  p1<-ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Peptides)/as.numeric(Scans),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Peptide/Scans")+scale_y_continuous(limits = c(0,0.8),breaks = c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0))
  print(p1)
  dev.off()
  png(paste0(outputpath,"Peptide_Spectra.png"),width=width,height = 450)
  p2<-ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Peptides)/as.numeric(Spectra),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Peptide/Spectra")+scale_y_continuous(limits = c(0,0.8),breaks = c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0))
  print(p2)
  dev.off()
  png(paste0(outputpath,"Peptide_Count.png"),width=width,height = 450)
  p3<-ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Peptides),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Peptide")+scale_y_continuous(expand = c(0,0))
  print(p3)
  dev.off()
  png(paste0(outputpath,"Protein_Group.png"),width=width,height = 450)
  p4<-ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(`Protein Groups`),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Protein Group")+scale_y_continuous(expand = c(0,0))
  print(p4)
  dev.off()
  png(paste0(outputpath,"Avg_Seq_Cov.png"),width=width,height = 450)
  p5<-ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Avg_Seq_Cov),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Avg_Seq_Cov")+scale_y_continuous(expand = c(0,0))
  print(p5)
  dev.off()
  all_summary.sub<-all_summary[,c(1,14:17)]
  for (i in 2:5) {
    all_summary.sub[,i]<-as.numeric(all_summary.sub[,i])
  }
  all_summary.sub$Sum<-apply(all_summary.sub[,2:5],1,sum)
  for (i in 2:5) {
    all_summary.sub[,i]<-as.numeric(all_summary.sub[,i])/as.numeric(all_summary.sub[,6])
  }
  all_summary.sub<-all_summary.sub[,-6]
  all_summary.sub<-reshape2::melt(all_summary.sub,id.vars = "File")
  png(paste0(outputpath,"Enzymatic.png"),width=width,height = 450)
  p6<-ggplot(all_summary.sub,mapping = aes(x=factor(File,levels=unique(File)),y=value,fill=factor(variable,levels = c("Non-Specific","N-term Specific","C-term Specific","Specific"))))+geom_bar(stat="identity",position = "stack")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Enzymatic")+scale_y_continuous(expand = c(0,0))
  print(p6)
  dev.off()
  ##mod
  candidate<-c()
  for (m in 1:nrow(all_summary)) {
    candidate<-c(candidate,as.character(arrang(all_summary$Modification[m])[,1]))
  }
  candidate<-unique(candidate)
  candidate<-str_replace(candidate,"\\[","\\\\[")
  candidate<-str_replace(candidate,"\\]","\\\\]")
  all_summary<-data.frame()
  for (i in 1:length(waitread)) {
    all_summary<-rbind(all_summary,file_summary2(waitread[i]))
  }
  colnames(all_summary)<-c("Spectra","Scans","Peptides","Sequence","Proteins","Protein Groups","Avg_Seq_Cov","Decoy spectra","T/D Distribution","Decoy Peptides","Decoy Proteins","Decoy Protein Groups","Specific","C-term Specific","N-term Specific","Non-Specific","Modification","Charge","Misscleavage","Mixture spectra","Precursor Mean","Precursor SD","Fragment Mean","Fragment SD")
  all_summary<-cbind(File=str_extract(waitread,paste0("(?<=",path,"/).*(?=_pFindResults/pFind.summary)")),all_summary)
  
  mod_10<-data.frame()
  mod_sub<-data.frame()
  mod_10<-as.data.frame(arrang(all_summary$Modification[1]))
  colnames(mod_10)<-c("Modification","Spectra","Rate","Sites","Relative Ratio")
  mod_10_ratio<-mod_10[,c(1,3)]
  mod_10_occu<-mod_10[,c(1,5)]
  mod_10_spectra<-mod_10[,c(1,2)]
  for (i in 2:nrow(all_summary)) {
    mod_sub<-data.frame(arrang(all_summary$Modification[i]))
    colnames(mod_sub)<-c("Modification","Spectra","Rate","Sites","Relative Ratio")
    mod_10_ratio<-merge(mod_10_ratio,mod_sub[,c(1,3)],by.x = "Modification",by.y = "Modification",all = TRUE)
    mod_10_occu<-merge(mod_10_occu,mod_sub[,c(1,5)],by.x = "Modification",by.y = "Modification",all = TRUE)
    mod_10_spectra<-merge(mod_10_spectra,mod_sub[,c(1,2)],by.x = "Modification",by.y = "Modification",all = TRUE)
  }
  colnames(mod_10_ratio)[2:ncol(mod_10_ratio)]<-all_summary$File
  colnames(mod_10_occu)[2:ncol(mod_10_occu)]<-all_summary$File
  colnames(mod_10_spectra)[2:ncol(mod_10_spectra)]<-all_summary$File
  for (i in 2:ncol(mod_10_ratio)) {
    mod_10_ratio[,i]<-str_remove_all(mod_10_ratio[,i],"^\\(")
    mod_10_ratio[,i]<-as.numeric(str_remove_all(mod_10_ratio[,i],"%\\)$"))
    mod_10_occu[,i]<-as.numeric(str_remove_all(mod_10_occu[,i],"%$"))
    mod_10_spectra[,i]<-as.numeric(mod_10_spectra[,i])
  }
  mod_10_ratio2<-mod_10_ratio[,-1]
  rownames(mod_10_ratio2)<-mod_10_ratio$Modification
  library(pheatmap)
  png(paste0(outputpath,"mod_10_ratio.png"),width=(380+10*length(waitread)),height = 500)
  pp<-pheatmap(mod_10_ratio2,cluster_cols = FALSE,na_col = "gray",cluster_rows = FALSE)
  print(pp)
  dev.off()
  mod_10_occu2<-mod_10_occu[,-1]
  rownames(mod_10_occu2)<-mod_10_occu$Modification
  png(paste0(outputpath,"mod_10_occupancy.png"),width=(380+10*length(waitread)),height = 500)
  pp1<-pheatmap(mod_10_occu2,cluster_cols = FALSE,na_col = "gray",cluster_rows = FALSE)
  print(pp1)
  dev.off()
  mod_10_spectra2<-mod_10_spectra[,-1]
  rownames(mod_10_spectra2)<-mod_10_spectra$Modification
  png(paste0(outputpath,"mod_10_spectra.png"),width=(380+10*length(waitread)),height = 500)
  pp2<-pheatmap(mod_10_spectra2,cluster_cols = FALSE,na_col = "gray",cluster_rows = FALSE,scale="row")
  print(pp2)
  dev.off()
  #charge
  charge_all<-data.frame()
  charge_all<-as.data.frame(arrang(all_summary$Charge[1]))
  colnames(charge_all)<-c("Charge","Spectra Count","Ratio")
  charge_all<-charge_all[,c(1,3)]
  for (i in 2:nrow(all_summary)) {
    charge_sub<-data.frame(arrang(all_summary$Charge[i]))
    colnames(charge_sub)<-c("Charge","Spectra Count","Ratio")
    charge_all<-merge(charge_all,charge_sub[,c(1,3)],by.x = "Charge",by.y = "Charge",all = TRUE)
  }
  colnames(charge_all)[2:ncol(charge_all)]<-all_summary$File
  charge_all<-charge_all[-1,]
  for (i in 2:ncol(charge_all)) {
    charge_all[,i]<-str_replace(charge_all[,i],"%$","") %>% as.numeric()
  }
  charge_all<-reshape2::melt(charge_all,id.vars = "Charge")
  png(paste0(outputpath,"Misscleavage.png"),width=width,height = 450)
  p9<-ggplot(charge_all,mapping = aes(x=factor(variable,levels =unique(variable)),y=value,fill=factor(Charge,levels=c("6","5","4","3","2"))))+geom_bar(stat="identity",position = "stack")+ggtheme
  print(p9)
  dev.off()
  ###Misscleavage
  Misscleavage<-data.frame()
  Misscleavage<-as.data.frame(arrang(all_summary$Misscleavage[1]))
  colnames(Misscleavage)<-c("Misscleavage","Spectra Count","Ratio")
  Misscleavage<-Misscleavage[,c(1,3)]
  for (i in 2:nrow(all_summary)) {
    Misscleavage_sub<-data.frame(arrang(all_summary$Misscleavage[i]))
    colnames(Misscleavage_sub)<-c("Misscleavage","Spectra Count","Ratio")
    Misscleavage<-merge(Misscleavage,Misscleavage_sub[,c(1,3)],by.x = "Misscleavage",by.y = "Misscleavage",all = TRUE)
  }
  colnames(Misscleavage)[2:ncol(Misscleavage)]<-all_summary$File
  for (i in 2:ncol(Misscleavage)) {
    Misscleavage[,i]<-str_replace(Misscleavage[,i],"%$","") %>% as.numeric()
  }
  Misscleavage<-reshape2::melt(Misscleavage,id.vars = "Misscleavage")
  png(paste0(outputpath,"Misscleavage.png"),width=width,height = 450)
  p8<-ggplot(Misscleavage,mapping = aes(x=factor(variable,levels =unique(variable)),y=value,fill=factor(Misscleavage,levels=c("4","3","2","1","0"))))+geom_bar(stat="identity",position = "stack")+ggtheme
  print(p8)
  dev.off()
}




###pfind结果读取
summary<-fread("Z:/ICU/-a1/pFindTask5/result/pFind.summary",sep="\n")

file_summary<-function(path){
  #  path<-
  summary<-fread(path,sep="\n",data.table = FALSE)
  nn<-which(summary[,1]=="Modifications:")+2
  mod<-paste0(summary[(nn+1):(nn+9),1],collapse = "|")
  mod<-str_replace_all(mod,"\\\t"," ")
  nncharge<-which(summary=="Charge Distribution of Peptides:")
  charge<-paste0(summary[(nncharge+1):(nncharge+5),1],collapse = "|")
  charge<-str_replace_all(charge,"\\\t"," ")
  nnmiss<-which(summary=="Missed Cleavages of Peptides:")
  Miss<-paste0(summary[(nnmiss+1):(nnmiss+5),1],collapse = "|")
  Miss<-str_replace_all(Miss,"\\\t"," ")
  nnmix<-which(summary=="Mixed Spectra of MS/MS Scans:")
  mix<-paste0(summary[(nnmix+1):(nnmix+6),1],collapse = "|")
  mix<-str_replace_all(mix,"\\\t"," ") 
  error<-which(summary=="Average Ion Errors:")
  c<-c(as.numeric(str_extract(summary[1,],"(?<=Spectra: )[0-9]{1,}")),
       as.numeric(str_extract(summary[2,],"(?<=Scans: )[0-9]{1,}")),
       as.numeric(str_extract(summary[3,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[4,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[5,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[6,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[7,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[8,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[9,],"(?<=\\=)[0-9]{1,}\\.[0-9]{1,}")),
       as.numeric(str_extract(summary[10,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[11,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[12,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[15,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[16,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[17,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[18,],"[0-9]{1,}")),
       mod,
       charge,
       Miss,
       mix,
       str_extract(summary[error+2,1],"(?<=Precursor_Mean:	)[0-9]{1,}\\.[0-9]{1,}"),
       str_extract(summary[error+2,1],"(?<=Precursor_Sd:	)[0-9]{1,}\\.[0-9]{1,}"),
       str_extract(summary[error+3,1],"(?<=Fragment_Mean:	)[0-9]{1,}\\.[0-9]{1,}"),
       str_extract(summary[error+3,1],"(?<=Fragment_Sd:	)[0-9]{1,}\\.[0-9]{1,}")
  )
  return(c)
}
file_summary2<-function(path){
  #  path<-
  summary<-fread(path,sep="\n",data.table = FALSE)
  nn<-which(summary[,1]=="Modifications:")+2
  mod<-c()
  for(can in 1:length(candidate)){
    mod<-c(mod,summary[grep(candidate[can],summary[,1]),1])
  }
  mod<-paste0(mod,collapse = "|")
  mod<-str_replace_all(mod,"\\\t"," ")
  nncharge<-which(summary=="Charge Distribution of Peptides:")
  charge<-paste0(summary[(nncharge+1):(nncharge+5),1],collapse = "|")
  charge<-str_replace_all(charge,"\\\t"," ")
  nnmiss<-which(summary=="Missed Cleavages of Peptides:")
  Miss<-paste0(summary[(nnmiss+1):(nnmiss+5),1],collapse = "|")
  Miss<-str_replace_all(Miss,"\\\t"," ")
  nnmix<-which(summary=="Mixed Spectra of MS/MS Scans:")
  mix<-paste0(summary[(nnmix+1):(nnmix+6),1],collapse = "|")
  mix<-str_replace_all(mix,"\\\t"," ") 
  error<-which(summary=="Average Ion Errors:")
  c<-c(as.numeric(str_extract(summary[1,],"(?<=Spectra: )[0-9]{1,}")),
       as.numeric(str_extract(summary[2,],"(?<=Scans: )[0-9]{1,}")),
       as.numeric(str_extract(summary[3,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[4,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[5,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[6,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[7,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[8,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[9,],"(?<=\\=)[0-9]{1,}\\.[0-9]{1,}")),
       as.numeric(str_extract(summary[10,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[11,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[12,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[15,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[16,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[17,],"[0-9]{1,}")),
       as.numeric(str_extract(summary[18,],"[0-9]{1,}")),
       mod,
       charge,
       Miss,
       mix,
       str_extract(summary[error+2,1],"(?<=Precursor_Mean:	)[0-9]{1,}\\.[0-9]{1,}"),
       str_extract(summary[error+2,1],"(?<=Precursor_Sd:	)[0-9]{1,}\\.[0-9]{1,}"),
       str_extract(summary[error+3,1],"(?<=Fragment_Mean:	)[0-9]{1,}\\.[0-9]{1,}"),
       str_extract(summary[error+3,1],"(?<=Fragment_Sd:	)[0-9]{1,}\\.[0-9]{1,}")
  )
  return(c)
}
LPS<-paste0("Y:/LH/LH/20220201_LPS/",list.files("Y:/LH/LH/20220201_LPS"),"/pFind.summary")
LPS<-LPS[-48]
LPS<-c(LPS,paste0("Y:/LH/LH/20220827_LPS0-6h/",list.files("Y:/LH/LH/20220827_LPS0-6h"),"/pFind.summary"))
LPS<-LPS[-90]
all_summary<-data.frame()
for (i in 1:length(LPS)) {
  all_summary<-rbind(all_summary,file_summary(LPS[i]))
}
colnames(all_summary)<-c("Spectra","Scans","Peptides","Sequence","Proteins","Protein Groups","Avg_Seq_Cov","Decoy spectra","T/D Distribution","Decoy Peptides","Decoy Proteins","Decoy Protein Groups","Specific","C-term Specific","N-term Specific","Non-Specific","Modification","Charge","Misscleavage","Mixture spectra","Precursor Mean","Precursor SD","Fragment Mean","Fragment SD")
all_summary<-cbind(File=str_extract(LPS,"(?<=[LPS|LPS0-6h]{1}/).*(?=_pFindResults/pFind.summary)"),all_summary)
all_summary$Group<-"New"
all_summary$Group[1:47]<-"Old"
library(ggplot2)
library(ggsci)
library(ggthemes)
library(reshape2)
ggtheme=theme_light()+theme(legend.title = element_blank(),
                            axis.text.y = element_text(colour="black", size=14,face = "bold",family = "serif"),
                            axis.text.x = element_text(colour="black", size=14,face = "bold",family = "serif",angle = 60,hjust = 1),
                            axis.title = element_text(colour="black", size=14,face = "bold",family = "serif"),
                            legend.position = "right",title =element_text(size=14,face = "bold",family = "serif"),
                            legend.text = element_text(size=14,face = "bold",family = "serif"))

ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Peptides)/as.numeric(Scans),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Peptide/Scans")+scale_y_continuous(limits = c(0,0.8),breaks = c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0))
ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Peptides)/as.numeric(Spectra),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Peptide/Spectra")+scale_y_continuous(limits = c(0,0.8),breaks = c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0))
ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Peptides),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Peptide")+scale_y_continuous(expand = c(0,0))
ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(`Protein Groups`),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Protein Group")+scale_y_continuous(expand = c(0,0))
ggplot(all_summary,mapping = aes(x=factor(File,levels=File),y=as.numeric(Avg_Seq_Cov),fill=Group))+geom_bar(stat="identity")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Avg_Seq_Cov")+scale_y_continuous(expand = c(0,0))
all_summary.sub<-all_summary[,c(1,14:17)]
for (i in 2:5) {
  all_summary.sub[,i]<-as.numeric(all_summary.sub[,i])
}
all_summary.sub$Sum<-apply(all_summary.sub[,2:5],1,sum)
for (i in 2:5) {
  all_summary.sub[,i]<-as.numeric(all_summary.sub[,i])/as.numeric(all_summary.sub[,6])
}
all_summary.sub<-all_summary.sub[,-6]
all_summary.sub<-melt(all_summary.sub,id.vars = "File")
ggplot(all_summary.sub,mapping = aes(x=factor(File,levels=unique(File)),y=value,fill=factor(variable,levels = c("Non-Specific","N-term Specific","C-term Specific","Specific"))))+geom_bar(stat="identity",position = "stack")+ggtheme+scale_fill_tableau()+xlab("")+ylab("Enzymatic")+scale_y_continuous(expand = c(0,0))

arrang<-function(x){
  tab<-data.frame()
  a<-unlist(str_split(x,"\\|"))
  for (n in 1:length(a)) {
    tab<-rbind(tab,unlist(str_split(a[n]," ")))
  }
  return(tab)
}


####mod summary
mod_10<-data.frame()
mod_sub<-data.frame()
mod_10<-as.data.frame(arrang(all_summary$Modification[1]))
colnames(mod_10)<-c("Modification","Spectra","Rate","Sites","Relative Ratio")
mod_10_ratio<-mod_10[,c(1,3)]
mod_10_occu<-mod_10[,c(1,5)]
mod_10_spectra<-mod_10[,c(1,2)]
for (i in 2:nrow(all_summary)) {
  mod_sub<-data.frame(arrang(all_summary$Modification[i]))
  colnames(mod_sub)<-c("Modification","Spectra","Rate","Sites","Relative Ratio")
  mod_10_ratio<-merge(mod_10_ratio,mod_sub[,c(1,3)],by.x = "Modification",by.y = "Modification",all = TRUE)
  mod_10_occu<-merge(mod_10_occu,mod_sub[,c(1,5)],by.x = "Modification",by.y = "Modification",all = TRUE)
  mod_10_spectra<-merge(mod_10_spectra,mod_sub[,c(1,2)],by.x = "Modification",by.y = "Modification",all = TRUE)
}
colnames(mod_10_ratio)[2:90]<-all_summary$File
colnames(mod_10_occu)[2:90]<-all_summary$File
colnames(mod_10_spectra)[2:90]<-all_summary$File


for (i in 2:ncol(mod_10_ratio)) {
  mod_10_ratio[,i]<-str_remove_all(mod_10_ratio[,i],"^\\(")
  mod_10_ratio[,i]<-as.numeric(str_remove_all(mod_10_ratio[,i],"%\\)$"))
  mod_10_occu[,i]<-as.numeric(str_remove_all(mod_10_occu[,i],"%$"))
  mod_10_spectra[,i]<-as.numeric(mod_10_spectra[,i])
}
mod_10_ratio2<-mod_10_ratio[,-1]
rownames(mod_10_ratio2)<-mod_10_ratio$Modification
library(pheatmap)
pheatmap(mod_10_ratio2,cluster_cols = FALSE,na_col = "gray",cluster_rows = FALSE)

mod_10_occu2<-mod_10_occu[,-1]
rownames(mod_10_occu2)<-mod_10_occu$Modification
pheatmap(mod_10_occu2,cluster_cols = FALSE,na_col = "gray",cluster_rows = FALSE)

mod_10_spectra2<-mod_10_spectra[,-1]
rownames(mod_10_spectra2)<-mod_10_spectra$Modification
pheatmap(mod_10_spectra2,cluster_cols = FALSE,na_col = "gray",cluster_rows = FALSE,scale="row")



candidate<-rownames(mod_10_ratio2)


###charge summary
charge_all<-data.frame()
charge_all<-as.data.frame(arrang(all_summary$Charge[1]))
colnames(charge_all)<-c("Charge","Spectra Count","Ratio")
charge_all<-charge_all[,c(1,3)]
for (i in 2:nrow(all_summary)) {
  charge_sub<-data.frame(arrang(all_summary$Charge[i]))
  colnames(charge_sub)<-c("Charge","Spectra Count","Ratio")
  charge_all<-merge(charge_all,charge_sub[,c(1,3)],by.x = "Charge",by.y = "Charge",all = TRUE)
}
colnames(charge_all)[2:90]<-all_summary$File
charge_all<-charge_all[-1,]
for (i in 2:ncol(charge_all)) {
  charge_all[,i]<-str_replace(charge_all[,i],"%$","") %>% as.numeric()
}
charge_all<-melt(charge_all,id.vars = "Charge")
ggplot(charge_all,mapping = aes(x=factor(variable,levels =unique(variable)),y=value,fill=factor(Charge,levels=c("6","5","4","3","2"))))+geom_bar(stat="identity",position = "stack")+ggtheme

###Misscleavage
Misscleavage<-data.frame()
Misscleavage<-as.data.frame(arrang(all_summary$Misscleavage[1]))
colnames(Misscleavage)<-c("Misscleavage","Spectra Count","Ratio")
Misscleavage<-Misscleavage[,c(1,3)]
for (i in 2:nrow(all_summary)) {
  Misscleavage_sub<-data.frame(arrang(all_summary$Misscleavage[i]))
  colnames(Misscleavage_sub)<-c("Misscleavage","Spectra Count","Ratio")
  Misscleavage<-merge(Misscleavage,Misscleavage_sub[,c(1,3)],by.x = "Misscleavage",by.y = "Misscleavage",all = TRUE)
}
colnames(Misscleavage)[2:90]<-all_summary$File
for (i in 2:ncol(Misscleavage)) {
  Misscleavage[,i]<-str_replace(Misscleavage[,i],"%$","") %>% as.numeric()
}
Misscleavage<-melt(Misggplot(Misscleavage,mapping = aes(x=factor(variable,levels =unique(variable)),y=value,fill=factor(Misscleavage,levels=c("4","3","2","1","0"))))+geom_bar(stat="identity",position = "stack")+ggtheme
                   scleavage,id.vars = "Misscleavage")
