#trend cluster
library(Mfuzz)
mfuzz1<-data.frame(Control=apply(heat.Lysine.ALB[,13:47],1,function(x){log2(mean(as.numeric(x^2)))}),
                   Sepsis=apply(heat.Lysine.ALB[,48:108],1,function(x){log2(mean(as.numeric(x^2)))}),
                   NonSurvivor=apply(heat.Lysine.ALB[,109:161],1,function(x){log2(mean(as.numeric(x^2)))}))
mfuzz_class1<-new('ExpressionSet',exprs=as.matrix(mfuzz1))
mfuzz_class1<-standardise(mfuzz_class1)
mfuzz_cluster1 <- mfuzz(mfuzz_class1, c = 3, m = mestimate(mfuzz_class1))
mfuzz.plot2(mfuzz_class1, cl = mfuzz_cluster1, mfrow = c(2, 2),centre=TRUE, centre.col="black",centre.lwd=2,time.labels = colnames(mfuzz_class1),x11 = FALSE,ylab="log2(Ratio)",ylim.set=c(-1.5,1.5))

mfuzz.p<-rbind(
  heat.Lysine.ALB[match(names(which(mfuzz_cluster1$cluster==1)),heat.Lysine.ALB$Info),],
  heat.Lysine.ALB[match(names(which(mfuzz_cluster1$cluster==2)),heat.Lysine.ALB$Info),],
  heat.Lysine.ALB[match(names(which(mfuzz_cluster1$cluster==3)),heat.Lysine.ALB$Info),]
)
pheatmap(mfuzz.p[,13:161],scale = "row",color = colorRampPalette(c("#004b7e","#004b7e","#004b7e","#0079c3","gray95","#f38600","darkred","darkred","darkred"))(60),cluster_cols = FALSE,clustering_distance_rows = "correlation",clustering_method = "ward.D2",cluster_rows = FALSE,annotation = Annotation.r,gaps_col = c(35,96,149),annotation_colors = ann_colors[2],show_colnames = FALSE)

plot<-as.data.frame(table(Mod_Info.sub2.ALB$Mod_pos[which(Mod_Info.sub2.ALB$AA=="K")]))
plot<-merge(plot,as.data.frame(table(Mod_Info.sub2.ALB$Mod_pos[which(Mod_Info.sub2.ALB$AA=="K" & correlation.P2.ALB$Trend=="Different")])),by="Var1",all=T)
plot<-merge(plot,as.data.frame(table(Mod_Info.sub2.ALB$Mod_pos[which(Mod_Info.sub2.ALB$AA=="K" & correlation.P2.ALB$Trend=="Same")])),by="Var1",all=T)

colnames(plot)<-c("Pos","All","Different","Same")
plot<-plot[which(plot$All>=3),]
plot<-plot[order(plot$All,decreasing = T),]
plot$Different[which(is.na(plot$Different))]<-0
plot$Same[which(is.na(plot$Same))]<-0
plot$Position<-paste0("Lys",as.numeric(as.character(plot$Pos))-24)
plot$Pos<-as.character(plot$Pos)
plot1<-plot[,c(5:3)]
plot1<-reshape2::melt(plot1,id.vars="Position")
ggplot(plot1,mapping = aes(x=variable,y=factor(Position,levels=unique(Position)),size=value,colour=value))+
  geom_point()+ggtheme+scale_size_continuous(limits = c(1,8),range = c(1,8))+
  scale_color_gradientn(colors = c("#fcd475", "#f08300", "firebrick"))+xlab("")+ylab("Sites")

plot.C<-as.data.frame(table(Mod_Info.sub2.ALB.C$Mod_pos))
plot.C$Var1<-as.numeric(as.character(plot.C$Var1))
plot.C$Info<-paste0("ALBU_",plot.C$Var1,"C$")
plot.C$Count<-0
for (i in 1:nrow(plot.C)) {
  plot.C$Count[i]<-length(which(str_detect(Mod_Info.sub2$Info,plot.C$Info[i])))
}

plot.C$Var1<-paste0("Cys",plot.C$Var1-24)
sub<-as.data.frame(table(Mod_Info.sub2.ALB.Csub$Mod_pos))
sub$Var1<-as.numeric(as.character(sub$Var1))
sub$Var1<-paste0("Cys",sub$Var1-24)
plot.C<-merge(plot.C,sub,by="Var1",all=T)
plot.C$Freq.y[is.na(plot.C$Freq.y)]<-0
plot.C<-plot.C[which(plot.C$Freq.x>=3),]
colnames(plot.C)<-c("Site","Total","Info","Count","Different")
plot.C$Same<-plot.C$Total-plot.C$Different
plot.C$Pos<-as.numeric(str_remove(plot.C$Site,"Cys"))
plot.C<-plot.C[order(plot.C$Pos),]
plot.C<-plot.C[,c(1,4:6)]
rownames(plot.C)<-plot.C$Site
pheatmap(plot.C[,c(4:6)],scale = "column",cluster_cols = FALSE,cluster_rows = FALSE,border_color = "transparent",gaps_row = 1,gaps_col = 1)
plot.C<-reshape2::melt(plot.C,id.vars=c("Site","Pos","Total"))
ggplot(plot.C,mapping = aes(x=variable,y=factor(Site,levels=unique(Site)),size=value,colour=value))+
  geom_point()+ggtheme+scale_size_continuous(limits = c(1,11),range = c(1,11))+
  scale_color_gradientn(colors = c("#fcd475", "#f08300", "firebrick"))+xlab("")+ylab("Sites")

C.Mod<-as.data.frame(table(correlation.P2.ALB.C$Mod))
sub<-as.data.frame(table(correlation.P2.ALB.C$Mod[which(correlation.P2.ALB.C$Trend=="Different")]))
C.Mod$Var1<-as.character(C.Mod$Var1)
C.Mod<-C.Mod[order(C.Mod$Freq,decreasing = T),]
sub$Var1<-as.character(sub$Var1)
C.Mod<-merge(C.Mod,sub,by="Var1",all=T)
colnames(C.Mod)<-c("Modification","Total","Different")
C.Mod$Different[is.na(C.Mod$Different)]<-0
C.Mod$Same<-C.Mod$Total-C.Mod$Different
C.Mod<-C.Mod[which(C.Mod$Total>=3),]
C.Mod<-reshape2::melt(C.Mod,id.vars=c("Modification","Total"))
C.Mod<-C.Mod[-c(16,47),]
ggplot(C.Mod,mapping = aes(x=variable,y=factor(Modification,levels=unique(Modification)),size=value,colour=value))+
  geom_point()+ggtheme+scale_size_continuous(limits = c(0,26),range = c(0,8))+
  scale_color_gradientn(colors = c("#fcd475", "#f08300", "firebrick"))+xlab("")+ylab("Modifications")
