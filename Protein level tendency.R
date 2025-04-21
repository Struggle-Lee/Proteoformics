#####Modification_deal
Final_G0<-list()
for (i in 1:6) {
  Final_G0[[i]]<-fread(path,data.table = FALSE)
}
Final_G0_MEMBER<-list()
Final_G0_LogP<-list()
Index.Info<-c()
for (i in 1:6) {
  #if(i>1){
  #  remove.colname<-match(colnames(Final_G0[[i-1]])[1:46],colnames(Final_G0[[i]])[1:48])
  #  remove.colname<-remove.colname[which(!is.na(remove.colname))]
  #  Final_G0[[i]]<-Final_G0[[i]][,-remove.colname]
  #}
  if(i<6){
    Final_G0[[i]]<-Final_G0[[i]][,-c(67:73)]
  }else{
    Final_G0[[i]]<-Final_G0[[i]][,-c(51:57)]
  }
}
for (i in 1:6) { 
  Final_G0[[i]]$Info<-paste0(Final_G0[[i]]$GO,"_",Final_G0[[i]]$Description)
  Index.Info<-c(Index.Info,Final_G0[[i]]$Info)
  if(i!=6){
    Final_G0_MEMBER[[i]]<-data.frame(Info=Final_G0[[i]]$Info,Final_G0[[i]][,1:23])
    Final_G0_LogP[[i]]<-data.frame(Info=Final_G0[[i]]$Info,Final_G0[[i]][,24:46])
    Final_G0[[i]]<-Final_G0[[i]][-c(1:46)]
  }else{
    Final_G0_MEMBER[[i]]<-data.frame(Info=Final_G0[[i]]$Info,Final_G0[[i]][,1:15])
    Final_G0_LogP[[i]]<-data.frame(Info=Final_G0[[i]]$Info,Final_G0[[i]][,16:30])
    Final_G0[[i]]<-Final_G0[[i]][-c(1:30)]  
  }
}
Index.Info<-unique(Index.Info)

Arrange<-data.frame(Info=Index.Info,GO=str_extract(Index.Info,".*[0-9]{2,}(?=_)"),Description=str_extract(Index.Info,"(?<=_).*"))  
for (i in 1:6) {
  Arrange<-cbind(Arrange,Final_G0_LogP[[i]][match(Arrange$Info,Final_G0_LogP[[i]]$Info),-1])    
}  
for (i in 1:6) {
  Arrange<-cbind(Arrange,Final_G0_MEMBER[[i]][match(Arrange$Info,Final_G0_MEMBER[[i]]$Info),-1])    
}  
Final_GO1<-Final_G0[[1]]
for (i in 2:6) {
  Final_GO1<-rbind(Final_GO1,Final_G0[[i]])
}

Arrange$Rank<-0
Arrange$TotalGeneInLibrary<-30293
Arrange$GeneinGO<-0
Arrange$GeneInHitList<-0
Arrange$GeneInGOAndHitList<-0
Arrange$'%InGO'<-0
Arrange$LogP<-0
Arrange$LogQ<-0
Arrange$Enrichment<-0
Arrange$Category<-""
Arrange$GeneID<-""
Arrange$Hits<-""
for (i in 1:length(Index.Info)) {
  sub<-Final_GO1[which(Final_GO1$Info==Index.Info[i]),]
  Arrange$GeneinGO[i]<-sub$`#GeneInGO`[1]
  Arrange$Category[i]<-sub$Category[1]
  Arrange$GeneID[i]<-paste0(unique(unlist(str_split(sub$GeneID,"\\|"))),collapse = "|")
  Arrange$Hits[i]<-paste0(unique(unlist(str_split(sub$Hits,"\\|"))),collapse = "|")
}
Arrange$GeneInHitList<-length(unique(unlist(str_split(Arrange$GeneID,"\\|"))))
Arrange$GeneInGOAndHitList<-apply(Arrange[,273:275],1,function(x){str_count(x[2],"\\|")+1})
Arrange$Enrichment<-(Arrange$GeneInGOAndHitList/Arrange$GeneInHitList)/(Arrange$GeneinGO/Arrange$TotalGeneInLibrary)
Arrange$Rank<-apply(Arrange[,134:263],1,function(x){length(which(x==1))})
Arrange$`%InGO`<-(Arrange$GeneInGOAndHitList/Arrange$GeneinGO)*100
Arrange$LogP<-apply(Arrange[,265:268],1,function(x){log10(phyper(x[4]-1,x[2],x[1]-x[2],x[3],lower.tail = FALSE))})
Arrange$LogQ<-log10(p.adjust(10^Arrange$LogP,"fdr"))

Arrange1<-Arrange[which(Arrange$LogQ<log10(0.01)),]

# 数据是列表，每个元素是一个向量，包含了一个术语的所有目标基因
# list_of_terms 是你的数据
list_of_terms <- str_split(Arrange1$Hits,"\\|")

# 计算共现矩阵
co_occurrence_matrix <- matrix(nrow = length(list_of_terms), ncol = length(list_of_terms))
for (i in 1:length(list_of_terms)) {
  for (j in i:length(list_of_terms)) {
    co_occurrence_matrix[i, j] <- length(intersect(list_of_terms[[i]], list_of_terms[[j]]))
    co_occurrence_matrix[j, i] <- co_occurrence_matrix[i, j]
  }
}

# 计算kappa统计量
library(vcd)
kappa_matrix <- matrix(1,nrow = length(list_of_terms), ncol = length(list_of_terms))
for (i in 1:(length(list_of_terms) - 1)) {
  for (j in (i + 1):length(list_of_terms)) {
    contingency_table <- matrix(c(co_occurrence_matrix[i, j],
                                  length(list_of_terms[[i]]) - co_occurrence_matrix[i, j],
                                  length(list_of_terms[[j]]) - co_occurrence_matrix[i, j],
                                  length(unique(unlist(list_of_terms))) - length(list_of_terms[[i]]) - length(list_of_terms[[j]]) + co_occurrence_matrix[i, j]),
                                nrow = 2)
    kappa_matrix[i, j] <- Kappa(contingency_table,weights ="Fleiss-Cohen")$Weighted[1]
    kappa_matrix[j, i] <- kappa_matrix[i, j]
  }
  print(i)
}

dist_matrix <- 1 - kappa_matrix

# 进行层次聚类
hc <- cutree(hclust(as.dist(dist_matrix),method = "average"),h=0.7)
Arrange1$GroupID<-cutree(hclust(as.dist(dist_matrix),method = "average"),h=0.7)  
length(unique(Arrange1$GroupID))
Arrange.sub<-data.frame()  
for (i in 1:215) {
  sub<-Arrange1[which(Arrange1$GroupID==i),]
  sub<-sub[order(sub$LogP),]
  Arrange.sub<-rbind(Arrange.sub,sub[1,])
}  
for (i in 1:ncol(Arrange1)) {
  Arrange1[which(is.na(Arrange1[,i])),i]<-0
}
for (i in 4:263) {
  Arrange.sub[which(is.na(Arrange.sub[,i])),i]<-0
}

Arrangeforheat<-Arrange.sub[,4:133]
rownames(Arrangeforheat)<-Arrange.sub$Description
for (i in 1:ncol(Arrangeforheat)) {
  Arrangeforheat[is.na(Arrangeforheat[,i]),i]<-0
  Arrangeforheat[,i]<-Arrangeforheat[,i]/min(Arrangeforheat[,i])
}
library(pheatmap)  
pheatmap(Arrangeforheat,scale = "none")  

Arrangeforheat.sub<-Arrangeforheat[-which(apply(Arrangeforheat,1,function(x){length(which(x>=0.25))})<=1),]
rownames(Arrangeforheat.sub)[6]<-"Regulation of IGF transport and uptake by IGFBPs"  
Arrangeforheat.sub<-Arrangeforheat.sub[which(apply(Arrangeforheat.sub,1,max)>=0.25),]

p<-pheatmap::pheatmap(Arrangeforheat.sub[-1,-2],scale = "none",color = colorRampPalette(c("gray95","#0079c3","#004b7e","#fcd475","#f08300","firebrick","darkred"))(60),clustering_method = "ward.D",clustering_distance_cols = "correlation",cutree_rows = 6,cutree_cols = 7)  
p<-pheatmap::pheatmap(Arrangeforheat.sub,scale = "none",color = colorRampPalette(c("gray95","#0079c3","#004b7e","#fcd475","#f08300","firebrick","darkred"))(60),clustering_method = "ward.D",clustering_distance_cols = "manhattan",cutree_rows = 6,cutree_cols = 7)

library(circlize)
library(ComplexHeatmap)
library(dendextend)
library(dendsort)
library(gridBase)
mycol1=colorRamp2(c(0,0.1,0.3,0.5,0.7,0.9,1),c("gray95","#0079c3","#004b7e","#fcd475","#f08300","firebrick","darkred"))
cir.ptm<-t(Arrangeforheat.sub[,-2])
cir.ptm<-cir.ptm[,p$tree_row$order[c(11:5,59:12,4:1)]]
cir.ptm<-cir.ptm[p$tree_col$order,]
circos.par(start.degree = 180, gap.degree = 90)
split<-factor(cutree(p$tree_col,k = 6)[p$tree_col$order],levels=1:6)
circos.heatmap(cir.ptm,col=mycol1,
               dend.side="inside",# dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               rownames.side="outside",# rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               track.height = 0.4,
               rownames.col="black",
               bg.border="transparent",
               #split = ann_row,
               #show.sector.labels = FALSE,
               rownames.cex=0.6, # 字体大小
               rownames.font=1, # 字体粗细
               cluster=TRUE,distance.method = "manhattan",clustering.method = "ward.D",
               dend.track.height=0.05,#调整行聚类树的高度
               dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，
                 #或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
                 color_branches(dend,k=6,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
               }
)
circos.clear()
lg=Legend(title="Legend",col_fun=mycol1,direction = c("horizontal"))
grid.draw(lg)