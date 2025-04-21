###Mod_protein level crosstalk
protein_list <- strsplit(as.character(mod_Info_enrich.sub$Uniprot), " ")
all_proteins <- unique(Mod_Info.sub2$UniprotID)
binary_matrix <- sapply(protein_list, function(proteins) {
  as.integer(all_proteins %in% proteins)
})
rownames(binary_matrix) <- all_proteins
colnames(binary_matrix)<-mod_Info_enrich.sub$MOD

kappa_matrix <- matrix(NA, nrow = ncol(binary_matrix), ncol = ncol(binary_matrix))
rownames(kappa_matrix) <- mod_Info_enrich.sub$Modification
colnames(kappa_matrix) <- mod_Info_enrich.sub$Modification

for (i in 1:(ncol(binary_matrix) - 1)) {
  for (j in (i + 1):ncol(binary_matrix)) {
    table_data <- table(binary_matrix[, i], binary_matrix[, j])
    kappa_value <- Kappa(table_data)$Unweighted[1]
    kappa_matrix[i, j] <- kappa_value
    kappa_matrix[j, i] <- kappa_value
  }
}
for (i in 1:(ncol(kappa_matrix))) {
  kappa_matrix[i,i]<-0
}
rownames(kappa_matrix)<-mod_Info_enrich.sub$MOD
colnames(kappa_matrix)<-mod_Info_enrich.sub$MOD
kappa_matrix[abs(kappa_matrix) < 0.35] <- 0
non_zero_indices <- rowSums(abs(kappa_matrix) > 0) > 0
kappa_matrix <- kappa_matrix[non_zero_indices, non_zero_indices]
library(igraph)
# 将 kappa_matrix 转换为图
#kappa_graph <- graph_from_adjacency_matrix(kappa_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
# 绘制图形
#plot(kappa_graph, 
#     edge.width = E(kappa_graph)$weight * 5, # 根据权重调整边宽
#     vertex.label = V(kappa_graph)$name,
#     vertex.size = 15,
#     main = "Modification Interaction Network")

edge_list <- which(kappa_matrix != 0, arr.ind = TRUE)
edge_list <- data.frame(
  source = rownames(kappa_matrix)[edge_list[, 1]],
  target = colnames(kappa_matrix)[edge_list[, 2]],
  weight = kappa_matrix[edge_list])
# 导出为 CSV 文件
write.csv(edge_list, "kappa_edges_0.35.csv", row.names = FALSE)

nodes <- unique(c(edge_list$source, edge_list$target))
# 创建节点属性数据框
node_attributes <- data.frame(
  node = nodes,
  attribute = "YourAttribute"  # 例如，修饰类型或其他信息
)

# 导出为 CSV 文件
write.csv(node_attributes, "node_attributes.csv", row.names = FALSE)




###domain infomation
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
protein_list <- unique(c(Sepsis.mod$UniprotID,str_extract(Sepsis.unmod$Protein,"(?<=sp\\|)[A-Z0-9]{4,10}(?=\\|)")))
domains <- getBM(attributes = c("uniprotswissprot", "interpro", "interpro_description","interpro_short_description","interpro_start","interpro_end"),
                 filters = "uniprotswissprot",
                 values = protein_list,
                 mart = mart)

domains.i <- getBM(attributes = c("uniprotswissprot", "interpro","superfamily","superfamily_start","superfamily_end"),
                   filters = "uniprotswissprot",
                   values = protein_list,
                   mart = mart)
domains$interpro_description1<-str_remove(domains$interpro_description,"(?<=,) .*")
domains$interpro_description1<-str_remove(domains$interpro_description1,",$")
domains$interpro_description1<-str_remove(domains$interpro_description1," superfamily$")
domains$interpro_description1<-str_remove(domains$interpro_description1,".containing domain$")
domains$interpro_description1<-str_remove(domains$interpro_description1,".like domain$")

domains.sub<-data.frame(Name=unique(domains$interpro_description1))
write.csv(domains.sub,"Domain1.csv",row.names = FALSE)

library(openxlsx)
domains.name<-read.xlsx("Domain1.xlsx",sheet = 1)
domains$Name<-domains.name$X2[match(domains$interpro_description1,domains.name$Domain)]

fasta.sep1$length<-apply(fasta.sep1,1,function(x){nchar(x[3])})
domains$length<-fasta.sep1$length[match(domains$uniprotswissprot,fasta.sep1$UniprotID)]
domains$length.real<-domains$interpro_end-domains$interpro_start+1
domains.sub<-domains[which(domains$length.real<=domains$length*0.9 & !is.na(domains$Name)),]
domains.family<-data.frame(Name=unique(domains.sub$Name))
domains.family$totaln<-0
domains.family$Pos<-""
domains.family$totalinter<-0

for (i in 1:nrow(domains.family)) {
  dsub<-domains.sub[which(domains.sub$Name==domains.family$Name[i]),]
  dsub.c<-c()
  for (x in 1:nrow(dsub)) {
    dsub.c<-c(dsub.c,paste0(dsub$uniprotswissprot[x],"_",seq(dsub$interpro_start[x],dsub$interpro_end[x],1))) 
  }
  dsub.c<-unique(dsub.c)
  domains.family$totaln[i]<-length(dsub.c)
  domains.family$Pos[i]<-paste0(dsub.c,collapse = ";")
  print(i)
}


Mod_Info.sub2$Domain<-apply(Mod_Info.sub2,1,function(x){
  sub<-domains.sub[which(domains.sub$uniprotswissprot==x[3] & domains.sub$interpro_end<=x[10] & domains.sub$interpro_start>=x[10]),]
  return(paste0(unique(sub$Name),collapse = ";"))
})


length(unique(c(paste0(Mod_Info.sub2$UniprotID,"_",Mod_Info.sub2$Sequence),paste0(str_extract(Sepsis.unmod$Protein,"(?<=sp\\|)[A-Z0-9]{4,10}(?=\\|)"),"_",Sepsis.unmod$Protein))))
list<-unique(c(paste0(Mod_Info.sub2$UniprotID,"_",Mod_Info.sub2$Sequence),paste0(str_extract(Sepsis.unmod$Protein,"(?<=sp\\|)[A-Z0-9]{4,10}(?=\\|)"),"_",Sepsis.unmod$Info)))
domain.identotal<-data.frame(ID=str_extract(list,".*(?=\\_)"),Sequence=str_extract(list,"(?<=\\_).*"))
domain.identotal$Start<-0
domain.identotal$End<-0
for (i in 1:nrow(domain.identotal)) {
  domain.identotal$Start[i]<-str_locate(str_replace_all(fasta.sep1$Seq[which(fasta.sep1$UniprotID==domain.identotal$ID[i])],"L","I"),str_replace_all(domain.identotal$Sequence[i],"L","I"))[1]
  domain.identotal$End[i]<-str_locate(str_replace_all(fasta.sep1$Seq[which(fasta.sep1$UniprotID==domain.identotal$ID[i])],"L","I"),str_replace_all(domain.identotal$Sequence[i],"L","I"))[2]
  print(i)
}
identotal<-c()
for (i in 1:nrow(domain.identotal)) {
  identotal<-c(identotal,paste0(domain.identotal$ID[i],"_",seq(domain.identotal$Start[i],domain.identotal$End[i],1)))
  identotal<-unique(identotal)
  print(i)
}
for (i in 1:nrow(domains.family)) {
  domains.family$totalinter[i]<-length(intersect(unlist(str_split(domains.family$Pos[i],";")),identotal))
}
domains.family$total<-163171
domains.family$mod.total<-14376
Mod_Info.sub2$Infopos<-paste0(Mod_Info.sub2$UniprotID,"_",Mod_Info.sub2$Mod_pos)
domains.family$modinter<-0
for (i in 1:nrow(domains.family)) {
  domains.family$modinter[i]<-length(intersect(unlist(str_split(domains.family$Pos[i],";")),Mod_Info.sub2$Infopos))
}
domains.family$ES<-0
domains.family$p_value<-0
# 计算超几何检验的P值
p_value <- phyper(
  domains.family$modinter - 1,
  domains.family$totalinter,
  domains.family$total-domains.family$totalinter,
  domains.family$mod.total,
  lower.tail = FALSE
)
domains.family$p_value<-p_value

# 计算富集分数
enrichment_score <- (domains.family$modinter/domains.family$mod.total)/(domains.family$totalinter/domains.family$total)
domains.family$ES<-enrichment_score

domains.family$p_adj<-p.adjust(domains.family$p_value,"fdr")
plot<-domains.family[1:5,]
plot$logP<--log10(plot$p_adj)
plot$ratio<-plot$modinter/plot$totalinter
ggplot(plot,mapping = aes(x=logP,y=ES,size=ratio*10,color=ratio,label=Name))+geom_point()+ggtheme+scale_colour_gradient(low = "#f08300",high="#b22222")+geom_text_repel(size=4,colour="black")+scale_size_continuous(limits = c(1,5))

plot<-domains.family[1:5,]
plot$logP<--log10(plot$p_adj)
plot$ratio<-plot$modinter/plot$totalinter
ggplot(plot,mapping = aes(x=logP,y=ES,size=ratio*10,color=ratio,label=Name))+geom_point()+ggtheme+scale_colour_gradient(low = "#f08300",high="#b22222")+geom_text_repel(size=4,colour="black")+scale_size_continuous(limits = c(1,5))

