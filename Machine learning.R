library(EFS)
library(stringr)
library(glmnet)

vect.t3<-c(sample(36:96,43),sample(97:149,37))
SEP_EFS.t3<-as.data.frame(t(Protein.matrix.Normed.re[,vect.t3]))
SEP_EFS.t3<-cbind(SEP_EFS.t3,t(Sepsis.mod.final1[match(SEP_EFS_result_Mod$Info[1:4000],Sepsis.mod.final1$Info),vect.t3+12]))
SEP_EFS.t3<-cbind(Group=c(rep(0,43),rep(1,37)),SEP_EFS.t3)
colnames(SEP_EFS.t3)[2:ncol(SEP_EFS.t3)]<-str_extract(colnames(SEP_EFS.t3)[2:ncol(SEP_EFS.t3)],"(?<=\\|).*(?=\\|)")
colnames(SEP_EFS.t3)[1200:5199]<-SEP_EFS_result_Mod$Info[1:4000]
Info.t3<-colnames(SEP_EFS.t3)
colnames(SEP_EFS.t3)[1200:5199]<-paste0("Feature",1:4000)
MOD_EFS_result_t3<-ensemble_fs(SEP_EFS.t3,1,runs = 500,selection = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,FALSE))
MOD_EFS_result_t3<-clean(MOD_EFS_result_t3,Info.t3[-1])
MOD_EFS_result_t3<-MOD_EFS_result_t3[order(MOD_EFS_result_t3$SUM,decreasing = T),]

x=as.matrix(SEP_253[,2:ncol(SEP_253)])
y=as.matrix(SEP_253[,1])
f1 = glmnet(x, y, family="binomial", nlambda=200, alpha=0.7)
print(f1)
plot(f1, xvar="lambda", label=TRUE)
cvfit=cv.glmnet(x,y,nfolds = 3)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se  #0.1164763 0.1278326 0.1061288 0.08028258
l.coef2<-coef(cvfit$glmnet.fit,s=0.02885206,exact = F)
l.coef1<-coef(cvfit$glmnet.fit,s= cvfit$lambda.1se,exact = F)
l.coef1

coxph(Surv(time, status) ~ `Carbamyl@ALBU_438K` , data = SEP_retained.f)
coxph(Surv(time, status) ~ `Methyl@ACTB_73H` , data = SEP_retained.f)
coxph(Surv(time, status) ~ `sp|P04004|VTNC_HUMAN` , data = SEP_retained.f)
coxph(Surv(time, status) ~ `sp|P48304|REG1B_HUMAN` , data = SEP_retained.f)

cox_model1 <- coxph(Surv(time, status) ~ `sp|P04004|VTNC_HUMAN` + `sp|P48304|REG1B_HUMAN` + `Methyl@ACTB_73H` , data = SEP_retained.f);concordance(cox_model1)