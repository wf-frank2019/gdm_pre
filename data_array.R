#!/usr/local/bin/Rscript
library(data.table)
library(optparse)
option_list = list(
         make_option(c("-p","--prefix"), type = "character", default = "GDM", help = "prefix of data"),
         make_option(c("-i", "--infoscore"), type = "numeric", default = 0.8, help = "imputation score R2"),
         make_option(c("-o","--object"), type = "character", default = "T2D", help = "object you want to study"),
	 make_option(c("-n", "--pcanum"), type = "numeric", default = 10, help = "PCA number"),
	 make_option(c("-c","--control"), type = "numeric", help = "flowchat control,0-10,order irreversible")
)
opt = parse_args(OptionParser(option_list = option_list))
pre = opt$prefix
infos = opt$infoscore
pheno = opt$object
pcan = opt$pcanum
cont = opt$control

if(cont==0){
print(paste0("get final base data Start at : ",date()))
options(scipen=200)
zt = fread(paste0(pheno,".QC1"),header=T)
zf = fread(paste0(pheno,".QC2"),header=T)
zt = as.data.frame(zt)
zf = as.data.frame(zf)
zf$BETA = zf$BETA*(-1)
z = rbind(zt,zf)
z$OR = exp(z$BETA)
colnames(z) = c("CHR","BP","SNP","A1","A2","N","SE","P","BETA","INFO","MAF","OR")
fwrite(z,paste0(pheno,".QC"),sep="\t",quote=F,row.names=F,col.names=T)
print(paste0("get final base data Finish at : ",date()))

if(FALSE){
dat = read.table(gzfile(paste0(pheno,".38.QC1.gz")), header=T)
for (i in 1:nrow(dat))
{
dat[i,3] = paste0(dat[i,1],":",dat[i,2],":",dat[i,5],":",dat[i,4])
}
dat2 = dat
dat2$OR = 1/dat2$OR
for (i in 1:nrow(dat2))
{
dat2[i,3] = paste0(dat2[i,1],":",dat2[i,2],":",dat2[i,4],":",dat2[i,5])
}
data = rbind(dat,dat2)
write.table(data,paste0(pheno,".QC"),sep="\t",quote=F,row.names=F)
}
}

if(cont==1){
print(paste0("heterozygosity check Start at : ",date()))
data = read.table(paste0(pre,".",infos,".QC.het"),header=T)
valid = subset(data,F <= mean(data$F)+3*sd(data$F)  & F >= mean(data$F)-3*sd(data$F) )
invalid = subset(data,F > mean(data$F)+3*sd(data$F)  |  F < mean(data$F)-3*sd(data$F) )
invalid$Heterozygosity = ifelse(invalid$F>0, "Too high","Too low")
write.table(valid[,c(1,2)], paste0(pre,".",infos,".valid.sample"), quote=FALSE, row.names=FALSE) 
write.table(invalid[,c(1,2,ncol(invalid))], paste0(pre,".",infos,".invalid.sample"), quote=FALSE, row.names=FALSE)
print(paste0("heterozygosity check Finish at : ",date()))
}

if(cont==2){
print(paste0("mismatch identify Start at : ",date()))
bim = fread(paste0(pre,".bim"))
bim = as.data.frame(bim)
print("read bim done")
colnames(bim) = c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
qc = read.table(paste0(pre,".",infos,".QC.snplist"), header = F, stringsAsFactors = F)
print("read QC.snplist done")
#height = read.table(gzfile(paste0(pheno,".QC.gz")),header = T,stringsAsFactors = FALSE, sep="\t")
height = fread(paste0(pheno,".QC.gz"))
height = as.data.frame(height)
print("read QCed-base done")
height$A1 = toupper(height$A1);height$A2 = toupper(height$A2);bim$B.A1 = toupper(bim$B.A1);bim$B.A2 = toupper(bim$B.A2)

info = merge(bim, height, by = c("SNP", "CHR", "BP"))
info = info[info$SNP %in% qc$V1,]
complement = function(x) {switch(x,"A" = "T","C" = "G","T" = "A","G" = "C",return(NA))}
info.match = subset(info, A1 == B.A1 & A2 == B.A2)
info$C.A1 = sapply(info$B.A1, complement);info$C.A2 <- sapply(info$B.A2, complement)
info.complement = subset(info, A1 == C.A1 & A2 == C.A2)
complement.snps = bim$SNP %in% info.complement$SNP
bim[complement.snps,]$B.A1 = sapply(bim[complement.snps,]$B.A1, complement)
bim[complement.snps,]$B.A2 = sapply(bim[complement.snps,]$B.A2, complement)

info.recode = subset(info, A1 == B.A2 & A2 == B.A1)
recode.snps = bim$SNP %in% info.recode$SNP
tmp = bim[recode.snps,]$B.A1
bim[recode.snps,]$B.A1 = bim[recode.snps,]$B.A2 
bim[recode.snps,]$B.A2 = tmp
info.crecode = subset(info, A1 == C.A2 & A2 == C.A1)
com.snps = bim$SNP %in% info.crecode$SNP
tmp = bim[com.snps,]$B.A1
bim[com.snps,]$B.A1 = as.character(sapply(bim[com.snps,]$B.A2, complement))
bim[com.snps,]$B.A2 = as.character(sapply(tmp, complement))
write.table(bim[,c("SNP", "B.A1")],paste0(pre,".",infos,".a1"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
mismatch = bim$SNP[!(bim$SNP %in% info.match$SNP | bim$SNP %in% info.complement$SNP | bim$SNP %in% info.recode$SNP | bim$SNP %in% info.crecode$SNP)] 
write.table(mismatch,paste0(pre,".",infos,".mismatch"),quote = FALSE,row.names = FALSE,col.names = FALSE)
print(paste0("mismatch identify Finish at : ",date()))
}

if(cont==3){
print(paste0("sex check Start at : ",date()))
valid = read.table(paste0(pre,".",infos,".valid.sample"), header=T)
dat = read.table(paste0(pre,".",infos,".QC.sexcheck"), header=T)
valid = subset(dat, STATUS=="OK" & FID %in% valid$FID)
invalid = subset(dat, STATUS=="PROBLEM" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], paste0(pre,".",infos,".QC.valid"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE) 
write.table(invalid[,c("FID", "IID")], paste0(pre,".",infos,".QC.Sexinvalid"), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
print(paste0("sex check Finish at : ",date()))
}

if(cont==4){
#If your quantitative traits only contain OR value:
dat = read.table(gzfile(paste0(pheno,".QC.gz")), header=T)
dat$BETA = log(dat$OR)
fwrite(dat,paste0(pheno,".QC.Transformed"), quote=F, row.names=F)
}

if(cont==5){
#make Plink clumping interval
a = seq(from = 0, to = 0.5, by = 0.01)
a[1] = 0.001
tmp = matrix(nrow=length(a),ncol=3)
tmp[,1] = a;tmp[,2] = 0;tmp[,3] = a
write.table(tmp,"range_list",sep=" ",quote=F,row.names=F,col.names=F)
}

if(cont==6){
#for binary traits results
print(paste0("plink best fit Start at : ",date()))
library(ggpubr)
library(ggpmisc)
library(ggplot2)
p.threshold = seq(from = 0, to = 0.5, by = 0.01);p.threshold[1]=0.001
phenotype = read.table(paste0(pre,".pheno"), header=T)
colnames(phenotype) = c("FID","IID","Phenotype")
pcs = read.table(paste0(pre,".",infos,".eigenvec"), header=F)
colnames(pcs) = c("FID", "IID", paste0("PC",1:pcan))
covariate = read.table(paste0(pre,".cov"), header=T)
phenos = merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))#merge the files, first height+sex；then +PCAs
#calculate the null model (model with PRS) using a logistic regression  (as is Qualitative)
null.model = glm(Phenotype~.,family = binomial(),data=phenos[,!colnames(phenos)%in%c("FID","IID")])
null.model1 = glm(Phenotype~1,family = binomial(),data=phenos[,!colnames(phenos)%in%c("FID","IID")])
null.r2 = as.numeric(1-logLik(null.model)/logLik(null.model1))
null.aic = summary(null.model)$aic 

plink.result = NULL
try = 0
prs.result = NULL
for(i in p.threshold){
	prs = read.table(paste0(pre,".",infos,".",i,".profile"), header=T)
	pheno.prs =  merge(phenos, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
	model = glm(Phenotype~.,family = binomial(),data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
	model1 = glm(Phenotype~1,family = binomial(),data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
	model.r2 = as.numeric(1-logLik(model)/logLik(model1))	
	model.aic = summary(model)$aic
        p.chisq = anova(object = model,test = "Chisq")["SCORE",5]
	prs.r2 = model.r2 - null.r2
	prs.result = rbind(prs.result, data.frame(Threshold=i, McFaddenR2=prs.r2, AIC=model.aic, p=p.chisq))
	try = try + 1
	if(try == 1){
	plink.result = as.data.frame(matrix(nrow=nrow(prs),ncol=length(p.threshold)+3))
	plink.result[,1:3] = prs[,1:3]
	plink.result[,4] = prs[,6]
	colnames(plink.result) = c("FID","IID","PEHNO",paste0("SCORE",i))
	}
	else{
	plink.result[,try+3] = prs[,6] 
	colnames(plink.result)[try+3] = paste0("SCORE",i)}
}
max = prs.result[which.max(prs.result$McFaddenR2),1]
print(paste0("best p value: ",max))
best = read.table(paste0(pre,".",infos,".",max,".profile"), header=T)
best$z.SCORE = (best$SCORE-mean(best$SCORE))/sd(best$SCORE)
write.table(prs.result,"plinkprs.result",sep="\t",quote=FALSE,row.names=FALSE)
write.table(best,"plinkprs.best",row.names=FALSE,quote=FALSE,sep="\t")
write.table(plink.result,"plinkprs.all_score",sep="\t",quote=FALSE,row.names=FALSE)
print(paste0("plink best fit Finish at : ",date()))

plinkbest = merge(best,phenotype,by=c("FID","IID"))
prsplink1 = ggplot(data=plinkbest,aes(x=z.SCORE))+geom_histogram(fill="#FFA07A",color="#696969",bins=20) +
             labs(x="Polygenic Risk Score",y="Number of individuals")+theme_bw()
ggsave("PRSdistrib_plink.pdf",prsplink1,width=10,height=8,limitsize = FALSE)
prsplink2 = ggplot(plinkbest,aes(x=as.factor(Phenotype),y=z.SCORE,fill=as.factor(Phenotype))) +geom_boxplot() + 
             scale_fill_brewer(palette="Dark2")+theme_classic() +
             geom_signif(comparison=list(c("0","1")),step_increase = 0.1,map_signif_level = F,test = wilcox.test) +
             labs(x="GDM class",y="Polygenic Risk Score") + theme(legend.position="none")
ggsave("PRS_plink_boxplot.pdf",prsplink2,width=10,height=8,limitsize = FALSE)
}

if(cont==7){
#for quantitative traits results
print(paste0("plink best fit Start at : ",date()))
library(ggpubr)
library(ggpmisc)
library(ggplot2)
p.threshold = seq(from = 0, to = 0.5, by = 0.01);p.threshold[1]=0.001
phenotype = read.table(paste0(pre,".pheno"), header=T)
colnames(phenotype) = c("FID","IID","Phenotype")
pcs = read.table(paste0(pre,".",infos,".eigenvec"), header=F)
colnames(pcs) = c("FID", "IID", paste0("PC",1:pcan))
covariate = read.table(paste0(pre,".cov"), header=T)
phenos = merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
#calculate the null model (model with PRS) using a lm regression  (as is Quantitative)
null.model = lm(Phenotype~.,data=phenos[,!colnames(phenos)%in%c("FID","IID")])
null.r2 = summary(null.model)$r.squared

plink.result = NULL
try = 0
prs.result = NULL
for(i in p.threshold){
   prs = read.table(paste0(pre,".",infos,".",i,".profile"), header=T)
   pheno.prs =  merge(phenos, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
   model = lm(Phenotype~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
   model.r2 = summary(model)$r.squared
   prs.r2 = model.r2-null.r2
   #Adding the PRS score on the basis of the original regression fitting, how much can the regression effect be improved
   prs.coef = summary(model)$coeff["SCORE",]
   prs.beta = as.numeric(prs.coef[1])
   prs.se = as.numeric(prs.coef[2])
   prs.p = as.numeric(prs.coef[4])
   prs.result = rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
   try = try + 1
   if(try == 1){
     plink.result = as.data.frame(matrix(nrow=nrow(prs),ncol=length(p.threshold)+3))
     plink.result[,1:3] = prs[,1:3]
     plink.result[,4] = prs[,6]
     colnames(plink.result) = c("FID","IID","PEHNO",paste0("SCORE",i))
     }
   else{
     plink.result[,try+3] = prs[,6] 
     colnames(plink.result)[try+3] = paste0("SCORE",i)}
}
max = prs.result[which.max(prs.result$R2),1]
print(paste0("best p value: ",max))
best = read.table(paste0(pre,".",infos,".",max,".profile"), header=T)
best$z.SCORE = (best$SCORE-mean(best$SCORE))/sd(best$SCORE)
write.table(prs.result,"plinkprs.result",sep="\t",quote=FALSE,row.names=FALSE)
write.table(best,"plinkprs.best",row.names=FALSE,quote=FALSE,sep="\t")
write.table(plink.result,"plinkprs.all_score",sep="\t",quote=FALSE,row.names=FALSE)
print(paste0("plink best fit Finish at : ",date()))

plinkbest = merge(best,phenotype,by=c("FID","IID"))
prsplink1 = ggplot(data=plinkbest,aes(x=z.SCORE))+geom_histogram(fill="#FFA07A",color="#696969",bins=20) +
             labs(x="Polygenic Risk Score",y="Number of individuals")+theme_bw()
ggsave("PRSdistrib_plink.pdf",prsplink1,width=10,height=8,limitsize = FALSE)

getLabel <- function(dat,group=10){
  dat = dat[order(dat["z.SCORE"]),]
  dat$label = 0 
 for(i in 1:group){
     if(i!=group){
     dat$label[((i-1)*round(nrow(dat)/group)+1):(i*round(nrow(dat)/group))] = i }
     else{
     dat$label[((i-1)*round(nrow(dat)/group)+1):(nrow(dat))] = i }}
  return(dat)
  }
plinkbest = getLabel(plinkbest)
prsplink2 = ggplot(plinkbest,aes(x=as.factor(label),y=Phenotype)) + geom_boxplot()+
	    theme_bw() + theme(legend.title=element_blank()) +
	    theme(panel.grid = element_blank()) + labs(x="PRS ranked group",y="phenotype value")
ggsave("PRS_plink_boxplot.pdf",prsplink2,width=10,height=8,limitsize = FALSE)

}

if(cont==8){
#get all covariate for Regression
covariate = read.table(paste0(pre,".cov"), header=T)
pcs = read.table(paste0(pre,".",infos,".eigenvec"), header=F)
colnames(pcs) = c("FID","IID", paste0("PC",1:pcan))
cov = merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,paste0(pre,".",infos,".covariate"), quote=FALSE, row.names=FALSE)
}

if(cont==9){
#prsice2 output best score
prsbest = read.table(paste0(pre,".",infos,".best"), header=T) 
phenotype = read.table(paste0(pre,".pheno"),header=T)
prsbest$z.PRS = (prsbest$PRS-mean(prsbest$PRS))/sd(prsbest$PRS)
best = merge(prsbest,phenotype,by=c("FID","IID"))
write.table(best,"prsice2prs.best",row.names=FALSE,quote=FALSE,sep="\t")
}

if(cont==10){
print(paste0("PRSice2 best fit Start at : ",date()))
library(pROC)
library(ggpubr)
library(ggpmisc)
library(ggplot2)
phenotype = read.table(paste0(pre,".pheno"),header=T)
colnames(phenotype) = c("FID","IID","Phenotype")
pcs = read.table(paste0(pre,".",infos,".eigenvec"), header=F)
colnames(pcs) = c("FID", "IID", paste0("PC",1:pcan))
covariate = read.table(paste0(pre,".cov"), header=T)
phenos = merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))

null.model = glm(Phenotype~.,family = binomial(),data=phenos[,!colnames(phenos)%in%c("FID","IID")])
null.model1 = glm(Phenotype~1,family = binomial(),data=phenos[,!colnames(phenos)%in%c("FID","IID")])
null.r2 = as.numeric(1-logLik(null.model)/logLik(null.model1))

score = fread(paste0(pre,".",infos,".all_score"), header=T)
score = as.data.frame(score)
p.threshold = ncol(score)-2
prs.result =  NULL
for(i in 1:p.threshold){
	score1 = score[,c(1,2,i+2)]
	str = unlist(strsplit(colnames(score1)[3],split = "Pt_"))[2]
	colnames(score1)[3] = "SCORE"
	pheno.prs = merge(phenos, score1, by=c("FID", "IID"))
	model = glm(Phenotype~.,family = binomial(),data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
	model1 = glm(Phenotype~1,family = binomial(),data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
	model.r2 = as.numeric(1-logLik(model)/logLik(model1))
	prs.r2 = model.r2 - null.r2
	model.aic = summary(model)$aic
	p.chisq = anova(object = model,test = "Chisq")["SCORE",5]
	logist_pre = predict(model,pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")],type='response')	
	logist_roc = roc(pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")]$Phenotype,as.numeric(logist_pre),levels = c(0,1),direction='<')
	auc = as.numeric(logist_roc$auc)
	prs.result = rbind(prs.result, data.frame(Threshold=str, McFaddenR2=prs.r2, AUC=auc,AIC=model.aic, Pvalue=p.chisq, count=i))
}
max = prs.result[which.max(prs.result$AUC),1]
print(paste0("best p value: ",max))
count = prs.result[which.max(prs.result$AUC),c("count")]+2
score1 = score[,c(1,2,count)]
colnames(score1)[3] = "SCORE"
pheno.prs = merge(phenos, score1, by=c("FID", "IID"))
model = glm(Phenotype~.,family = binomial(),data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
logist_pre = predict(model,pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")],type='response')
logist_roc = roc(pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")]$Phenotype,as.numeric(logist_pre),levels = c(0,1),direction='<')
print(paste0("AUC of clumping model is: ",as.numeric(logist_roc$auc)))
plot.roc(logist_roc,col="#FF0000",smooth=TRUE,print.auc=TRUE,cex.lab=1.5,percent=TRUE,print.auc.cex=2.5,legacy.axes=TRUE)

pheno.prs$z.SCORE = (pheno.prs$SCORE-mean(pheno.prs$SCORE))/sd(pheno.prs$SCORE)
prsicebest = pheno.prs[,c(1,2,3,ncol(pheno.prs)-1,ncol(pheno.prs))]
write.table(prsicebest,"prsice2PRS.best",row.names=FALSE,quote=FALSE,sep="\t")
prs.result = prs.result[,-ncol(prs.result)]
write.table(prs.result,"prsice2PRS.result",sep="\t",quote=FALSE,row.names=FALSE)

log_list=list()
for(i in 1:1000){
  train_sub = sample(nrow(pheno.prs),9/10*nrow(pheno.prs))
  train_data = pheno.prs[train_sub,]
  test_data = pheno.prs[-train_sub,]
  train_data$Phenotype = as.factor(train_data$Phenotype)
  test_data$Phenotype = as.factor(test_data$Phenotype)
  model = glm(Phenotype~.,family = binomial(),data=train_data[,!colnames(train_data)%in%c("FID","IID")])
  logist_pre = predict(model,test_data[,!colnames(test_data)%in%c("FID","IID")],type='response')
  logist_roc = roc(test_data$Phenotype,as.numeric(logist_pre),levels = c(0,1),direction='<')
  log_list[[i]] = as.numeric(logist_roc$auc)
}
log_list_numeric = unlist(log_list)
print(paste0("AUC of clumping model(leave-one-out,1000 random sampling) is: ",mean(log_list_numeric)))
print(paste0("PRSice2 best fit Finish at : ",date()))

if(FALSE){
prs = read.table("prsice2prs.best", header=T)
pheno.prs = merge(pheno, prs[,c("FID","IID", "PRS")], by=c("FID", "IID"))
model = glm(Phenotype~.,family = binomial(),data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
logist_pre = predict(model,pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")],type='response')
logist_roc = roc(pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")]$Phenotype,as.numeric(logist_pre),levels = c(0,1),direction='<')
print(paste0("AUC of clumping model is: ",as.numeric(logist_roc$auc)))
plot.roc(logist_roc,col="#FF0000",smooth=T,print.auc=T,cex.lab=1.5,percent=T,print.auc.cex=2.5,legacy.axes=TRUE)
#logist2 = as.data.frame(t(model$fitted.values))
#plot.roc(pheno.prs$Phenotype,as.numeric(logist2),col="red",print.auc=T,cex.lab=1.5,percent=T,print.auc.cex=2.5,legacy.axes=TRUE)
log_list=list()
for(i in 1:1000){
  train_sub = sample(nrow(pheno.prs),9/10*nrow(pheno.prs))
  train_data = pheno.prs[train_sub,]
  test_data = pheno.prs[-train_sub,]
  train_data$Phenotype = as.factor(train_data$Phenotype)
  test_data$Phenotype = as.factor(test_data$Phenotype)
  model = glm(Phenotype~.,family = binomial(),data=train_data[,!colnames(train_data)%in%c("FID","IID")])
  logist_pre = predict(model,test_data[,!colnames(test_data)%in%c("FID","IID")],type='response')
  logist_roc = roc(test_data$Phenotype,as.numeric(logist_pre),levels = c(0,1),direction='<')
  #Bayes_obs = data.frame(prob=Bayes_pre ,obs=test_data$type)
  #table(test_data$type,Bayes.pred ,dnn = c('real','predict'))
  log_list[[i]] = as.numeric(logist_roc$auc)
}
log_list_numeric = unlist(log_list)
print(paste0("AUC of clumping model(leave-one-out,random 1000) is: ",mean(log_list_numeric)))
}

print(paste0("plot Start at : ",date()))
prsPlot <- function(data,
		    group1 = 10,
		    group2 = 5,
		    group3 = 20,
		    widths = 10,
		    heights = 8,
		    dpis = 300){
  if(widths > 30 || heights > 25){
    stop("Please reduce your picture size")
  }
  if(widths > 1.5*heights || heights > 1.5*widths){
    stop("Please balance your picture size")
  }
  if(dpis > 900){
    stop("dpi must set inside (0,900]")
  }
  if(group1*group2*group3*widths*heights*dpis<=0){
    stop("Please input correct positive integer!")
  }
  if(group2%%2 == 0){
    stop("Please input odd number of group2!")
  }
theme_sam = theme_classic()+theme(axis.title=element_text(face="bold", size=18),
                              axis.text=element_text(size=14),
                              legend.title=element_text(face="bold", size=18),
                              legend.text=element_text(size=14),
                              axis.text.x=element_text(angle=45, hjust=1))
theme_frank = theme_bw()+theme(text = element_text(size=20),axis.title.x=element_text(size=20),
		axis.title.y=element_text(size=20))+theme(panel.grid = element_blank())
print(paste0("random seed: ",sample(1:100000000, 1)))
prsprsice1 = ggplot(data,aes(x=z.SCORE))+geom_histogram(fill="skyblue",color="#696969",bins=20)+
		labs(x="Polygenic Risk Score",y="Number of individuals") + theme_frank
ggsave("PRSdistrib_prsice.pdf",prsprsice1,width = widths,height = heights,dpi=dpis,limitsize = FALSE)
print("People risk score distribution plot done")
prsprsice2 = ggplot(data,aes(x=as.factor(Phenotype),y=z.SCORE,fill=as.factor(Phenotype))) +
		geom_boxplot()+scale_fill_brewer(palette="Dark2")+labs(x="Phenotype class",y="Polygenic Risk Score") +
		geom_signif(comparison=list(c("0","1")),step_increase = 0.1,map_signif_level = F,test = wilcox.test) +
		theme_frank + theme(legend.position="none")
ggsave("PRS_prsice_boxplot.pdf",prsprsice2,width = widths,height = heights,dpi=dpis,limitsize = FALSE)
print("PRS comparision among phenotypes done")
data2 = data[data$Phenotype == 1,]
data3 = data[data$Phenotype == 0,]

getLabel <- function(dat,group){
  dat = dat[order(dat["z.SCORE"]),]
  dat$label = 0 
 for(i in 1:group){
     if(i!=group){
     dat$label[((i-1)*round(nrow(dat)/group)+1):(i*round(nrow(dat)/group))] = i }
     else{
     dat$label[((i-1)*round(nrow(dat)/group)+1):(nrow(dat))] = i }}
  return(dat)
  }
data2 = getLabel(data2,group=group1); data3 = getLabel(data3,group=group1)
plot = rbind(data2,data3)
plot$Phenotype = as.factor(plot$Phenotype)
p1 = ggplot(plot,aes(x=as.factor(label),y=z.SCORE)) + geom_boxplot(aes(fill=Phenotype)) + 
	scale_fill_manual(values = c("darkgreen", "#CD5555")) + theme_bw() + 
	labs(x="PRS percentile",y="SCORE")
plotname1 = paste0(group1,"percentile.pdf")
ggsave(plotname1,p1,width = widths,height= heights,dpi=dpis,limitsize = FALSE)
print("group percentile plot done")
dat2 = data[order(data["z.SCORE"]),]
dat2 = getLabel(dat2,group=group2)
tmp1 = list()
odbase = ceiling(group2/2)
 for(i in 1:group2){
  targetP = dim(dat2[dat2$label == i & dat2$Phenotype == 1,])[1]
  baseN = dim(dat2[dat2$label == odbase & dat2$Phenotype == 0,])[1]
  targetN = dim(dat2[dat2$label == i & dat2$Phenotype == 0,])[1]
  baseP = dim(dat2[dat2$label == odbase & dat2$Phenotype == 1,])[1]
  or = (targetP*baseN)/(targetN*baseP)
  tmp1 = rbind(tmp1,data.frame(strata = i,OR = or)) 
}
yn = paste0("Odds Ratio(based on group",odbase,")")
p2 = ggplot(tmp1,aes(x=strata,y=OR,ymin=0,ymax=3))+geom_point(colour = "#D55E00",size = 4)+
	labs(x="PRS rank group",y=yn)+geom_pointrange(colour = "#D55E00",size = 0.5)+theme_frank
plotname2 = paste0(group2,"OddsRatio.pdf")
ggsave(plotname2,p2,width = widths,height = heights,dpi=dpis,limitsize = FALSE)
print("group Odds Ratio plot done") 
dat3 = data[order(data["z.SCORE"]),]
dat3 = getLabel(dat3,group=group3)
tmp2 = list()
 for(i in 1:group3){
    case = dim(dat3[dat3$label == i & dat3$Phenotype == 1,])[1]
    contr = dim(dat3[dat3$label == i & dat3$Phenotype == 0,])[1]
    frac = (case*dim(dat3[dat3$Phenotype == 0,])[1])/(contr*dim(dat3[dat3$Phenotype == 1,])[1])
    tmp2 = rbind(tmp2,data.frame(class = i,fraction = frac)) }
p3 = ggplot(tmp2,aes(x=as.factor(class),y=fraction))+geom_point(shape = 16,size = 8)+
	labs(x = "PRS rank group", y = "case/control fraction") + theme_frank 
plotname3 = paste0(group3,"caseproportion.pdf")
ggsave(plotname3,p3,width = widths,height = heights,dpi=dpis,limitsize = FALSE)
print("case proportion in each group plot done")
}

prsPlot(prsicebest) #output all plots of results
print(paste0("plot Finish at : ",date()))

}

if(cont==11){
library(qqman)
library(dplyr)
library(ggplot2)
library(ggrepel)
theme_wf = theme_bw() + theme(legend.position="none",panel.border=element_blank(),
		panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
		axis.line=element_line(colour="black"))
p1 = 5e-8
clump = read.table(paste0(pre,".",infos,".snp"),header=T)
gwas = fread("base.plot",header=T,stringsAsFactors=F)
colnames(gwas) = c("SNP","CHR","BP","P")
gwas = as.data.frame(gwas)
chro = seq(1,22,by=1)
gwas_auto = gwas[gwas$CHR %in% chro,]
gwas_auto$CHR = as.numeric(gwas_auto$CHR)
sig = gwas_auto[gwas_auto$P < 1e-6,]
nosig = gwas_auto[gwas_auto$P >= 1e-6,]
for(i in 1:22){
 samples = nosig[nosig$CHR==i,]
 samples_index = samples[sample(nrow(samples),nrow(samples)/10,replace=F),]
 sig = rbind(sig,samples_index)
}
sig$BP =  as.numeric(sig$BP)
sig$is_highlight = ifelse(sig$SNP %in% clump$SNP,"yes","no")
sig$is_annotate=ifelse(sig$is_highlight=="yes" & sig$P < 1e-6,"yes","no")
wf <- sig %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>% left_join(sig, ., by=c("CHR"="CHR")) %>%arrange(CHR, BP) %>% mutate( BPcum=BP+tot)
axisdf = wf %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pa = ggplot(wf, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
	scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
	scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) + scale_y_continuous(expand = c(0, 0)) +
	geom_point(data=subset(wf, is_highlight=="yes"), color="orange", size=2) +
	geom_label_repel(data=subset(wf, is_annotate=="yes"), aes(label=SNP), size=2) + theme_wf +
	labs(x="chromosome")+geom_hline(aes(yintercept= -log10(p1)), colour="#990000", linetype="dashed")
ggsave("snpsAfterClumpHL.pdf",pa,width=16,height=10,limitsize = FALSE)
pb = ggplot(wf, aes(x=BPcum, y=-log10(P))) + geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
	scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
	scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) + scale_y_continuous(expand = c(0, 0)) +
	geom_point(data=subset(wf, is_highlight=="yes"), color="#FFD9BB", size=2) + theme_wf +
	labs(x="chromosome")+geom_hline(aes(yintercept= -log10(p1)), colour="#990000", linetype="dashed")
ggsave("snpsAfterClump.pdf",pb,width=16,height=10,limitsize = FALSE)
}


if(FALSE){
data = read.table("/home/wangfan819/PGD/BRI002/gdmPRS/prsTest/Final3/GDM.0.8.eigenvec",header=F)
pheno = read.table("/home/wangfan819/PGD/BRI002/gdmPRS/prsTest/Final3/GDM.pheno",header=T)
colnames(data)[1:2]= c("FID","IID")
merge = merge(data,pheno,by=c("FID","IID"))
merge$T2D = as.factor(merge$T2D)
tmp = list()
for(i in 1:9){
dat = merge[,c(3,i+3,23)]
colnames(dat) = c("X","Y","GDM")
yn = paste0("pc",i+1)
tmp[[i]] = ggplot(dat,aes(x=X,y=Y,color=GDM))+geom_point()+theme_bw()+labs(x="pc1",y=yn)
}
tmp2 = list()
for(i in 1:8){
dat = merge[,c(4,i+4,23)]
colnames(dat) = c("X","Y","GDM")
yn = paste0("pc",i+2)
tmp2[[i]] = ggplot(dat,aes(x=X,y=Y,color=GDM))+geom_point()+theme_bw()+labs(x="pc2",y=yn)
}
}

