rm(list = ls())
setwd("/Users/liqun/Desktop/Qun/4_TADA_Denovo/Code/")
source("./src/TADA.v.1.2.R")
tada.file = "Denovogene_Mut_TADA.txt"
tada.data=read.table(tada.file,header=T)
n.family = 549
n = data.frame(dn=n.family, ca=NA, cn=NA)
sample.counts <- list(cls1=n, cls2=n)

# dn.cls1(lof)
# dn.cls2(mis)

cls1.counts=data.frame(dn=tada.data$dn.cls1, ca=NA, cn=NA)
rownames(cls1.counts)=tada.data$gene.id
cls2.counts=data.frame(dn=tada.data$dn.cls2, ca=NA, cn=NA)
rownames(cls2.counts)=tada.data$gene.id
tada.counts=list(cls1=cls1.counts,cls2=cls2.counts)

mu=data.frame(cls1=tada.data$mut.cls1,cls2=tada.data$mut.cls2)
denovo.only=data.frame(cls1=TRUE,cls2=TRUE)

###############
# lambda=2.0 # the burden
# 501/17226 refseq Hg19 genes = 0.02908394
# pi=0.02908394
# gamma.mean.dn=(lambda-1)/pi+1
###############
# risk genes
pi=223/17226

# cls1-gamma.mean.dn-lof
# num(lof-case) = 72
# num(lof-control) = 10
# num(syn-case) = 136
# num(syn-control) = 31

cls1_λ.dn_lof = (72*31)/(10*136) 

# cls2-gamma.mean.dn-mis
# num(mis-case) = 356
# num(mis-control) = 53
cls2_λ.dn_mis = (356*31)/(53*136)

cls1_gamma.mean.dn_lof= 1+((cls1_λ.dn_lof-1)/pi)
cls2_gamma.mean.dn_mis= 1+((cls2_λ.dn_mis-1)/pi)


cls1= data.frame(gamma.mean.dn=cls1_gamma.mean.dn_lof,beta.dn=0.2,gamma.mean.CC=NA,beta.CC=NA,rho1=NA,nu1=NA,rho0=NA,nu0=NA)
cls2= data.frame(gamma.mean.dn= cls2_gamma.mean.dn_mis,beta.dn=0.2,gamma.mean.CC=NA,beta.CC=NA,rho1=NA,nu1=NA,rho0=NA,nu0=NA)
hyperpar=list(cls1=cls1,cls2=cls2)

re.TADA <- do.call(cbind.data.frame, TADA(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only))
re.TADA$qval=Bayesian.FDR(re.TADA$BF.total, pi0 = 0.95)
re.TADA.null=do.call(cbind.data.frame, TADAnull(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only, nrep=1000))
re.TADA$pval=bayesFactor.pvalue(re.TADA$BF.total,re.TADA.null$BFnull.total)
hist(re.TADA$pval)

write.table(re.TADA,"re.TADA.20211208.txt",sep = "\t")

library(qqman)
pdf("qqman20211209.pdf",width = 5.416667, height = 3.125000)
data = read.table("manhattan.plot.20211209.txt",header = T)
colnames(data) = c("SNP","CHR","BF1","BF2","BFTotal","Q","P","BP")
data = data[!duplicated(data$SNP),]
manhattan(data, genomewideline = -log10(0.05/nrow(data)), annotatePval = 0.05/nrow(data),suggestiveline = F)
abline(h = -log10(0.05), col = "red")
# qq(data$P)
dev.off()
