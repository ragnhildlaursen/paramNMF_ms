rm(list=ls())
setwd("~/projects/paramNMF_ms")
##############################################################
## Plot robustness of signatures
##############################################################


library(ggplot2)
library(ggpubr)


#load("BRCA/BRCArun/SimilaritySignaturesBRCA119.RData")
#load("BRCA/BRCArun/SimilaritySignaturesBRCA21boot.RData")
load("BRCA/result/SimilaritySignaturesBRCA214v2.RData")
noSig = nrow(ResCosineMatTri)
# Transforming data
datTri = data.frame(s = t(ResCosineMatTri[,c(1:12)]), type = "tri")
#m = cosMatch(TriRes$Signatures, MixRes$Signatures)$match
datMix = data.frame(s = t(ResCosineMatMix[,c(1:12)]), type = "mix")
#m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
datDi = data.frame(s = t(ResCosineMatDi[,c(1:12)]), type = "di")
#m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
datMono = data.frame(s = t(ResCosineMatMono[,c(1:12)]), type = "mono")
dat1 = rbind(datTri,datDi,datMix,datMono)
dat2 = reshape(dat1, varying = colnames(dat1)[1:noSig], direction = "long")
dat2$type = factor(dat2$type, levels = c("mono","di","mix","tri"))
# remove mix
#dat2 = dat2[dat2$type != "mix",]

##------------------------------------------------------------------
## Compare true signatures and estimated signatures
##------------------------------------------------------------------
par(mfrow = c(2,2))
plots = list()
for(sig in 1:noSig){
  plots[[sig]] = ggplot(dat2[dat2$time == sig,], aes(x=type, y=s, fill=type)) +
    geom_violin(position=position_dodge(1), scale = "width") +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    ggtitle(paste("Signature",sig))
}

ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],
          plots[[5]],plots[[6]],plots[[7]],plots[[8]],
          nrow = 4, ncol = 2, common.legend = TRUE)
ggplot(dat2, aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1), scale = "width") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  ggtitle("Reconstruction with bootstrap BRCA214")

ggplot(dat2, aes(x=type, y=s, col = type)) +
  geom_boxplot()+
  geom_jitter(size = 0.2) +
  
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  ggtitle("Reconstruction with bootstrap BRCA214")




load("DiRes.RData")

m0 = cosMatch(DiRes$Signatures, ResultsDi[[1]]$Signatures)$match
m1 = cosMatch(DiRes$Signatures, ResultsDi[[2]]$Signatures)$match
m2 = cosMatch(DiRes$Signatures, ResultsDi[[3]]$Signatures)$match
m3 = cosMatch(DiRes$Signatures, ResultsDi[[8]]$Signatures)$match

par(mfrow = c(5,1), mar = c(2,2,3,3))
for(i in 1:4){
  plot(DiRes$Signatures[i,], type = "h", main = DiRes$gkl)
  plot(ResultsDi[[1]]$Signatures[m0[i],], type = "h", main = ResultsDi[[1]]$gkl)
  plot(ResultsDi[[2]]$Signatures[m1[i],], type = "h", main = ResultsDi[[2]]$gkl)
  plot(ResultsDi[[3]]$Signatures[m2[i],], type = "h", main = ResultsDi[[3]]$gkl)
  plot(ResultsDi[[8]]$Signatures[m3[i],], type = "h", main = ResultsDi[[8]]$gkl)
}
