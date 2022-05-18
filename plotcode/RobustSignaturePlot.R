rm(list=ls())
setwd("~/projects/paramNMF_ms")
##############################################################
## Plot robustness of signatures
##############################################################
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
##################
## cosine similarity 
load("BRCA/result/SimilarityBRCA214w4sig500iparamboot50.RData")
load("BRCA/result/SimilarityBRCA214w4sig500iparamboot50Mix.RData")
#load("UCUT/SimilaritySignaturesUCUTparamboot.RData")
order21 = c(1,4,3,2)
order214 = c(4,2,3,1)
order = order214
noSig = nrow(ResCosineMatMono)
# Transforming data
#datPenta = data.frame(s = t(ResCosineMatPenta[order,]), type = "Penta")
datTri = data.frame(s = t(ResCosineMatTri[order,]), type = "Tri")
#m = cosMatch(TriRes$Signatures, MixRes$Signatures)$match
datMix = data.frame(s = t(ResCosineMatMix[order,]), type = "Mix")
#m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
datDi = data.frame(s = t(ResCosineMatDi[order,]), type = "Di")
#m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
datMono = data.frame(s = t(ResCosineMatMono[order,]), type = "Mono")
dat1 = rbind(datTri,datDi,datMix,datMono)
#dat1 = rbind(datPenta,datDi,datMono)
dat2 = reshape(dat1, varying = colnames(dat1)[1:noSig], direction = "long")
dat2$type = factor(dat2$type, levels = c("Mono","Di","Mix","Tri"))
#dat2$type = factor(dat2$type, levels = c("Mono","Di","Penta"))
##------------------------------------------------------------------
## Compare true signatures and estimated signatures
##------------------------------------------------------------------
colors3 = c("#AC0136", "#FF9912","#A895CD", "#27408B")
colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")

#colors3 = brewer.pal(n = 6, name = "Dark2")[c(1,3,6)]

# plot in one plot
sig.labs <- paste("Signature",1:noSig)
names(sig.labs) <- c(1:noSig)

ggplot(dat2, aes(x=type, y=s, fill=type, col = type)) +
  #geom_violin(position=position_dodge(1), scale = "width") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  #geom_jitter()+
  geom_boxplot()+
  facet_wrap(time~., nrow = 1, scales = "free_y", labeller = labeller(time = sig.labs))+
  scale_color_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                     name = "Interaction model")+
  scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")+
  xlab("Interaction model")+
  ylab("Cosine similarity")+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"))




par(mfrow = c(2,2))
plots = list()
for(sig in 1:noSig){
  plots[[sig]] = ggplot(dat2[dat2$time == sig,], aes(x=type, y=s, fill=type, col = type)) +
    #geom_violin(position=position_dodge(1), scale = "width") +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    boxplot()
    #geom_jitter()+
    ggtitle(paste("Signature",sig))+
    scale_color_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                      name = "Interaction model")+
    scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                       name = "Interaction model")+
    theme(legend.position = "bottom")
}

ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],
          #plots[[5]],plots[[6]],plots[[7]],plots[[8]],
          nrow = 2, ncol = 2, common.legend = TRUE)





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
