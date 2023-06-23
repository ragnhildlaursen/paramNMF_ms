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
load("BRCA/result/SimilaritySignaturesBRCA214optimal8v2.RData")

order = c(7,8,3,2,5,1,4,6)
noSig = nrow(ResCosineMatMono)
# Transforming data
datTri = data.frame(s = t(ResCosineMatTri[order,]), type = "Tri")
#m = cosMatch(TriRes$Signatures, MixRes$Signatures)$match
datMix = data.frame(s = t(ResCosineMatMix[order,]), type = "Mix")
#m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
datDi = data.frame(s = t(ResCosineMatDi[order,]), type = "Di")
#m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
datMono = data.frame(s = t(ResCosineMatMono[order,]), type = "Mono")
dat1 = rbind(datTri,datDi,datMix,datMono)
dat2 = reshape(dat1, varying = colnames(dat1)[1:noSig], direction = "long")
dat2$type = factor(dat2$type, levels = c("Mono","Di","Mix","Tri"))
##------------------------------------------------------------------
## Compare true signatures and estimated signatures
##------------------------------------------------------------------
#colors3 = c("#AC0136", "#FF9912","#A895CD", "#27408B")
colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")


# plot in one plot
sig.labs <- paste("Signature",1:noSig)
names(sig.labs) <- c(1:noSig)

ggplot(dat2, aes(x=type, y=s, fill=type, col = type)) +
  #geom_jitter()+
  geom_boxplot()+
  facet_wrap(time~., nrow = 2, scales = "free_y", labeller = labeller(time = sig.labs))+
  scale_color_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                     name = "Interaction model")+
  scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")+
  xlab("Interaction model")+
  ylab("Cosine similarity")+
  theme_bw()+
  theme(text = element_text(size = 15),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"))


#### UCUT

load("UCUT/result/SimilaritySignaturesUCUTparamboot.RData")
load("UCUT/result/SimilaritySignaturesUCUTmixmodel.RData")

order = c(1,2)
noSig = nrow(ResCosineMatMono)
# Transforming data
datPenta = data.frame(s = t(ResCosineMatPenta[order,]), type = "Penta")
#m = cosMatch(TriRes$Signatures, MixRes$Signatures)$match
datMix = data.frame(s = t(ResCosineMatMix[order,]), type = "Mix")
#m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
datDi = data.frame(s = t(ResCosineMatDi[order,]), type = "Di")
#m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
datMono = data.frame(s = t(ResCosineMatMono[order,]), type = "Mono")
dat1 = rbind(datPenta,datDi,datMix, datMono)
dat2 = reshape(dat1, varying = colnames(dat1)[1:noSig], direction = "long")
dat2$type = factor(dat2$type, levels = c("Mono","Di","Mix","Penta"))

##------------------------------------------------------------------
## Compare true signatures and estimated signatures
##------------------------------------------------------------------

colors3 = c(brewer.pal(n = 6, name = "Dark2")[c(1,3,6)],"#436EEE")[c(1,2,4,3)]

# plot in one plot
sig.labs <- paste("Signature",1:noSig)
names(sig.labs) <- c(1:noSig)

ggplot(dat2, aes(x=type, y=s, fill=type, col = type)) +
  #geom_jitter()+
  geom_boxplot()+
  facet_wrap(time~., nrow = 2, scales = "free_y", labeller = labeller(time = sig.labs))+
  scale_color_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                     name = "Interaction model")+
  scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")+
  xlab("Interaction model")+
  ylab("Cosine similarity")+
  theme_bw()+
  theme(text = element_text(size = 15),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"))
