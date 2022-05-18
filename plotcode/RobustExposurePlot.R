rm(list=ls())
setwd("~/projects/paramNMF_ms")
#################################################
## Plot robustness of exposures
#################################################
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
##################
## cosine similarity 
order21 = c(1,4,3,2)
order214 = c(4,2,3,1)
order = order214

datalist = list()
for(dwn in c(1,2,5)){
  load(paste0("BRCA/result/ExposureBRCA214dwn",dwn,".RData"))
  
  datTri = data.frame(s = t(ResCosineMatTri), type = "Tri")
  #m = cosMatch(TriRes$Signatures, MixRes$Signatures)$match
  datMix = data.frame(s = t(ResCosineMatMix), type = "Mix")
  #m = cosMatch(TriRes$Signatures, DiRes$Signatures)$match
  datDi = data.frame(s = t(ResCosineMatDi), type = "Di")
  #m = cosMatch(TriRes$Signatures, MonoRes$Signatures)$match
  datMono = data.frame(s = t(ResCosineMatMono), type = "Mono")
  
  dat1 = rbind(datTri,datMix,datDi,datMono)
  dat2 = reshape(dat1, varying = colnames(dat1)[1:nrow(ResCosineMatMono)], direction = "long")
  dat2$type = factor(dat2$type, levels = c("Mono","Di","Mix","Tri"))
  dat2$time = factor(dat2$time)
  dat2$downsample = dwn
  datalist[[dwn]] = dat2
}


datall = rbind(datalist[[1]],datalist[[2]],datalist[[5]])
datall$downsample = factor(datall$downsample, levels = c(1,2,5), labels = c("1%","2%","5%"))
##------------------------------------------------------------------
## Compare true exposures and estimated exposures
##------------------------------------------------------------------
#datall = datall[datall$type != "Mix",]
colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")

ggplot(datall, aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1), scale = "width") +
  facet_grid(rows = vars(time), cols = vars(downsample), scales = "free") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +
  ylab("Cosine similarity")+xlab("Interaction model")+
  scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")

ggplot(datall, aes(x=type, y=s, fill=type, col = type)) +
  #geom_violin(position=position_dodge(1), scale = "width") +
  geom_jitter()+
  facet_grid(cols = vars(downsample), scales = "free") +
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank())+
  ggtitle("Downsamling the mutation counts")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +
  ylab("Cosine similarity")+xlab("Interaction model")+
  scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")+
  scale_color_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")

## UCUT exposures --------------------------------------------
datalist = list()
for(dwn in c(1,2,5)){
  load(paste0("UCUT/ExposureUCUTdwn",dwn,".RData"))
  
  datPenta = data.frame(s = t(ResCosineMatPenta), type = "Penta")
  datDi = data.frame(s = t(ResCosineMatDi), type = "Di")
  datMono = data.frame(s = t(ResCosineMatMono), type = "Mono")
  
  dat1 = rbind(datPenta,datDi,datMono)
  dat2 = reshape(dat1, varying = colnames(dat1)[1:nrow(ResCosineMatMono)], direction = "long")
  dat2$type = factor(dat2$type, levels = c("Mono","Di","Penta"))
  dat2$time = factor(dat2$time)
  dat2$downsample = dwn
  datalist[[dwn]] = dat2
}

datall = rbind(datalist[[1]],datalist[[2]],datalist[[5]])
datall$downsample = factor(datall$downsample, levels = c(1,2,5), labels = c("1%","2%","5%"))
##------------------------------------------------------------------
## Compare true exposures and estimated exposures
##------------------------------------------------------------------

colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")

ggplot(datall, aes(x=type, y=s, fill=type)) +
  #geom_violin(position=position_dodge(1), scale = "width") +
  geom_boxplot()+
  facet_grid(rows = vars(time), cols = vars(downsample), scales = "free") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +
  ylab("Signatures")+
  scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")




ggplot(datall, aes(x=type, y=s, fill=type, col = type)) +
  #geom_violin(position=position_dodge(1), scale = "width") +
  geom_jitter()+
  facet_grid(cols = vars(downsample), scales = "free") +
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank())+
  ggtitle("Downsamling the mutation counts")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +
  ylab("Cosine similarity")+xlab("Interaction model")+
  scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                    name = "Interaction model")+
  scale_color_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                     name = "Interaction model")


par(mfrow = c(2,2))
p1 = ggplot(dat2[dat2$time == 1,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 1")
p2 = ggplot(dat2[dat2$time == 2,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 2")
p3 = ggplot(dat2[dat2$time == 3,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 3")
p4 = ggplot(dat2[dat2$time == 4,], aes(x=type, y=s, fill=type)) +
  geom_violin(position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.9) +
  ggtitle("Signature 4")
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2, common.legend = TRUE)





## Mean
mnTri <- rowMeans(ResCosineMatTri)
mnMix <- rowMeans(ResCosineMatMix)
mnDi <- rowMeans(ResCosineMatDi)
mnMono <- rowMeans(ResCosineMatMono)
# plot mean
plot(mnTri,ylim=c(0.5,1),col="red",pch=19,cex=1)
points(mnMix,col="blue",pch=19,cex=0.8)
points(mnDi,col="green",pch=19,cex=0.6)
points(mnMono,col="orange",pch=19,cex=0.5)
## Var
vrTri <- apply(ResCosineMatTri,1,var)
vrMix <- apply(ResCosineMatMix,1,var)
vrDi <- apply(ResCosineMatDi,1,var)
vrMono <- apply(ResCosineMatMono,1,var)
plot(vrTri,col="red",pch=19,cex=0.8)
points(vrMix,col="blue",pch=19,cex=0.6)
points(vrMix,col="blue",pch=19,cex=0.6)
## Quantiles
qTri <- apply(ResCosineMatTri,1,quantile,probs=c(0.05,0.95))
qMix <- apply(ResCosineMatMix,1,quantile,probs=c(0.05,0.95))
## Mean with quantiles

plot(1:21,mnTri,ylim=c(0,1),col="red",pch=19,cex=0.9,
     main=paste("Ratio of down sampling:",dwn),
     xlab="Patient",ylab="Cosine Similarity",cex.lab=1.2)
axis(side=1, at=0:21, labels = TRUE)
points(1:21+0.3,mnMix,col="blue",pch=19,cex=0.9)
for (i in 1:21){
  points(rep(i,2),qTri[,i],col="red",type="l",lwd=2)
  points(rep(i+0.3,2),qMix[,i],col="blue",type="l",lwd=2)
}
lines(c(14,15),rep(0.5,2),col="red",lwd=2)
text(15,0.5,"Sole tri-nucleotide",pos=4,cex=1.2)
lines(c(14,15),rep(0.4,2),col="blue",lwd=2)
text(15,0.4,"Mixture",pos=4,lwd=2,cex=1.2)
points(14.5,0.3,pch=19,col="black",cex=1.2)
text(15,0.3,"Mean",pos=4,cex=1.2)
lines(c(14,15),rep(0.2,2),lty=1,lwd=2)
text(15,0.2,"Quantile interval",pos=4,cex=1.2)
text(1:21,rep(0.05,15),nSimMt,cex=0.8)
text(0.5,0.12,"Number of sampled mutations",pos=4,cex=1.2)

ResultFixDi <- NMFglmSQR(Data = V, NoSignatures = 4, Signatures = DiRes$Signatures,tol=0.2,Seeds=3)
cosMatch(DiRes$Signatures,ResultFixDi$Signatures)