rm(list=ls())
setwd("~/projects/paramNMF_ms")
###################################
## Plotting signatures
###################################
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source("~/projects/paramNMF_ms/ModelSelection.R")

################################################################
## BRCA214 signature plot
################################################################

source("BRCA/loadBRCA214models.R")
load("C:/Users/au543194/Documents/projects/paramNMF_ms/BRCA/result/notUseResults/BRCA214sig8/BRCA214optimal8sig.RData")

resMono = resFactors[[1]]
resDi = resFactors[[2]]
resTri = resFactors[[3]]
resMix = resFactors[[4]]

order = c(7,8,3,2,5,1,4,6)
resTri$signatures = resTri$signatures[order,]
mMono = cosMatch(resTri$signatures, resMono$signatures)$match
mMix = cosMatch(resTri$signatures, resMix$signatures)$match
mDi = cosMatch(resTri$signatures, resDi$signatures)$match
mTri = cosMatch(resTri$signatures, resTri$signatures)$match

noSig = 8
#sum(resMix$exposures == 0)
###################################################################################################
## plot signatures

sub = factor(rep(c("C > A","C > G","C > T","T > A","T > C","T > G"), each = 16))
# dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S = cbind(res1$Signatures[,1],
#                                                                                       res2$Signatures[,3],
#                                                                                       res3$Signatures[,3]))

#col.model = c(brewer.pal(7, name = "Blues")[c(4)], brewer.pal(7, name = "Oranges")[5], brewer.pal(7, name = "Blues")[c(7)])
# mono
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMono$signatures[mMono,]))
dat$type = "mono"
datmono = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")
# mix
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMix$signatures[mMix,]))
dat$type = "mix"
datmix = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")
# di
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resDi$signatures[mDi,]))
dat$type = "di"
datdi = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")
# tri
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resTri$signatures[mTri,]))
dat$type = "tri"
dattri = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")

datall = rbind(datmono,datdi,datmix,dattri)


col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")

time.labs <- c("Signature 1   (SBS2)", "Signature 2   (SBS5)", "Signature 3   (SBS8)", "Signature 4   (SBS13)", "Signature 5   (SBS17b)", "Signature 6   (SBS18)", "Signature 7   (SBS1)", "Signature 8   (SBS39)")
names(time.labs) <- c("1", "2", "3", "4", "5", "6", "7", "8")

cossim = matrix(c(0.80,0.90,0.94,0.86,0.83,0.94,0.99,0.90,
                  0.92,0.91,0.94,0.91,0.87,0.88,1.00,0.82,
                  0.92,0.98,0.90,0.87,0.87,0.91,1.00,0.86,
                  0.14,0.98,0.77,0.94,0.88,0.79,0.99,0.66),byrow = T, nrow = 4)


cossim = cossim[,order]

datall$comb = paste0(datall$time,datall$sub)
datall$type = factor(datall$type, levels = c("mono","di","mix","tri"))
colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")
plots = list()
for(sig in 1:8){
  
  annotations <- data.frame(
    xpos = rep(Inf,384),
    ypos =  rep(Inf,384),
    annotateText = c(c(rep("",95),paste0("(",cossim[4,sig],")")),c(rep("",95),paste0("(",cossim[2,sig],")")),c(rep("",95),paste0("(",cossim[3,sig],")")),c(rep("",95),paste0("(",cossim[1,sig],")"))),
    hjustvar = rep(1,384),
    vjustvar = rep(1,384), type = c(rep("mono",96),rep("di",96),rep("mix",96),rep("tri",96)),
    sub = rep(c("C > A","C > G","C > T","T > A","T > C","T > G"), each = 16))
  
  annotations$type = factor(annotations$type, levels = c("mono","di","mix","tri"))
  
  plots[[sig]] = ggplot(datall[datall$Signature == sig, ], aes(x = m, y = true, fill = type))+
    geom_bar(stat = "identity", width = 0.5)+
    facet_grid(rows = vars(type), cols = vars(sub), scale = "free_x", switch = "x",labeller = labeller(time = time.labs))+
    theme_bw()+
    theme(text = element_text(size=9, face = "bold"), 
          axis.text.x=element_blank(),
          strip.text = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          legend.text = element_text(size = 10, face = "bold"),
          #legend.position = "none", 
          #strip.text.y.right = element_text(angle = -90, size = 8),
          strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
          strip.text.x = element_text(size = 9), panel.spacing.x = unit(0.2,"line"))+ 
    ylab("")+xlab("")+
    ggtitle(time.labs[sig])+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)+
    scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                      name = "Interaction model")+
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), cex = 3, fontface = "plain")
}

p3 = ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]], nrow = 4, ncol = 2, common.legend = T, legend = "bottom")



p3 +
  theme(plot.margin = margin(15,0.1,2,0.1)) +
  geom_text(aes(x = 0.19, y = 0.98, hjust = 0, vjust = 0, label = "Corresponding COSMIC signature"), cex = 3.5)+
  geom_text(aes(x = 0.17, y = 0.98, hjust = 0, vjust = 0, label = sprintf('\u2190')), cex = 4.5) +
  geom_text(aes(x = 0.3, y = 0.955, hjust = 0, vjust = 0, label = "Cosine similarity to COSMIC"), cex = 3.3)+
  geom_text(aes(x = 0.45, y = 0.955, hjust = 0, vjust = 0, label = sprintf('\u2192')), cex = 4.5)


######################################################
## UCUT 
#####################################################
  
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
source("UCUT/loadUCUTmodels.R")

resMono = resFactors[[1]]
resDi = resFactors[[19]]
resPenta = resFactors[[21]]

load("UCUT/result/UCUTditriMIXmodel.RData")

resMix = resFactors[[5]]

mMono = cosMatch(resPenta$signatures, resMono$signatures)$match
mDi = cosMatch(resPenta$signatures, resDi$signatures)$match
mMix = cosMatch(resPenta$signatures, resMix$signatures)$match
noSig = 2

###################################################################################################
## plot signatures

sub = factor(rep(c("C > A","C > G","C > T","T > A","T > C","T > G"), each = 256))

# mono
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMono$signatures[mMono,]))
dat$type = "mono"
datmono = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")
# di
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resDi$signatures[mDi,]))
dat$type = "di"
datdi = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")
# penta
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resPenta$signatures))
dat$type = "penta"
datpenta = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")
# mix 
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMix$signatures[mMix,]))
dat$type = "mix"
datmix = reshape(dat, varying = paste0("S.",c(1:noSig)), direction = "long", v.names = "true", timevar = "Signature")


datall = rbind(datmono,datdi,datmix,datpenta)


col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")


datall$comb = paste0(datall$time,datall$sub)
datall$type = factor(datall$type, levels = c("mono","di","mix","penta"))
colors3 = c(brewer.pal(n = 6, name = "Dark2")[c(1,3,6)],"#436EEE")[c(1,2,4,3)]

plots = list()
for(sig in 1:noSig){
  plots[[sig]] = ggplot(datall[datall$Signature == sig, ], aes(x = m, y = true, fill = type))+
    geom_bar(stat = "identity", width = 0.5)+
    facet_grid(rows = vars(type), cols = vars(sub), scale = "free_x", switch = "x")+
    #theme_bw()+
    theme(text = element_text(size=8, face = "bold"), 
          axis.text.x=element_blank(),
          strip.text = element_blank(),
          axis.ticks.x = element_blank(), 
          #legend.position = "none", 
          #strip.text.y.right = element_text(angle = -90, size = 8),
          strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
          strip.text.x = element_text(size = 9), panel.spacing.x = unit(0.2,"line"))+ 
    ylab("")+xlab("")+
    #ylim(0,0.03)+
    ggtitle(paste("Signature",sig))+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)+
    scale_fill_manual(values = colors3, labels = c("Mono","Di","Mix","Penta"), 
                      name = "Interaction model")
}
ggarrange(plots[[1]],plots[[2]], ncol = 1, common.legend = T, legend = "bottom")



################################################################
## BRCA21 signature plot
################################################################ 
source("BRCA/loadBRCA21models.R")
resMono = resFactors[[1]]
resMix = resFactors[[8]]
resDi = resFactors[[11]]
resTri = resFactors[[15]]
resTri$signatures = resTri$signatures[c(1,4,3,2),]
mMono = cosMatch(resTri$signatures, resMono$signatures)$match
mMix = cosMatch(resTri$signatures, resMix$signatures)$match
mDi = cosMatch(resTri$signatures, resDi$signatures)$match
mTri = cosMatch(resTri$signatures, resTri$signatures)$match


# plot signatures

sub = factor(rep(c("C > A","C > G","C > T","T > A","T > C","T > G"), each = 16))
# dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S = cbind(res1$Signatures[,1],
#                                                                                       res2$Signatures[,3],
#                                                                                       res3$Signatures[,3]))

col.model = c(brewer.pal(7, name = "Blues")[c(4)], brewer.pal(7, name = "Oranges")[5], brewer.pal(7, name = "Blues")[c(7)])
# mono
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMono$signatures[mMono,]))
dat$type = "mono"
datmono = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")
# mix
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMix$signatures[mMix,]))
dat$type = "mix"
datmix = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")
# di
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resDi$signatures[mDi,]))
dat$type = "di"
datdi = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")
# tri
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resTri$signatures[mTri,]))
dat$type = "tri"
dattri = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")

datall = rbind(datmono,datdi,datmix,dattri)


time.labs <- c("Signature 1", "Signature 2", "Signature 3", "Signature 4")
names(time.labs) <- c("1", "2", "3", "4")

col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")
#col.model <- brewer.pal(5, name = "Set1")[c(1,3,3)]

#col.comb = c(brewer.pal(7, name = "Reds")[-1],brewer.pal(7, name = "Blues")[-1],brewer.pal(7, name = "Greens")[-1])
#col.comb = c(rep(brewer.pal(7, name = "Greens")[-1],3))
col.comb = c(rep(col.model[3],6),rep(col.model[3],6),rep(col.model[3],6),rep(col.model[3],6))

datall$comb = paste0(datall$time,datall$sub)
datall$type = factor(datall$type, levels = c("mono","di","mix","tri"))

plots = list()
colors2 = c("#438EEE", "#6A0136","#FF9912", "#698F3F")
colors = c("#800080", "#436EEE","#FF9912", "#27408B")
colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")
col.list = list(colors[c(1,1,3,4)],colors[c(1,4,3,4)],colors[c(1,3,3,4)],colors[c(1,3,3,4)])
for(sig in 1:4){
  plots[[sig]] = ggplot(datall[datall$Signature == sig, ], aes(x = m, y = true, fill = type))+
    geom_bar(stat = "identity", width = 0.5)+
    facet_grid(rows = vars(type), cols = vars(sub), scale = "free_x", switch = "x",labeller = labeller(time = time.labs))+
    theme_bw()+
    theme(text = element_text(size=8, face = "bold"), 
          axis.text.x=element_blank(),
          strip.text = element_blank(),
          axis.ticks.x = element_blank(),
          #legend.position = "none", 
          #strip.text.y.right = element_text(angle = -90, size = 8),
          strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
          strip.text.x = element_text(size = 8), 
          panel.spacing.x = unit(0.2,"line"),
          panel.border = element_blank())+ 
    ylab("")+xlab("")+
    ggtitle(paste("Signature",sig))+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)+
    scale_fill_manual(values = colors3, labels = c("Mono","Di",'Mix',"Tri"), 
                      name = "Interaction model")
}
ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]], nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
