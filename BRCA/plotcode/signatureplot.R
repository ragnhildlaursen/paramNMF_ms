rm(list=ls())
###################################
## Plotting signatures
###################################
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source("~/projects/paramNMF_ms/GitModelSelection.R")

# BRCA21
source("BRCA/loadBRCA21models.R")
resMono = resFactors[[1]]
resMix = resFactors[[9]]
resDi = resFactors[[11]]
resTri = resFactors[[15]]

mMono = cosMatch(resTri$Signatures, resMono$Signatures)$match
mMix = cosMatch(resTri$Signatures, resMix$Signatures)$match
mDi = cosMatch(resTri$Signatures, resDi$Signatures)$match


###################################################################################################
## plot signatures

sub = factor(rep(c("C > A","C > G","C > T","T > A","T > C","T > G"), each = 16))
# dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S = cbind(res1$Signatures[,1],
#                                                                                       res2$Signatures[,3],
#                                                                                       res3$Signatures[,3]))

col.model = c(brewer.pal(7, name = "Blues")[c(4)], brewer.pal(7, name = "Oranges")[5], brewer.pal(7, name = "Blues")[c(7)])
# mono
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMono$Signatures[mMono,]))
dat$type = "mono"
datmono = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")
# mix
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMix$Signatures[mMix,]))
dat$type = "mix"
datmix = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")
# di
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resDi$Signatures[mDi,]))
dat$type = "di"
datdi = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")
# tri
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resTri$Signatures))
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
datall$type = factor(datall$type, levels = c("mono","mix","di","tri"))

plots = list()
colors2 = c("#438EEE", "#6A0136","#FF9912", "#698F3F")
colors = c("#800080", "#436EEE","#FF9912", "#27408B")
colors3 = c("#AC0136","#436EEE", "#FF9912", "#27408B")
col.list = list(colors[c(1,1,3,4)],colors[c(1,4,3,4)],colors[c(1,3,3,4)],colors[c(1,3,3,4)])
for(sig in 1:4){
plots[[sig]] = ggplot(datall[datall$Signature == sig, ], aes(x = m, y = true, fill = type))+
  geom_bar(stat = "identity", width = 0.5)+
  facet_grid(rows = vars(type), cols = vars(sub), scale = "free", switch = "x",labeller = labeller(time = time.labs))+
  theme(text = element_text(size=10, face = "bold"), 
        axis.text.x=element_blank(),
        strip.text = element_blank(),
        axis.ticks = element_blank(), 
        #legend.position = "none", 
        #strip.text.y.right = element_text(angle = -90, size = 8),
        strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
        strip.text.x = element_text(size = 9), panel.spacing.x = unit(0,"line"))+ 
  ylab("")+xlab("")+
  ggtitle(paste("Signature",sig))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_fill_manual(values = colors3)
}
ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]], nrow = 2, ncol = 2, common.legend = T, legend = "bottom")

################################################################
## BRCA214 signature plot
################################################################

source("BRCA/loadBRCA214models.R")

resMono = resFactors[[1]]
resMix = resFactors[[13]]
resDi = resFactors[[36]]
resTri = resFactors[[45]]

mMono = cosMatch(resTri$signatures, resMono$signatures)$match
mMix = cosMatch(resTri$signatures, resMix$signatures)$match
mDi = cosMatch(resTri$signatures, resDi$signatures)$match


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
datmono = reshape(dat, varying = paste0("S.",c(1:8)), direction = "long", v.names = "true", timevar = "Signature")
# mix
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resMix$signatures[mMix,]))
dat$type = "mix"
datmix = reshape(dat, varying = paste0("S.",c(1:8)), direction = "long", v.names = "true", timevar = "Signature")
# di
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resDi$signatures[mDi,]))
dat$type = "di"
datdi = reshape(dat, varying = paste0("S.",c(1:8)), direction = "long", v.names = "true", timevar = "Signature")
# tri
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resTri$signatures))
dat$type = "tri"
dattri = reshape(dat, varying = paste0("S.",c(1:8)), direction = "long", v.names = "true", timevar = "Signature")

datall = rbind(datmono,datdi,datmix,dattri)


col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")


datall$comb = paste0(datall$time,datall$sub)
datall$type = factor(datall$type, levels = c("mono","mix","di","tri"))

plots = list()
for(sig in 1:8){
  plots[[sig]] = ggplot(datall[datall$Signature == sig, ], aes(x = m, y = true, fill = type))+
    geom_bar(stat = "identity", width = 0.5)+
    facet_grid(rows = vars(type), cols = vars(sub), scale = "free", switch = "x")+
    theme(text = element_text(size=8, face = "bold"), 
          axis.text.x=element_blank(),
          strip.text = element_blank(),
          axis.ticks = element_blank(), 
          #legend.position = "none", 
          #strip.text.y.right = element_text(angle = -90, size = 8),
          strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
          strip.text.x = element_text(size = 7), 
          panel.spacing.x = unit(0,"line"))+ 
    ylab("")+xlab("")+
    ggtitle(paste("Signature",sig))+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    scale_fill_manual(values = colors3)
}
ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],
          plots[[5]],plots[[6]],plots[[7]],plots[[8]],
          nrow = 2, ncol = 4, common.legend = T, legend = "bottom")




