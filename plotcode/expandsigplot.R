rm(list=ls())
setwd("~/projects/paramNMF_ms")
###################################
## Plotting signatures
###################################
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source("GitModelSelection.R")

# BRCA21
source("BRCA/loadBRCA21models.R")

resTri = resFactors[[15]]
#resTri$signatures = resTri$signatures[c(1,4,3,2),]

###################################################################################################
## plot signatures

sub = factor(rep(c("C > A","C > G","C > T","T > A","T > C","T > G"), each = 16))

col.model = c(brewer.pal(7, name = "Blues")[c(4)], brewer.pal(7, name = "Oranges")[5], brewer.pal(7, name = "Blues")[c(7)])

# tri
dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(resTri$signatures))
dat$type = "tri"
dattri = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")
dattri$Signature = factor(dattri$Signature)
time.labs <- c("Signature 1", "Signature 2", "Signature 3", "Signature 4")
names(time.labs) <- c("1", "2", "3", "4")

col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")
#col.model <- brewer.pal(5, name = "Set1")[c(1,3,3)]

#col.comb = c(brewer.pal(7, name = "Reds")[-1],brewer.pal(7, name = "Blues")[-1],brewer.pal(7, name = "Greens")[-1])
#col.comb = c(rep(brewer.pal(7, name = "Greens")[-1],3))
col.comb = c(rep(col.model[3],6),rep(col.model[3],6),rep(col.model[3],6),rep(col.model[3],6))

#datall$comb = paste0(datall$time,datall$sub)

plots = list()
colors2 = c("#438EEE", "#6A0136","#FF9912", "#698F3F")
colors = c("#800080", "#436EEE","#FF9912", "#27408B")
colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")
col.list = list(colors[c(1,1,3,4)],colors[c(1,4,3,4)],colors[c(1,3,3,4)],colors[c(1,3,3,4)])
for(sig in 1:4){
  plots[[sig]] = ggplot(dattri[dattri$Signature == sig, ], aes(x = m, y = true, fill = sub))+
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
          panel.border = element_blank(), 
          plot.title = element_text(vjust = - 5, hjust = 1), 
          legend.position = element_blank())+ 
    ylab("")+xlab("")+
    ggtitle(paste("Signature",sig))+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)+
    scale_fill_manual(values = col.sub)
}
ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]], ncol = 1, legend = "none",widths = c(1, 0.05, 1))


ggplot(dattri, aes(x = m, y = true, fill = Signature))+
  geom_bar(stat = "identity", width = 0.5)+
  facet_grid(rows = vars(Signature), cols = vars(sub), scale = "free", switch = "x",labeller = labeller(time = time.labs))+
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
        panel.border = element_blank(), 
        plot.title = element_text(vjust = - 5, hjust = 1), 
        legend.position = "none")+ 
  ylab("")+xlab("")+
  #geom_text(mapping = aes(x = 0.9, y = c(0.1,0.2,0.3,0.4), labels = paste("Signature",c(1:4)) ))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)
  #scale_fill_manual(values = col.sub)


##################################################
## plot exposures 

# tri
dat = data.frame(m = factor(rownames(V), levels=unique(rownames(V))), S =  resTri$exposures)
dattri = reshape(dat, varying = c("S.1","S.2","S.3","S.4"), direction = "long", v.names = "true", timevar = "Signature")

time.labs <- c("Signature 1", "Signature 2", "Signature 3", "Signature 4")
names(time.labs) <- c("1", "2", "3", "4")

col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")

#datall$comb = paste0(datall$time,datall$sub)

colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")

dattri$Signature = factor(dattri$Signature)
dattri$m = factor(dattri$m, levels = unique(rownames(V)[21:1]))

ggplot(dattri, aes(x = m, y = true, fill = Signature))+
  geom_bar(stat = "identity", width = 0.8, position = "fill")+
  coord_flip()+
  theme_bw()+
  theme(text = element_text(size=8, face = "bold"), 
        #axis.text.x=element_blank(),
        strip.text = element_blank(),
        #axis.ticks.x = element_blank(),
        #legend.position = "none", 
        #strip.text.y.right = element_text(angle = -90, size = 8),
        strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
        strip.text.x = element_text(size = 8), 
        panel.spacing.x = unit(0.2,"line"),
        panel.border = element_blank(), 
        plot.title = element_text(vjust = - 5, hjust = 1), 
        legend.position = "none")+
  ylab("Number of mutations")+xlab("")


  facet_grid(rows = vars(Signature), cols = vars(sub), scale = "free", switch = "x",labeller = labeller(time = time.labs))+
  + 
  ylab("")+xlab("")+
  #geom_text(mapping = aes(x = 0.9, y = c(0.1,0.2,0.3,0.4), labels = paste("Signature",c(1:4)) ))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)+
  scale_fill_manual(values = col.sub)


###########################################################################
## Plot data 

  sub = factor(rep(c("C > A","C > G","C > T","T > A","T > C","T > G"), each = 16))
  
  col.model = c(brewer.pal(7, name = "Blues")[c(4)], brewer.pal(7, name = "Oranges")[5], brewer.pal(7, name = "Blues")[c(7)])
  
  # tri
  dat = data.frame(m = factor(colnames(V), levels=unique(colnames(V))), sub, S =  t(V))
  dattri = reshape(dat, varying = colnames(dat)[-c(1,2)], direction = "long", v.names = "true", timevar = "Signature")
  dattri$Signature = factor(dattri$Signature)
  #time.labs <- c("Signature 1", "Signature 2", "Signature 3", "Signature 4")
  #names(time.labs) <- c("1", "2", "3", "4")
  
  col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")
 
  col.comb = c(rep(col.model[3],6),rep(col.model[3],6),rep(col.model[3],6),rep(col.model[3],6))
  
  #datall$comb = paste0(datall$time,datall$sub)
  
  plots = list()
  colors2 = c("#438EEE", "#6A0136","#FF9912", "#698F3F")
  colors = c("#800080", "#436EEE","#FF9912", "#27408B")
  colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")
  
  
  ggplot(dattri, aes(x = m, y = true))+
    geom_bar(stat = "identity", width = 0.5)+
    facet_grid(rows = vars(Signature), cols = vars(sub), scale = "free", switch = "x")+
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
          panel.border = element_blank(), 
          plot.title = element_text(vjust = - 5, hjust = 1), 
          legend.position = "none")+ 
    ylab("")+xlab("")
    #geom_text(mapping = aes(x = 0.9, y = c(0.1,0.2,0.3,0.4), labels = paste("Signature",c(1:4)) ))+
    #scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)
  #scale_fill_manual(values = col.sub) 

#########################################################################
#### Signature plot of SFS area from sampling algorithm
##########################################################################
library(SFS)
dev1 = sampleSFS(P=t(resTri$signatures), E = t(resTri$exposures), maxIter = 10^5, beta = 0.5, check = 1000)
prob.min = dev1$Pminimum
prob.max = dev1$Pmaximum
N = 4
V = t(V)
# Converting into right data frame
sig = rep(c(1:N), dev1$totalIter)
mut <- c("C > A","C > G","C > T","T > A","T > C","T > G")
sub = rep(mut, each = 16)
dat1 = data.frame(m = factor(rownames(V), levels=unique(rownames(V))), sub = sub, S = t(prob.min))
datmin = reshape(dat1, varying = colnames(dat1)[-c(1,2)], direction = "long", v.names = "min")
dat1 = data.frame(m = factor(rownames(V), levels=unique(rownames(V))),S = t(prob.max))
datmax = reshape(dat1, varying = colnames(dat1)[-1], direction = "long", v.names = "max")
data2 = merge(datmin,datmax, by = c("m","time","id"))
data2$time = factor(data2$time)

equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}
## Plotting the signatures
col.sub = c("#800080", "#FF9912", "#436EEE", "#ffdf12", "#27408B", "#E066FF")
g1 = ggplot(data2, aes(x = m, y = min, fill = sub))+
  #geom_errorbar(aes(ymin = min, ymax = min, col = time), width = 1.05, size = 1)+
  #geom_errorbar(aes(ymin = max, ymax = max, col = time), width = 1.05, size = 1)+
  geom_bar(aes(x = m, y = max), stat = "identity", width = 0.78, fill = "black")+
  #geom_errorbar(aes(ymin = min, ymax = max, col = time), lwd = 1, size = 0.5, linetype = "11")+
  geom_bar(stat = "identity", width = 0.8)+
  facet_grid(rows = vars(time), cols = vars(sub), scales = "free", switch = "x")+theme_bw()+
  theme(text = element_text(size=12, face = "bold"), axis.text.x=element_blank(),axis.text.y = element_text(size = 8),axis.ticks = element_blank(), 
        legend.position = "none",strip.text.y.right = element_text(angle = 0, size = 15), panel.spacing.x = unit(0,"line"),
        strip.background.x = element_rect(color="black", fill="white",linetype="blank"),
        strip.text.x = element_text(size = 9))+ 
  ylab("Probability")+xlab("Mutation types")+ggtitle("Signatures")+
  scale_fill_manual(values = col.sub)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=equal_breaks(n=3, s=0.2), 
                     expand = c(0.05, 0))
#scale_colour_grey()+
#scale_fill_grey()
g1

# include colors on strips for N = 3
g11 <- ggplot_gtable(ggplot_build(g1))
stripr <- which(grepl('strip-r', g11$layout$name))
fills <- c("#E69F00", "#56B4E9", "#009E73")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g11$grobs[[i]]$grobs[[1]]$childrenOrder))
  g11$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


####################################################################
## exposure plot
##################################################################

## equivalent to signature plot
prob.min = dev1$Eminimum%*%diag(1/colSums(dev1$E_lastCheckResults[1:N,]))
prob.max = dev1$Emaximum%*%diag(1/colSums(dev1$E_lastCheckResults[1:N,]))
dat1 = data.frame(m = factor(colSums(M), levels=unique(colSums(M))), S = t(prob.min))
datmin = reshape(dat1, varying = colnames(dat1)[-1], direction = "long", v.names = "min")
dat1 = data.frame(m = factor(colSums(M), levels=unique(colSums(M))),S = t(prob.max))
datmax = reshape(dat1, varying = colnames(dat1)[-1], direction = "long", v.names = "max")

dat2 = merge(datmin,datmax, by = c("m","time","id"))
dat2$time = factor(dat2$time)

## different plots 
g2 = ggplot(dat2, aes(x = m, y = min))+
  #geom_errorbar(aes(ymin = min, ymax = min, col = time), width = 1, size = 1)+
  #geom_errorbar(aes(ymin = max, ymax = max, col = time), width = 1, size = 1)+
  geom_bar(aes(x = m, y = max), stat = "identity", width = 0.78, fill = "tomato2")+
  #geom_errorbar(aes(ymin = min, ymax = max, col = time), lwd = 1, size = 0.5, linetype = "11")+
  geom_bar(stat = "identity", width = 0.8, fill = "grey")+
  facet_grid(cols = vars(time))+
  theme_bw()+
  theme(text = element_text(size=12, face = "bold"), axis.text.x=element_text(angle = 315, size = 5, hjust = 0.7, vjust = -0.6), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "none",strip.text.y.right = element_text(angle = 0))+ 
  ylab("Probability")+xlab("Patients")+ggtitle("   Normalized Exposures")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0,0.5,1))+
  coord_flip()
#scale_colour_grey()+
#scale_fill_grey()

g2

# include colors on strips for N = 3
g22 <- ggplot_gtable(ggplot_build(g2))
stripr <- which(grepl('strip-t', g22$layout$name))
fills <- c("#E69F00", "#56B4E9", "#009E73")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g22$grobs[[i]]$grobs[[1]]$childrenOrder))
  g22$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# final plot
ggarrange(g11,g22, labels = c("(A)","(B)"), widths = c(2, 1))