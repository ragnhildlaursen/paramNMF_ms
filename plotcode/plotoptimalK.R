setwd("~/projects/paramNMF_ms/")
rm(list=ls())

####################################
## Find optimal models BRCA214
####################################
library(Rcpp)
library(RcppArmadillo)
library(bench)
library(gtools)
library(ggplot2)
library(reshape2)
library(dplyr)


source("BRCA/loadBRCA214models.R")


sig = c(1:15)


load("BRCA/result/BRCA214allmodels20init.RData")


# find model with minimum BIC
optimallist = list()

bicmin = c()
bicmono = c()
bicfull = c()
bicdi = c()

for(i in 1:15){
  optimallist[[i]] = allmodels[[i]][,which.min(allmodels[[i]][4+i,])]
  bicmin[i] = min(allmodels[[i]][4+i,])
  bicmono[i] = allmodels[[i]][4+i,1]
  bicfull[i] = allmodels[[i]][4+i,which(allmodels[[i]][2+i,] == 96*i)]
  if(i > 1){
    bicdi[i] = allmodels[[i]][4+i,which(colSums(allmodels[[i]][c(2:(1+i)),] == rep(42,i)) == i)]
  }else{
    bicdi[i] = allmodels[[i]][4+i,2]
  }
}

dat1 = data.frame(sig, flex = bicmin, mono = bicmono, di = bicdi, full = bicfull)
dat2 = reshape(dat1, varying = list(names(dat1)[-1]), times = names(dat1)[-1], direction = "long",v.names = "BIC", timevar = "Model")

dat2 %>%
  group_by(Model) %>%
  mutate(min = min(BIC),
         minsig = sig[which.min(BIC)]) -> dat4


dat5 = dat4

dat6 = dat5[dat5$sig != 1,]


# include the optimal flexible model
model = matrix(0, nrow = 15, ncol = 3)


for(i in 1:15){
  model[i,] = table(factor(optimallist[[i]][-c(1,(i+2):(i+4))], levels = c(12,42,96)))
}

model2 = model[-1,]

dat6$mono = rep(model2[,1],4)
dat6$di = rep(model2[,2],4)
dat6$full = rep(model2[,3],4)

dat7 = dat6[dat6$sig > 4,]
dat7$Model = factor(dat7$Model, levels = c("mono","di","full","flex"),labels = c("Mono-nucleotide","Di-nucleotide","Tri-nucleotide","Mixture"))

p2 = ggplot(dat7, aes(x = sig, y = BIC, group = Model, color = Model))+
  geom_line(lwd = 0.7)+
  geom_label(aes(x = sig, y = exp(11.7), label = full), fill = "#27408B", col = "white", fontface = "bold",cex = 3)+
  geom_label(aes(x = sig, y = exp(11.77), label = di), fill = "#FF9912", col = "white", fontface = "bold",cex = 3)+
  geom_label(aes(x = sig, y = exp(11.84), label = mono), fill = "#AC0136", col = "white", fontface = "bold",cex = 3)+
  ggtitle("BRCA")+
  xlab("Number of signatures")+
  ylab("Bayesian Information Criterion (BIC)")+
  geom_point(aes(x = minsig, y = min), cex = 2)+
  scale_color_manual(name = "Model for signatures",values = c("#AC0136", "#FF9912", "#27408B", "gray30"))+
  #theme(legend.title = element_text("Model for signatures"))+
  scale_y_log10()+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_x_discrete(limits = c(5:15))

p2 + annotate("text",x = 10, y = exp(11.93), label = "Optimal mixture", col = "gray30")+ 
  annotate("text",x = 4.3, y = exp(11.7), label = "Tri ", col = "gray30", cex = 3, fontface = "bold")+
  annotate("text",x = 4.3, y = exp(11.77), label = "Di  ", col = "gray30", cex = 3, fontface = "bold")+
  annotate("text",x = 4.3, y = exp(11.84), label = "Mono", col = "gray30", cex = 3, fontface = "bold")


##################################
## Find optimal model UCUT
##################################
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(gtools)

source("UCUT/loadUCUTmodels.R")


sig = c(1:7)

library(ggplot2)
library(reshape2)
library(dplyr)


load("UCUT/result/UCUTallmodels3init.RData")


# find model with minimum BIC
optimallist = list()

bicmin = c()
bicmono = c()
bicdi = c()

for(i in sig){
  optimallist[[i]] = allmodels[[i]][,which.min(allmodels[[i]][4+i,])]
  bicmin[i] = min(allmodels[[i]][4+i,])
  bicmono[i] = allmodels[[i]][4+i,1]
  if(i > 1){
    bicdi[i] = allmodels[[i]][4+i,which(colSums(allmodels[[i]][c(2:(1+i)),] == rep(66,i)) == i)]
  }else{
    bicdi[i] = allmodels[[i]][4+i,which(allmodels[[i]][c(2:(1+i)),] == 66)]
  }
}
#bicfull = as.vector(dat3[which(dat3[,2] == "full")[1:7],3])$BIC

dat1 = data.frame(sig, flex = bicmin, mono = bicmono, di = bicdi)
dat2 = reshape(dat1, varying = list(names(dat1)[-1]), times = names(dat1)[-1], direction = "long",v.names = "BIC", timevar = "Model")

dat2 %>%
  group_by(Model) %>%
  mutate(min = min(BIC),
         minsig = sig[which.min(BIC)]) -> dat4


dat5 = dat4
dat6 = dat5

# include the optimal flexible model
model = matrix(0, nrow = 7, ncol = 6)


for(i in 1:7){
  model[i,] = table(factor(optimallist[[i]][-c(1,(i+2):(i+4))], levels = c(18,48,66,96,102,120)))
}

model2 = model
#model2[model2 == 0] = NA
dat6$mono = rep(model2[,1],3)
dat6$dimono = rep(model2[,2],3)
dat6$di = rep(model2[,3],3)
dat6$ditri = rep(model2[,6],3)

dat7 = dat6[dat6$Model != "full",]
dat7$Model = factor(dat7$Model, levels = c("mono","di","flex"),labels = c("Mono-nucleotide","Di-nucleotide","Mixture"))


library(RColorBrewer)

col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
p3 = ggplot(dat7, aes(x = sig, y = BIC, group = Model, color = Model))+
  geom_line(lwd = 0.7)+
  geom_label(aes(x = sig, y = exp(10.28), label = mono), fill = "#1B9E77", col = "white", fontface = "bold",cex = 3)+
  geom_label(aes(x = sig, y = exp(10.25), label = di), fill = "#7570B3", col = "white", fontface = "bold",cex = 3)+
  geom_label(aes(x = sig, y = exp(10.22), label = dimono), fill = "#D95F02", col = "white", fontface = "bold",cex = 3)+
  geom_label(aes(x = sig, y = exp(10.19), label = ditri), fill = "#E6AB02", col = "white", fontface = "bold",cex = 3)+
  ggtitle("UCUT")+
  xlab("Number of signatures")+
  ylab("Bayesian Information Criterion (BIC)")+
  geom_point(aes(x = minsig, y = min), cex = 2)+
  scale_color_manual(name = "Model for signatures",values = c("#1B9E77", "#7570B3", "gray30"))+
  #theme(legend.title = element_text("Model for signatures"))+
  scale_y_log10()+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_x_discrete(limits = c(1:7))

p3 + annotate("text",x = 4, y = exp(10.31), label = "Optimal mixture", col = "gray30")+ 
  annotate("text",x = 0.5, y = exp(10.28), label = "Mono", col = "gray30", cex = 3, fontface = "bold")+
  annotate("text",x = 0.5, y = exp(10.25), label = "Di", col = "gray30", cex = 3, fontface = "bold")+
  annotate("text",x = 0.5, y = exp(10.22), label = "Di and mono", col = "gray30", cex = 3, fontface = "bold")+
  annotate("text",x = 0.5, y = exp(10.19), label = "Tri and di", col = "gray30", cex = 3, fontface = "bold")

# annotate("text",x = 0.1, y = exp(10.28), label = "Mono", col = "gray30", cex = 3, fontface = "bold")+
#   annotate("text",x = -0.03, y = exp(10.25), label = "Di", col = "gray30", cex = 3, fontface = "bold")+
#   annotate("text",x = 0.31, y = exp(10.22), label = "Di and mono", col = "gray30", cex = 3, fontface = "bold")+
#   annotate("text",x = 0.215, y = exp(10.19), label = "Tri and di", col = "gray30", cex = 3, fontface = "bold")
