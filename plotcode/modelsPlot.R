rm(list=ls())
#########################################
## Plot of models for BRCA 214 
#########################################
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(grid)
setwd("~/projects/paramNMF_ms/")
#
############## BRCA model plot ------------------------------------------------------


source("BRCA/loadBRCA214models.R")

load("BRCA/result/BRCA214allmodels20init.RData")

resMat = t(allmodels[[8]])

#load("UCUT/UCUTmodelsummary.RData")

noSig = 8
# colors for signatures
col.model = c("#AC0136", "#FF9912", "#27408B")

resMatnew = resMat
resMatnew[resMatnew == 12] = col.model[1]
resMatnew[resMatnew == 42] = col.model[2]
resMatnew[resMatnew == 96] = col.model[3]

# parameter order
prm_order = order(resMat[,10])

gkl = as.numeric(resMat[,11])

noModels = nrow(resMat)


diff = diff(range(gkl))/10
## ggplot
dat1 = data.frame(idx = c(1:length(gkl)), gkl = gkl[prm_order], sigty = resMat[prm_order,c(1:noSig)+1])
dat2 = reshape(dat1, varying = colnames(dat1)[-c(1,2)], direction = "long")
dat2$sigty = factor(dat2$sigty)

dat2$bic = 2*dat2$gkl + log(sum(V > 0.5))*(as.numeric(resMat[prm_order,10]) + nrow(V)*noSig) 


mC = (max(dat2$gkl) - min(dat2$gkl))/(max(dat2$bic) - min(dat2$bic))
dat2$bic = dat2$bic* mC
addC = min(dat2$gkl) - min(dat2$bic)
dat2$bic = dat2$bic + addC
datparam = data.frame(idx = c(1:nrow(resMat)), param = resMat[prm_order,10], y = max(gkl)+diff/1.2)
p = ggplot(dat2, aes(x = idx, y = gkl))+
  #geom_vline(xintercept = 16, lty = "dashed")+
  geom_point(size =2)+
  geom_line()+
  #geom_vline(xintercept = 26, color = "palegreen4")+
  geom_point(aes(x = idx, y = max(gkl)+diff+diff/1.8*rep(c(noSig:1),each = noModels), col = dat2$sigty), shape = 15, size = 4)+
  #geom_text(datparam, mapping = aes(x = idx, y = y, label = param), size = 2.5)+
  geom_point(aes(x = idx, y = bic), col = "palegreen4")+
  geom_line(aes(x = idx, y = bic), col = "palegreen4")+
  #geom_text(aes(x = idx, y = max(gkl)+diff/1.2, label = rep(resMat[prm_order,"nprmtot"],4)), size = 3.5)+
  xlab("Possible interaction models")+
  ylab("Generalized Kullback-Leibler (GKL)")+
  #ggtitle("Fit of data for all possible interaction models on BRCA21")+
  theme_bw()+
  theme(legend.position = c(0.75,0.53), 
        legend.text.align = 0.5,
        legend.title = element_text(face = "bold"),
        legend.background = element_blank()
        
  )+
  scale_color_manual(values = col.model,
                     labels = c(expression(L + M + R), expression(L%*%M + M%*%R), expression(L%*%M%*%R)),
                     name = "Signature factorization") +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Generalized Kullback-Leibler (GKL)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~(.-addC)/mC, name="Bayesian Information Criteria (BIC)")
  ) +
  theme(
    axis.title.y = element_text(color = "black", size=13),
    axis.title.y.right = element_text(color = "palegreen4", size=13)
  )+ 
  scale_x_discrete(breaks = c(1:45))+
  geom_point(aes(x = which.min(dat2$bic), y = dat2$bic[which.min(dat2$bic)]), color = "palegreen4", pch = 1, size = 6)
p

colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")

# BRCA214 
yval1 = 33800
yval2 = 41250
dattext = data.frame(xp = c(1,18,26,45), yp = 42000, type = c("Mono", "Di", "Mix","Tri"))

p2 = p + geom_rect(aes(xmin = 0.6, xmax = 1.4, ymin = yval1, ymax = yval2), 
                   fill = "white", alpha = 0, color = colors3[1])+ 
  geom_rect(aes(xmin = 17.6, xmax = 18.4, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[2])+ 
  geom_rect(aes(xmin = 25.6, xmax = 26.4, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[3])+ 
  geom_rect(aes(xmin = 44.6, xmax = 45.4, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[4])+
  # geom_segment(aes(x = 6, y = 71000, xend = 6, yend = 75000),
  #              arrow = arrow(length = unit(0.28, "cm")))+
  # geom_text(data = data.frame(x = 6, y = 70000), 
  #           mapping = aes(x = x, y = y, label = "Total parameters"), 
  #           fontface = "italic", col = "black")+
  
  geom_text(dattext,mapping = aes(x = xp, y = yp, label = type), size = 3.5, fontface = "plain")

p2

# plot with number of instead of many of the same color 

model = matrix(0, nrow = 45, ncol = 3)
resMatnew = resMat[prm_order,]
for(i in 1:45){
  model[i,] = table(factor(resMatnew[i,2:9], levels = c(12,42,96)))
}

dat2$mono = rep(model[,1],8)
dat2$di = rep(model[,2],8)
dat2$full = rep(model[,3],8)

dat2



p = ggplot(dat2, aes(x = idx, y = gkl))+
  #geom_vline(xintercept = 16, lty = "dashed")+
  geom_point(size =3)+
  geom_line()+
  #geom_vline(xintercept = 26, color = "palegreen4")+
  #geom_point(aes(x = idx, y = max(gkl)+diff+diff/1.8*rep(c(1:noSig),each = noModels), col = dat2$sigty), shape = 15, size = 4)+
  #geom_text(datparam, mapping = aes(x = idx, y = y, label = param), size = 3.5)+
  geom_label(aes(x = id, y = 35000, label = full), fill = "#27408B", col = "white", fontface = "bold",cex = 3)+
  geom_label(aes(x = id, y = 36200, label = di), fill = "#FF9912", col = "white", fontface = "bold",cex = 3)+
  geom_label(aes(x = id, y = 37400, label = mono), fill = "#AC0136", col = "white", fontface = "bold",cex = 3)+
  geom_point(aes(x = idx, y = bic), col = "palegreen4")+
  geom_line(aes(x = idx, y = bic), col = "palegreen4")+
  #geom_text(aes(x = idx, y = max(gkl)+diff/1.2, label = rep(resMat[prm_order,"nprmtot"],4)), size = 3.5)+
  xlab("Possible interaction models")+
  ylab("Generalized Kullback-Leibler (GKL)")+
  #ggtitle("Fit of data for all possible interaction models on BRCA21")+
  theme_bw()+
  theme(legend.position = c(0.75,0.53), 
        legend.text.align = 0.5,
        legend.title = element_text(face = "bold"),
        legend.background = element_blank()
        
  )+
  geom_text(aes(x = 0, y = 35000, label = "Mono"))+
  scale_color_manual(values = col.model,
                     labels = c(expression(L + M + R), expression(L%*%M + M%*%R), expression(L%*%M%*%R)),
                     name = "Signature factorization") +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Generalized Kullback-Leibler (GKL)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~(.-addC)/mC, name="Bayesian Information Criteria (BIC)")
  ) +
  theme(
    axis.title.y = element_text(color = "black", size=13),
    axis.title.y.right = element_text(color = "palegreen4", size=13)
  )+ 
  geom_point(aes(x = which.min(dat2$bic), y = dat2$bic[which.min(dat2$bic)]), color = "palegreen4", pch = 1, size = 6 )

p

colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")

# BRCA214 
yval1 = 34200
yval2 = 41000
dattext = data.frame(xp = c(1,18,26,45), yp = 42000, type = c("Mono", "Di", "Mix","Tri"))

p2 = p + geom_rect(aes(xmin = 0.78, xmax = 1.22, ymin = yval1, ymax = yval2), 
                   fill = "white", alpha = 0, color = colors3[1])+ 
  geom_rect(aes(xmin = 6.78, xmax = 7.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[2])+ 
  geom_rect(aes(xmin = 11.78, xmax = 12.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[3])+ 
  geom_rect(aes(xmin = 14.78, xmax = 15.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[4])+
  # geom_segment(aes(x = 6, y = 71000, xend = 6, yend = 75000),
  #              arrow = arrow(length = unit(0.28, "cm")))+
  # geom_text(data = data.frame(x = 6, y = 70000), 
  #           mapping = aes(x = x, y = y, label = "Total parameters"), 
  #           fontface = "italic", col = "black")+
  
  geom_text(dattext,mapping = aes(x = xp, y = yp, label = type), size = 3.5, fontface = "plain")

p2



########## UCUT model plot ---------------------------------------------------------
source("UCUT/loadUCUTmodels.R")
load("UCUT/result/UCUTmodelsummary.RData")

noSig = 2

# parameter order
prm_order = order(resMat[,"nprmtot"])
prm_order2 = order(resMat[,"nprm1"],resMat[,"nprm2"])
prm_order2 = order(resMat[,4], decreasing = T)
gkl = resMat[,"GKL"]

noModels = nrow(resMat)

diff = diff(range(gkl))/10


dat1 = data.frame(idx = c(1:length(gkl)), gkl = gkl[prm_order2], sigty = resMat[prm_order2,c(1:noSig)])
dat2 = reshape(dat1, varying = colnames(dat1)[-c(1,2)], direction = "long")
dat2$sigty = factor(dat2$sigty, levels = c(18,48,66,96,102,1536))

#dat2$bic = 2*dat2$gkl + log(nrow(V)*ncol(V))*(resMat[prm_order2,"nprmtot"]) 
dat2$bic = 2*dat2$gkl + log(sum(V > 0.5))*(resMat[prm_order2,"nprmtot"])

mC = (max(dat2$gkl) - min(dat2$gkl))/(max(dat2$bic) - min(dat2$bic))
dat2$bic = dat2$bic* mC
addC = min(dat2$gkl) - min(dat2$bic)
dat2$bic = dat2$bic + addC

dat2$sigty2 = dat2$sigty
dat2$sigty2[8] = 96
dat2$sigty2[29] = 102

dat2$sigty2[5] = 96
dat2$sigty2[26] = 48

datparam = data.frame(idx = c(1:nrow(resMat)), param = resMat[prm_order2,"nprmtot"], y = max(gkl)+diff/1.2)



p = ggplot(dat2, aes(x = idx, y = gkl))+
  #geom_vline(xintercept = 16, lty = "dashed")+
  geom_point(size =3)+
  geom_line()+
  geom_point(aes(x = idx, y = max(gkl)+diff+diff/1.8*rep(c(1:noSig),each = noModels), col = dat2$sigty2), shape = 15, size = 4)+
  geom_text(datparam, mapping = aes(x = idx, y = y, label = param), size = 3.5)+
  geom_point(aes(x = idx, y = bic), col = "palegreen4")+
  geom_line(aes(x = idx, y = bic), col = "palegreen4")+
  xlab("Possible interaction models")+
  ylab("Generalized Kullback-Leibler (GKL)")+
  # geom_segment(aes(x = 18, y = 10200, xend = 18, yend = 10500),
  #              arrow = arrow(length = unit(0.28, "cm")))+
  # geom_text(data = data.frame(x = 18, y = 10000), mapping = aes(x = x, y = y, label = "Total parameters"), fontface = "italic", col = "black")+
  theme_bw()+
  theme(legend.position = c(0.15,0.39), 
        legend.text.align = 0.5,
        legend.title = element_text(face = "bold"),
        #legend.background = element_blank()
        #legend.title.align = 0.5
        
  )+
  scale_color_manual(values = brewer.pal(n = 6, name = "Dark2"),
                     labels = c(expression(L[2] + L[1] + M + R[1] + R[2]), 
                                expression(L[2] + L[1]%*%M + M%*%R[1] + R[2]), 
                                expression(L[2]%*%L[1] + L[1]%*%M + M%*%R[1] + R[1]%*%R[2]), 
                                expression(L[1]%*%M%*%R[1]), 
                                expression(L[2] + L[1]%*%M%*%R[1] + R[2]),
                                expression(L[2]%*%L[1]%*%M%*%R[1]%*%R[2])),
                     name = "Signature factorization")+
  geom_vline(xintercept = c(1.5,3.5,6.5,10.5,15.5), linetype = "dashed", color = "darkgray")+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Generalized Kullback-Leibler (GKL)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~(.-addC)/mC, name="Bayesian Information Criteria (BIC)")
  ) +
  theme(
    axis.title.y = element_text(color = "black", size=13),
    axis.title.y.right = element_text(color = "palegreen4", size=13)
  )+ 
  scale_x_discrete(breaks = c(1:15))+
  geom_point(aes(x = which.min(dat2$bic), y = dat2$bic[which.min(dat2$bic)]), color = "palegreen4", pch = 1, size = 6 )


p 

colors3 = c(brewer.pal(n = 6, name = "Dark2"),"#436EEE")


dattext = data.frame(xp = c(1,15,21), yp = 11400, type = c("Mono", "Di", "Penta"))
#yval1 = 79200
#yval2 = 89000
yval1 = 10860
yval2 = 11250
p + geom_rect(aes(xmin = 0.78, xmax = 1.22, ymin = yval1, ymax = yval2), 
              fill = "white", alpha = 0, color = colors3[1])+ 
  geom_rect(aes(xmin = 14.78, xmax = 15.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[3])+
  geom_rect(aes(xmin = 20.78, xmax = 21.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[6])+
  #geom_rect(aes(xmin = 10.78, xmax = 11.22, ymin = yval1, ymax = yval2), 
  #          fill = "white", alpha = 0, color = colors3[7])+
  geom_text(dattext,mapping = aes(x = xp, y = yp, label = type), size = 3.5, fontface = "plain")

###### UCUT other models ---------------------------------------------------------

source("UCUT/loadUCUTmodels.R")
load("UCUT/result/UCUTmodelsummary.RData")
resMatold = resMat
load("UCUT/result/UCUTditriMIXmodel.RData")
resMat2 = rbind(resMatold,resMat)
load("UCUT/result/UCUTditrimodel.RData")
resMat3 = rbind(resMat2,c(120,120,240,resditri$gkl))

resMat = resMat3

resMat = resMat[-c(8,9,13,23,24,25,27,11,15,18,6,20,21),]

resMat = resMatold[c(1,5,6,19,20,21),]

noSig = 2

# parameter order

prm_order2 = order(resMat[,4], decreasing = T)
gkl = resMat[,"GKL"]

noModels = nrow(resMat)

diff = diff(range(gkl))/10


dat1 = data.frame(idx = c(1:length(gkl)), gkl = gkl[prm_order2], sigty = resMat[prm_order2,c(1:noSig)])
dat2 = reshape(dat1, varying = colnames(dat1)[-c(1,2)], direction = "long")
#dat2$sigty = factor(dat2$sigty, levels = c(18,96,48,102,66,120))
dat2$sigty = factor(dat2$sigty, levels = c(18,66,1536))

#dat2$bic = 2*dat2$gkl + log(nrow(V)*ncol(V))*(resMat[prm_order2,"nprmtot"]) 
dat2$bic = 2*dat2$gkl + log(sum(V > 0.5))*(resMat[prm_order2,"nprmtot"])

mC = (max(dat2$gkl) - min(dat2$gkl))/(max(dat2$bic) - min(dat2$bic))
dat2$bic = dat2$bic* mC
addC = min(dat2$gkl) - min(dat2$bic)
dat2$bic = dat2$bic + addC

dat2$sigty2 = dat2$sigty
dat2$sigty2[13] = 18
dat2$sigty2[14] = 66
 
dat2$sigty2[28] = 120
dat2$sigty2[29] = 120

datparam = data.frame(idx = c(1:nrow(resMat)), param = resMat[prm_order2,"nprmtot"], y = max(gkl)+diff/1.2)



p = ggplot(dat2, aes(x = idx, y = gkl))+
  #geom_vline(xintercept = 16, lty = "dashed")+
  geom_point(size =3)+
  geom_line()+
  geom_point(aes(x = idx, y = max(gkl)+diff+diff/1.8*rep(c(1:noSig),each = noModels), col = dat2$sigty2), shape = 15, size = 4)+
  geom_text(datparam, mapping = aes(x = idx, y = y, label = param), size = 3.5)+
  geom_point(aes(x = idx, y = bic), col = "palegreen4")+
  geom_line(aes(x = idx, y = bic), col = "palegreen4")+
  xlab("Possible interaction models")+
  ylab("Generalized Kullback-Leibler (GKL)")+
  # geom_segment(aes(x = 18, y = 10200, xend = 18, yend = 10500),
  #              arrow = arrow(length = unit(0.28, "cm")))+
  # geom_text(data = data.frame(x = 18, y = 10000), mapping = aes(x = x, y = y, label = "Total parameters"), fontface = "italic", col = "black")+
  theme_bw()+
  theme(legend.position = c(0.23,0.28), 
        #legend.position = "bottom",
        legend.text.align = 0.5,
        legend.title = element_text(face = "bold"),
        text = element_text(size = 12)
        #legend.background = element_blank()
        #legend.title.align = 0.5
        
  )+
  scale_color_manual(values = c(brewer.pal(n = 7, name = "Dark2")[c(1,4,2,5,3,7,6)])[c(1,5,7)],
                     labels = c(expression(L[2] + L[1] + M + R[1] + R[2]), 
                                #expression(L[1]%*%M%*%R[1]),
                                #expression(L[2] + L[1]%*%M + M%*%R[1] + R[2]),
                                #expression(L[2] + L[1]%*%M%*%R[1] + R[2]),
                                expression(L[2]%*%L[1] + L[1]%*%M + M%*%R[1] + R[1]%*%R[2]),
                                #expression(L[2]%*%L[1] + L[1]%*%M%*%R[1] + R[1]%*%R[2]),
                                expression(L[2]%*%L[1]%*%M%*%R[1]%*%R[2])),
                     name = "Signature factorization")+
  #geom_vline(xintercept = c(1.5,3.5,5.5,7.5,12.5), linetype = "dashed", color = "darkgray")+
  scale_y_continuous(
    
    # Features of the first axis
    name = "Generalized Kullback-Leibler (GKL)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~(.-addC)/mC, name="Bayesian Information Criteria (BIC)")
  ) +
  theme(
    axis.title.y = element_text(color = "black", size=13),
    axis.title.y.right = element_text(color = "palegreen4", size=13)
  )+ 
  scale_x_discrete(breaks = c(1:15))+
  geom_point(aes(x = which.min(dat2$bic), y = dat2$bic[which.min(dat2$bic)]), color = "palegreen4", pch = 1, size = 6 )


p 

colors3 = c(brewer.pal(n = 6, name = "Dark2"),"#436EEE")


dattext = data.frame(xp = c(1,12,14), yp = 10920, type = c("Mono", "Di", "Mix"))

yval1 = 10650
yval2 = 10850
p + geom_rect(aes(xmin = 0.78, xmax = 1.22, ymin = yval1, ymax = yval2), 
              fill = "white", alpha = 0, color = colors3[1])+ 
  geom_rect(aes(xmin = 11.78, xmax = 12.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[3])+
  #geom_rect(aes(xmin = 17.78, xmax = 18.22, ymin = yval1, ymax = yval2), 
  #          fill = "white", alpha = 0, color = colors3[6])+
  geom_rect(aes(xmin = 13.78, xmax = 14.22, ymin = yval1, ymax = yval2), 
           fill = "white", alpha = 0, color = colors3[7])+
  geom_text(dattext,mapping = aes(x = xp, y = yp, label = type), size = 3.5, fontface = "plain")
  #geom_text(data.frame(x = 13.9, y = 10300, label = c("prm")),mapping = aes(x = x, y = y, label = label), fontface = "bold")+
  #geom_text(data.frame(x = 13.9, y = 9800, label = c("18 \n 96 \n 48 \n 102 \n 66 \n 120")),mapping = aes(x = x, y = y, label = label))

dattext = data.frame(xp = c(1,3,6), yp = 11400, type = c("Mono", "Di", "Penta"))
yval1 = 10830
yval2 = 11280

p + geom_rect(aes(xmin = 0.88, xmax = 1.12, ymin = yval1, ymax = yval2), 
              fill = "white", alpha = 0, color = colors3[1])+ 
  geom_rect(aes(xmin = 2.88, xmax = 3.12, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[3])+
  #geom_rect(aes(xmin = 17.78, xmax = 18.22, ymin = yval1, ymax = yval2), 
  #          fill = "white", alpha = 0, color = colors3[6])+
  geom_rect(aes(xmin = 5.88, xmax = 6.12, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[6])+
  geom_text(dattext,mapping = aes(x = xp, y = yp, label = type), size = 3.5, fontface = "plain")
##################### BIC plot ---------------------------------------------------- 
#BIC = 2*resMat[,"GKL"] + log(nrow(V)*ncol(V))*resMat[,"nprmtot"]*8
BIC = 2*resMat[,"GKL"] + log(nrow(V) + ncol(V))*resMat[,"nprmtot"]*noSig
plot(BIC, ylim = c(range(BIC)[1],range(BIC)[2]+10000), main = "BIC")
for(i in 1:noSig){
  points(rep(max(BIC)+10000 - (i-1)*1000,noModels), pch = 15, col = resMatnew[,i])
}

abline(v = which.min(BIC), lty = "dashed")
abline(v = 9, lty = "dashed")

# gkl
plot(gkl, ylim = c(range(gkl)[1],range(gkl)[2]+6000), main = "GKL", xlab = "Model number")
for(i in 1:noSig){
  points(rep(max(gkl) + i*700,noModels), pch = 15, col = resMatnew[,i])
}

# gkl ordered
gkl_order = gkl[prm_order]
plot(gkl_order, ylim = c(range(gkl_order)[1],range(gkl)[2]+6000), main = "Model fit for BRCA214", xlab = "Model number")
for(i in 1:noSig){
  points(rep(max(gkl) + i*700,noModels), pch = 15, col = resMatnew[prm_order,i])
}
abline(v = 20.5, lty = "dashed")
abline(v = 16.5, lty = "dashed")
abline(v = 7, lty = "dashed")



  




