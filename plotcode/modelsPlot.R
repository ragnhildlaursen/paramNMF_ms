rm(list=ls())
#########################################
## Plot of models for BRCA 214 
#########################################
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(grid)
############## BRCA model plot ------------------------------------------------------

setwd("~/projects/paramNMF_ms/")
source("BRCA/loadBRCA214models.R")

load("BRCA/result/BRCA214modelsummary4sig500init.RData")


#load("UCUT/UCUTmodelsummary.RData")

noSig = 4
# colors for signatures
col.model = c("#AC0136", "#FF9912", "#27408B")

resMatnew = resMat
resMatnew[resMatnew == 12] = col.model[1]
resMatnew[resMatnew == 42] = col.model[2]
resMatnew[resMatnew == 96] = col.model[3]

# parameter order
prm_order = order(resMat[,"nprmtot"])

gkl = resMat[,"GKL"]

noModels = nrow(resMat)


diff = diff(range(gkl))/10
## ggplot
dat1 = data.frame(idx = c(1:length(gkl)), gkl = gkl[prm_order], sigty = resMat[prm_order,c(1:noSig)])
dat2 = reshape(dat1, varying = colnames(dat1)[-c(1,2)], direction = "long")
dat2$sigty = factor(dat2$sigty)

dat2$bic = 2*dat2$gkl + log(nrow(V)*ncol(V))*(resMat[prm_order,"nprmtot"]) 


mC = (max(dat2$gkl) - min(dat2$gkl))/(max(dat2$bic) - min(dat2$bic))
dat2$bic = dat2$bic* mC
addC = min(dat2$gkl) - min(dat2$bic)
dat2$bic = dat2$bic + addC
datparam = data.frame(idx = c(1:nrow(resMat)), param = resMat[prm_order,"nprmtot"], y = max(gkl)+diff/1.2)
p = ggplot(dat2, aes(x = idx, y = gkl))+
  #geom_vline(xintercept = 16, lty = "dashed")+
  geom_point(size =3)+
  geom_line()+
  geom_point(aes(x = idx, y = max(gkl)+diff+diff/1.8*rep(c(1:noSig),each = noModels), col = dat2$sigty), shape = 15, size = 4)+
  geom_text(datparam, mapping = aes(x = idx, y = y, label = param), size = 3.5)+
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
  scale_x_discrete(breaks = c(1:15))+
  geom_point(aes(x = which.min(dat2$bic), y = dat2$bic[which.min(dat2$bic)]), color = "palegreen4", pch = 1, size = 6 )
p

colors3 = c("#AC0136","#A895CD", "#FF9912", "#27408B")
colors3 = c("#AC0136","#FF9912", "#436EEE", "#27408B")

# BRCA21
yval1 = 5050
yval2 = 5770
dattext = data.frame(xp = c(1,7,8,15), yp = 5880, type = c("Mono", "Di", "Mix","Tri"))

p + geom_rect(aes(xmin = 0.78, xmax = 1.22, ymin = yval1, ymax = yval2), 
                fill = "white", alpha = 0, color = colors3[1])+ 
  geom_rect(aes(xmin = 6.78, xmax = 7.22, ymin = yval1, ymax = yval2), 
                fill = "white", alpha = 0, color = colors3[2])+ 
  geom_rect(aes(xmin = 7.78, xmax = 8.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[3])+ 
  geom_rect(aes(xmin = 14.78, xmax = 15.22, ymin = yval1, ymax = yval2), 
            fill = "white", alpha = 0, color = colors3[4])+
  geom_segment(aes(x = 6, y = 4400, xend = 6, yend = 4700),
               arrow = arrow(length = unit(0.28, "cm")))+
  geom_text(data = data.frame(x = 6, y = 4300), mapping = aes(x = x, y = y, label = "Total parameters"), fontface = "italic", col = "black")+
geom_text(dattext,mapping = aes(x = xp, y = yp, label = type), size = 3.5, fontface = "plain")
  


# BRCA214 
yval1 = 79200
yval2 = 89000
dattext = data.frame(xp = c(1,7,12,15), yp = 91000, type = c("Mono", "Di", "Mix","Tri"))

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

#### adding BIC values
p2 + scale_y_continuous(
  
  # Features of the first axis
  name = "Generalized Kullback Leibler",
  
  # Add a second axis and specify its features
  sec.axis = sec_axis(~.+ 10000, name="BIC values")
)
################################### UCUT model plot ---------------------------------------------------------
source("UCUT/loadUCUTmodels.R")
load("UCUT/UCUTmodelsummary.RData")

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

##



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



  



