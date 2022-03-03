rm(list=ls())
#########################################
## Plot of models for BRCA 214 
#########################################
setwd("~/projects/paramNMF_ms/")
source("BRCA/loadBRCA214models.R")
library(RColorBrewer)

# colors for signatures
col.model = c("#AC0136", "#FF9912", "#27408B")
resMatnew = resMat
resMatnew[resMatnew == 12] = col.model[1]
resMatnew[resMatnew == 42] = col.model[2]
resMatnew[resMatnew == 96] = col.model[3]

# parameter order
prm_order = order(resMat[,"nprmtot"])

# bic 
noModels = nrow(resMat)
#BIC = 2*resMat[,"GKL"] + log(nrow(V)*ncol(V))*resMat[,"nprmtot"]*8
BIC = 2*resMat[,"GKL"] + log(nrow(V))*resMat[,"nprmtot"]*8
plot(BIC, ylim = c(range(BIC)[1],range(BIC)[2]+10000), main = "BIC")
for(i in 1:8){
  points(rep(max(BIC)+10000 - (i-1)*1000,noModels), pch = 15, col = resMatnew[,i])
}

abline(v = which.min(BIC), lty = "dashed")
abline(v = 9, lty = "dashed")

# gkl
gkl = resMat[,"GKL"]
plot(gkl, ylim = c(range(gkl)[1],range(gkl)[2]+6000), main = "GKL", xlab = "Model number")
for(i in 1:8){
  points(rep(38000 - (i-1)*700,noModels), pch = 15, col = resMatnew[,i])
}
# gkl ordered
gkl_order = gkl[prm_order]
plot(gkl_order, ylim = c(range(gkl_order)[1],range(gkl)[2]+6000), main = "Model fit for BRCA214", xlab = "Model number")
for(i in 1:8){
  points(rep(38000 - (i-1)*700,noModels), pch = 15, col = resMatnew[prm_order,i])
}
abline(v = 20.5, lty = "dashed")
abline(v = 16.5, lty = "dashed")
abline(v = 7, lty = "dashed")

## ggplot
dat1 = data.frame(idx = c(1:length(gkl)), gkl = gkl[prm_order], sigty = resMat[prm_order,c(1:8)])
dat2 = reshape(dat1, varying = colnames(dat1)[-c(1,2)], direction = "long")
dat2$sigty = factor(dat2$sigty)
ggplot(dat2, aes(x = idx, y = gkl))+
  geom_vline(xintercept = 16, lty = "dashed")+
  geom_point()+
  geom_point(aes(x = idx, y = max(gkl)+500+500*rep(c(1:8),each = 45), col = dat2$sigty), shape = 15)+
  xlab("Model")+
  ylab("Generalized Kullback-Leibler (GKL)")+
  ggtitle("Fit of data for all possible interaction models on BRCA214")+
  #theme_bw()+
  theme(legend.position = c(0.7,0.5))+
  scale_color_manual(values = col.model, 
                     labels = c("Mono-nucleotide", "Di-nucleotide", "Tri-nucleotide"), 
                     name = "Signature parametrization")

  



# legend for models 
legend(x = c(2,6.5), y = c(5600,4200),  fill = col.model[c(1,2,3)], 
       legend = c("Mono-nucleotide","Di-nucleotide","Tri-nucleotide"), cex = 0.8, box.lty = 0)

# legend for GKL and BIC
legend(x = c(3.5,6.5), y = c(3300,2000), pch = rep(16,4), col = col.vec[c(1,3)], 
       legend = c("GKL","BIC"), cex = 0.8, box.lty = 0)
#dev.off()