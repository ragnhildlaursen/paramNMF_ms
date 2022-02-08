##-----------------------------------------
## UCUT26: Table and plots of GKL and BIC
##-----------------------------------------
## Read data: 
## Model for the two signatures (given by number of parameters) 
## and corresponding GKL
resMat <- matrix( scan(file="GitUCUT26results.txt"),ncol=4,byrow=21 )
colnames(resMat) <- c("nprm1","nprm2","nprmtot","GKL")
## Parameters and corresponding factor models
prm.str <- c(18,48,66,78,102,96)
col.str <- c("orange","skyblue","blue","brown","green","red")
mdl.str <- c(expression(L[2]+L[1]+M+R[1]+R[2]), 
             expression(L[2]+L[1]~x~M+M~x~R[1]+R[2]),
             expression(L[2]~x~L[1]+L[1]~x~M+M~x~R[1]+R[1]~x~R[2]),
             expression(L[2]~x~M+L[1]~x~M+M~x~R[1]+M~x~R[2]),
             expression(L[2]+L[1]~x~M~x~R[1]+R[2]),
             expression(L[1]~x~M~x~R[1]))
##------------------------------------------------------
## Table with summary statistics for all 21 models
##------------------------------------------------------
resMat <- cbind(resMat,rep(0,21),rep(0,21),rep(0,21))
colnames(resMat)[5] <- "Complexity"  
colnames(resMat)[6] <- "rawBIC"
colnames(resMat)[7] <- "diffBIC"   
number.nonz <- 5260
resMat[,"Complexity"] <- resMat[,"nprmtot"]*log( number.nonz )
resMat[,"rawBIC"] <- resMat[,"Complexity"]+2*resMat[,"GKL"]
resMat[,"diffBIC"] <- resMat[,"rawBIC"]-min(resMat[,"rawBIC"])
print(resMat[,c("nprm1","nprm2","Complexity","GKL","rawBIC","diffBIC")])
##-------------------------------------------------------
## Table with summary statistics for the six models
## where the two signatures have the same parametrization
##-------------------------------------------------------
sameIndx <- which(resMat[,"nprm1"]==resMat[,"nprm2"])
sameRes <- resMat[sameIndx,]
indx <- c(5,1,2,4,6,3)
sameRes <- sameRes[indx,]
print(sameRes)
##---------------------------------------------------
## Plot the GKL for the six models where the two signatures 
## have the same parametrization
## The plot is saved in the file UCUT26GKLsame.pdf 
##------------------------------------------------
#pdf(file="UCUT26GKLsame.pdf",width=8,height=6)
plot(1:6,sameRes[,"GKL"],ylim=(range(sameRes[indx,"GKL"])+c(0,800)),
     xlab="Model",ylab="Generalized Kullback-Leibler (GKL)",pch=19,cex.lab=1.3,cex.axis=1.3)
points(1:6,sameRes[,"GKL"],type="l")
mx <- max(resMat[,"GKL"])+800
for (i in indx){
  points(i,mx,pch=15,col=col.str[match(sameRes[i,"nprm1"],prm.str)],cex=1.5)
  points(i,mx-200,pch=15,col=col.str[match(sameRes[i,"nprm2"],prm.str)],cex=1.5)
  text(i,mx-400,sameRes[i,"nprmtot"])
}
## Add factor and number of parameters to plot
yval <- mx-600
text(2.8,yval,"Parametrization",pos=4)  
text(5.5,yval+70,"Number of",pos=NULL) ; text(5.5,yval,"parameters",pos=NULL)
j <- 0
for (i in c(1,6,2,4,5,3)){
  j <- j+1
  points(2.8,yval-100*j,pch=15,col=col.str[i],cex=1.5)
  text(3,yval-100*j,mdl.str[i],pos=4)
  text(5.5,yval-100*j,prm.str[i],pos=2)
}
#dev.off()
##-------------------------------------------------------------------------
## GKL plot of the 21 models
##---------------------------
## Index of plotting
indx <-     c(1,6,2,4,8,3, 19,11,17,20,15, 5,9,13,7, 14,18,12, 21,16, 10)
prm.indx <- c(1,1,1,1,1,1,  1, 2, 2, 2, 1, 1,2, 1,1,  1, 1, 1, 1, 1,  1)
mdl.indx <- c(1:6,8:12,14:17,19:21,23,24,26)
#pdf(file="UCUT26GKLnew.pdf",width=10,height=5)
plot(mdl.indx,resMat[indx,"GKL"],ylim=(range(resMat[,"GKL"])+c(0,800)),
     xlab="Model",ylab="Generalized Kullback-Leibler (GKL)",
     pch=19,cex.lab=1.2,cex.axis=1.2,xaxt='n')
points(mdl.indx,resMat[indx,"GKL"],type="l")
abline(v=c(7,13,18,22,25),col="gray")
axis(side=1,at=mdl.indx,labels=1:21,cex.axis=1.2)
mx <- max(resMatr[,"GKL"])+800
j <- 0
for (i in indx){
  j <- j+1
  points(mdl.indx[j],ifelse(prm.indx[j]==1,mx,mx-200),pch=15,col=col.str[match(resMat[i,"nprm1"],prm.str)],cex=1.5)
  points(mdl.indx[j],ifelse(prm.indx[j]==1,mx-200,mx),pch=15,col=col.str[match(resMat[i,"nprm2"],prm.str)],cex=1.5)
  text(mdl.indx[j],mx-400,resMat[i,"nprmtot"],cex=0.9)
}
## Text concerning poor fit
text(3,10600,"Poor fit to data",pos=4)
arrows(4,10500,1.3,10420,length=0.15,angle=20,lwd=2)
arrows(5,10500,7.7,10180,length=0.15,angle=20,lwd=2)
## Add factor and number of parameters to plot
xval <- 16 ; yval <- 10500
text(xval,yval,"Parametrization of signature",pos=4)  
text(xval+8,yval+90,"Number of",pos=4) ; text(xval+8,yval,"parameters",pos=4)
for (i in 1:6){
  points(xval-0.5,yval-120*i,pch=15,col=col.str[i],cex=1.5)
  text(xval,yval-120*i,mdl.str[i],pos=4)
  text(xval+10,yval-120*i,prm.str[i],pos=2)
}
#dev.off()
##--------------------------------------
## BIC plot of the 21 models
##--------------------------------------
#pdf(file="UCUT26BICnew.pdf",width=10,height=5)
plot(mdl.indx,resMat[indx,"rawBIC"],ylim=(range(resMat[,"rawBIC"])+c(0,800)),
     xlab="Model",ylab="Bayesian Information Criterion (BIC)",
     pch=19,cex.lab=1.2,cex.axis=1.2,xaxt='n')
points(mdl.indx,resMat[indx,"rawBIC"],type="l")
abline(v=c(7,13,18,22,25),col="gray")
axis(side=1,at=mdl.indx,labels=1:21,cex.axis=1.2)
mx <- max(resMat[,"rawBIC"])+800
j <- 0
for (i in indx){
  j <- j+1
  points(mdl.indx[j],ifelse(prm.indx[j]==1,mx,mx-200),pch=15,col=col.str[match(resMat[i,"nprm1"],prm.str)],cex=1.5)
  points(mdl.indx[j],ifelse(prm.indx[j]==1,mx-200,mx),pch=15,col=col.str[match(resMat[i,"nprm2"],prm.str)],cex=1.5)
  text(mdl.indx[j],mx-400,resMat[i,"nprmtot"],cex=0.9)
}
#abline(h=19650,lwd=2,col="black",lty="dotted")
#text(1,19500,"Appropriate models below dotted line",pos=4)
text(26,19147,pos=3,"1") ; text(17,19264,pos=3,"2") ; text(21,19353,pos=3,"3")
text(24,19573,pos=3,"4") ; text(12,19742,pos=3,"5") ; text(6,20060,pos=3,"6")
text(10,19700,pos=2,"BIC ranking")
arrows(10,19700,11.7,19900,length=0.15,angle=20,lwd=2)
#dev.off()

