##-------------------------------------------------
## Histogram of the number of non-zero entries 
##-------------------------------------------------
## Load UCUT data
load("UCUT_5_all.RData")
V5 = Vall
##-------------------------------------------------------
## Number of non-zero data points for each patient
n.pat <- dim(V5)[2]
n.nonz <- rep(0,n.pat)
for (i in 1:n.pat){
  n.nonz[i] <- length( which(V5[,i]>0.5) )
}
print(sort(n.nonz))
## 4^4*6=1536
brks <- 50*(0:14)
#pdf(file="UCUT26hist.pdf",width=8,height=6)
hist(n.nonz,breaks=brks,main="",
     xlab="Number of non-zero entries for UCUT patient",freq=TRUE,axes=FALSE,
     ylim=c(0,10),cex.lab=1.2)
axis(1,at=brks,labels=brks,las=1,cex.axis=1.2)
axis(2,at=0:10,labels=0:10,las=1,cex.axis=1.2)
text(200,8,expression(paste("Number of mutation types is ",4^4%.%6," = ",1536,sep="")),
     pos=4,cex=1.2)
text(200,7,"Total number of non-zero entries is 5260",pos=4,cex=1.2)
#dev.off()