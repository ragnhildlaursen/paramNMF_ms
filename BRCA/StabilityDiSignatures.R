#######################################################
## Stability of signatures
#######################################################

setwd("~/projects/paramNMF_ms")
source("BRCA/loadBRCAmodels.R")

init = 20   # number of initialisations
tol = 0.1

DiRes <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tol,Seeds=sample(1:1000,init))

res1 = DiRes$AllResults[[2]][[7]] # 2059
res2 = DiRes$AllResults[[2]][[12]] # 2087
res3 = DiRes$AllResults[[2]][[13]] # 2116
res4 = DiRes$AllResults[[2]][[15]] # 2120
res5 = DiRes$AllResults[[2]][[9]] # 2125
cosMatch(res5,res3)
plot(DiRes$intermed[-c(1:1000),469])
par(mfrow = c(4,4))
for(i in c(401:416)){
  plot(DiRes$intermed[-c(1:100),i])
}

par(mfrow = c(1,1))
DiRes2 <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tol,Seeds=c(25:35))

res21 = DiRes2$AllResults[[2]][[3]] # 2088
res22 = DiRes2$AllResults[[2]][[8]] # 2109
cosMatch(res1,res21)

DiRes3 <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = rep(list(Mdi),noSig),tol=tol,Seeds=c(36:50))

res31 = DiRes3$AllResults[[2]][[2]] # 2122
res32 = DiRes3$AllResults[[2]][[5]] # 2134
res33 = DiRes3$AllResults[[2]][[13]] # 2107
cosMatch(res22,res33)

# now changed to 5000 max iterations
DiRes4 <- NMFglmSQR(Data = V, NoSignatures = 4, DesignMatrix = rep(list(Mdi),4),tol=tol,Seeds=c(191:200))

res41 = DiRes4$AllResults[[2]][[7]] # 2116
res42 = DiRes4$AllResults[[2]][[9]] # 2140
res43 = DiRes4$AllResults[[2]][[12]] # 2113
res44 = DiRes4$AllResults[[2]][[17]] # 2124
res45 = DiRes4$AllResults[[2]][[20]] # 2107
cosMatch(res41,res43)$cossim
gklval = DiRes4$AllResults[[3]]
gklval[which(gklval<2140)]
for(i in which(gklval<2150)){
  print(cosMatch(DiRes4$Signatures,DiRes4$AllResults[[2]][[i]])$cossim)
}
cosMatch(res44,res42)

par(mfrow = c(1,1))
plot(DiRes4$intermed[-c(1:1000),469])
par(mfrow = c(4,4))
for(i in c(301:316)){
  plot(DiRes4$intermed[-c(1:100),i])
}

# Try Mix Result
MixRes4 <- NMFglmSQR(Data = V, NoSignatures = noSig, DesignMatrix = list(Mmono,Mdi,Mdi,Mfull),tol=tol,Seeds=c(53:70))

cosMatch(MixRes4$Signatures,MixRes4$AllResults[[2]][[5]])$cossim
gklval = MixRes4$AllResults[[3]]
gklval[which(gklval<1870)]
for(i in which(gklval<1870)){
print(cosMatch(MixRes4$AllResults[[2]][[3]],MixRes4$AllResults[[2]][[i]])$cossim)
}


