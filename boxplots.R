setwd("C:/Users/edbels/Documents/GitHub/mastat_thesis")

load("cbmc_stability_1.RData")

methods <- c("Merge","Seurat","Harmony","scTransform","scMerge")

par(mfrow=c(2,3))
for (i in 1:5){
  boxplot(ari[[i]], main = methods[i], names = c("batch","cluster"))
}


dev.off()
par(mfrow=c(2,3))
for (i in 1:5){
  boxplot(asw[[i]], main = methods[i], names = c("batch","cluster"))
}

dev.off()
par(mfrow=c(1,3))
for (i in c(1,2,4)){
  boxplot(hvg[[i]], main = methods[i])
}

dev.off()
par(mfrow=c(2,3))
for (i in 1:5){
  boxplot(kbet[[i]], main = methods[i])
}

dev.off()
par(mfrow=c(2,3))
for (i in 1:5){
  lisi[[i]][,1] <- lisi[[i]][,1]/lisi[[i]][,3]
  boxplot(lisi[[i]][,1:2], main = methods[i], names = c("batch","cluster"))
}

dev.off()
par(mfrow=c(2,3))
for (i in 1:5){
  boxplot(nmi[[i]], main = methods[i], names = c("batch","cluster"))
}
