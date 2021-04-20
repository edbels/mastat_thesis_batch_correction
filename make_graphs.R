# make boxplots over iterations of performance metrics

rm(list = ls())


### choose option of split 
# -> 1: equal split
# -> 2: unequal split
# -> 3: equal split, but with all "Naive CD4 T" in batch one
setwd("C:/Users/edbels/Documents/GitHub/mastat_thesis")

testset <- "pbmc"
tot_iterations <- 1


mmd <- readRDS(paste0("results/",testset,"/mmd.rds"))
m_mmd <- apply(data.frame(mmd),1, function(x) min(x[x>0]))
testvalues <- c(0, 0.1,0.25,0.5, 1.1,1.25,1.5, which.min(m_mmd),
                which.min(abs(m_mmd - quantile(m_mmd,0.25))),
                which.min(abs(m_mmd - median(m_mmd))),
                which.min(abs(m_mmd - quantile(m_mmd,0.75))),
                which.max(m_mmd))

# ARI1
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- ari[[1]][1]
  seurat[i] <- ari[[2]][1]
  harmony[i] <- ari[[3]][1]
  sctransform[i] <- ari[[4]][1]
  scMerge[i] <- ari[[5]][1]
  liger[i] <- ari[[6]][1]
  i <- i+1
}
fname <- "ARIbatch"
yname <- "ARI(batch)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(-0.001,0.005), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.005, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
         col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()
  
jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(-0.001,0.5), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,0.5, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()
  
jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(-0.001,0.001), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.001, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()



# ARI2
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- ari[[1]][2]
  seurat[i] <- ari[[2]][2]
  harmony[i] <- ari[[3]][2]
  sctransform[i] <- ari[[4]][2]
  scMerge[i] <- ari[[5]][2]
  liger[i] <- ari[[6]][2]
  i <- i+1
}
fname <- "ARIct"
yname <- "ARI(celltype)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(0,1.1), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.3, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,1.1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.1,1, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0,1.1), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.35, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()


# ARI2
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- ari[[1]][2]
  seurat[i] <- ari[[2]][2]
  harmony[i] <- ari[[3]][2]
  sctransform[i] <- ari[[4]][2]
  scMerge[i] <- ari[[5]][2]
  liger[i] <- ari[[6]][2]
  i <- i+1
}
fname <- "ARIct"
yname <- "ARI(celltype)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(0,1.1), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.3, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,1.1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.1,1, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0,1.1), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.35, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()



# NMI1
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- nmi[[1]][1]
  seurat[i] <- nmi[[2]][1]
  harmony[i] <- nmi[[3]][1]
  sctransform[i] <- nmi[[4]][1]
  scMerge[i] <- nmi[[5]][1]
  liger[i] <- nmi[[6]][1]
  i <- i+1
}
fname <- "NMIbatch"
yname <- "NMI(batch)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(0,0.005), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.004, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,0.5), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,0.4, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0,0.002), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.002, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()


# NMI2
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- nmi[[1]][2]
  seurat[i] <- nmi[[2]][2]
  harmony[i] <- nmi[[3]][2]
  sctransform[i] <- nmi[[4]][2]
  scMerge[i] <- nmi[[5]][2]
  liger[i] <- nmi[[6]][2]
  i <- i+1
}
fname <- "NMIct"
yname <- "NMI(celltype)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(0,1.1), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.4, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,1.1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.1,0.4, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0,1.1), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.35, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

# ASW1
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- asw[[1]][1]
  seurat[i] <- asw[[2]][1]
  harmony[i] <- asw[[3]][1]
  sctransform[i] <- asw[[4]][1]
  scMerge[i] <- asw[[5]][1]
  liger[i] <- asw[[6]][1]
  i <- i+1
}
fname <- "aswbatch"
yname <- "asw(batch)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(-0.05,0.015), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,-0.03, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,0.25), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,0.2, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(-0.005,0.005), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,-0.002, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()


# ASW2
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- asw[[1]][2]
  seurat[i] <- asw[[2]][2]
  harmony[i] <- asw[[3]][2]
  sctransform[i] <- asw[[4]][2]
  scMerge[i] <- asw[[5]][2]
  liger[i] <- asw[[6]][2]
  i <- i+1
}
fname <- "aswct"
yname <- "asw(celltype)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(0.2,0.8), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.7, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,1.1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,0.9, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0.3,0.5), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.38, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()


# LISI1
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- lisi[[1]][1]/lisi[[1]][3]
  seurat[i] <- lisi[[2]][1]/lisi[[2]][3]
  harmony[i] <- lisi[[3]][1]/lisi[[3]][3]
  sctransform[i] <- lisi[[4]][1]/lisi[[4]][3]
  scMerge[i] <- lisi[[5]][1]/lisi[[5]][3]
  liger[i] <- lisi[[6]][1]/lisi[[6]][3]
  i <- i+1
}
fname <- "LISIbatch"
yname <- "LISI(batch)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(0.8,1.2), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,1.15, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0.5,1.1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,0.7, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0.9,1.05), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,1.04, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()


# LISI2
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- lisi[[1]][2]
  seurat[i] <- lisi[[2]][2]
  harmony[i] <- lisi[[3]][2]
  sctransform[i] <- lisi[[4]][2]
  scMerge[i] <- lisi[[5]][2]
  liger[i] <- lisi[[6]][2]
  i <- i+1
}
fname <- "LISIct"
yname <- "LISI(celltype)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(0.8,1.5), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,1, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,2.4), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.1,0.8, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0.8,1.5), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,1, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

# kbet05
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- kbet[[1]][1]
  seurat[i] <- kbet[[2]][1]
  harmony[i] <- kbet[[3]][1]
  sctransform[i] <- kbet[[4]][1]
  scMerge[i] <- kbet[[5]][1]
  liger[i] <- kbet[[6]][1]
  i <- i+1
}
fname <- "kbet05"
yname <- "kBet (5%)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(-0.05,0.5), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.28, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,1, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0,0.9), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.8, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

# kbet1
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- kbet[[1]][2]
  seurat[i] <- kbet[[2]][2]
  harmony[i] <- kbet[[3]][2]
  sctransform[i] <- kbet[[4]][2]
  scMerge[i] <- kbet[[5]][2]
  liger[i] <- kbet[[6]][2]
  i <- i+1
}
fname <- "kbet1"
yname <- "kBet (10%)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(-0.05,0.5), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.28, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,1, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0,0.9), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.8, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

# kbet20
merge <- rep(0,length(testvalues))
seurat <- rep(0,length(testvalues))
harmony <- rep(0,length(testvalues))
sctransform <- rep(0,length(testvalues))
scMerge <- rep(0,length(testvalues))
liger <- rep(0,length(testvalues))
i <- 1
for (sb in testvalues) {
  fn <- paste0("results/",testset,"/",sb,"_it",tot_iterations,".RData")
  load(fn)
  merge[i] <- kbet[[1]][3]
  seurat[i] <- kbet[[2]][3]
  harmony[i] <- kbet[[3]][3]
  sctransform[i] <- kbet[[4]][3]
  scMerge[i] <- kbet[[5]][3]
  liger[i] <- kbet[[6]][3]
  i <- i+1
}
fname <- "kbet20"
yname <- "kBet (20%)"
jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(-0.05,0.6), xlab = "fraction in group 1", ylab = yname)
lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
legend(0.4,0.55, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(0,1), xlab = "sampling depth", ylab = yname)
lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
legend(0.8,1, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()

jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(0,0.6), xlab = "distance", ylab = yname)
lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
legend(0.8,0.5, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
       col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
dev.off()





















mg(m,s,h,st,sm,l,"ARIbatch")
mg <- function(merge,seurat,harmony,sctransform,scMerge,liger, fname){
  jpeg(paste0("results/",testset,"/",fname,"_unequalSplit.jpeg"))
  plot(c(0.1,0.25,0.5),merge[c(2,3,4)],'l',ylim=c(-0.001,0.005), xlab = "fraction in group 1", ylab = "ARI(batch)")
  lines(c(0.1,0.25,0.5),seurat[c(2,3,4)],col="blue",lty =2)
  lines(c(0.1,0.25,0.5),harmony[c(2,3,4)],col="red",lty = 3)
  lines(c(0.1,0.25,0.5),sctransform[c(2,3,4)],col="green")
  lines(c(0.1,0.25,0.5),scMerge[c(2,3,4)],col="purple", lty = 2)
  lines(c(0.1,0.25,0.5),liger[c(2,3,4)],col="brown", lty = 3)
  legend(0.8,0.005, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
         col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
  dev.off()
  
  jpeg(paste0("results/",testset,"/",fname,"_Depth.jpeg"))
  plot(c(0.1,0.25,0.5,1),merge[c(5,6,7,1)],'l',ylim=c(-0.001,0.5), xlab = "sampling depth", ylab = "ARI(batch)")
  lines(c(0.1,0.25,0.5,1),seurat[c(5,6,7,1)],col="blue",lty =2)
  lines(c(0.1,0.25,0.5,1),harmony[c(5,6,7,1)],col="red",lty =3)
  lines(c(0.1,0.25,0.5,1),sctransform[c(5,6,7,1)],col="green")
  lines(c(0.1,0.25,0.5,1),scMerge[c(5,6,7,1)],col="purple",lty =2)
  lines(c(0.1,0.25,0.5,1),liger[c(5,6,7,1)],col="brown",lty =3)
  legend(0.8,0.5, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
         col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
  dev.off()
  
  jpeg(paste0("results/",testset,"/",fname,"_SeparateCluster.jpeg"))
  plot(c(0,0.25,0.5,0.75,1),merge[c(8,9,10,11,12)],'l',ylim=c(-0.001,0.001), xlab = "distance", ylab = "ARI(batch)")
  lines(c(0,0.25,0.5,0.75,1),seurat[c(8,9,10,11,12)],col="blue",lty =2)
  lines(c(0,0.25,0.5,0.75,1),harmony[c(8,9,10,11,12)],col="red",lty =3)
  lines(c(0,0.25,0.5,0.75,1),sctransform[c(8,9,10,11,12)],col="green")
  lines(c(0,0.25,0.5,0.75,1),scMerge[c(8,9,10,11,12)],col="purple",lty =2)
  lines(c(0,0.25,0.5,0.75,1),liger[c(8,9,10,11,12)],col="brown",lty =3)
  legend(0.8,0.5, legend=c("Merge", "Seurat v3","harmony","scTransform","scMerge","liger"),
         col=c("black", "blue","red","green","purple","brown"), lty=1:3, cex=0.8)
}

