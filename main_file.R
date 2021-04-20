# the code below does the same as h0_cbmc, but n times -> to check the stability of performance metrics

rm(list = ls())


# install.packages("devtools")
library(devtools)
# install.packages("doSNOW")
# install.packages("BiocManager")
# library(BiocManager)
# devtools::install_github('theislab/kbet',force = TRUE)
# devtools::install_github("immunogenomics/harmony",'force = TRUE')
# devtools::install_github("immunogenomics/lisi")
# remotes::install_github("satijalab/seurat")
# devtools::install_github('MacoskoLab/liger')
# install.packages("rliger")
# devtools::install_github('satijalab/seurat-data')
 remotes::install_github('satijalab/seurat-wrappers', force = "TRUE")
# BiocManager::install("scater")
# BiocManager::install("scMerge")
# install.packages("kernlab")
# install.packages("NMI")
# BiocManager::install("CellMixS")

library(ggplot2)
library(SingleCellExperiment)
library(scMerge)
library(scater)
library(patchwork)
library(mclust)
library(lisi)
library(cluster)
library(NMI)
library(harmony)
library(Matrix)
library(sctransform)
library(BiocSingular)
library(rliger)
library(Seurat)
library(liger)
library(SeuratData)
library(SeuratWrappers)
library(CellMixS)

setwd("C:/Users/edbels/Documents/GitHub/mastat_thesis")

source("seurat_metrics.R")
source("DI_function.R")
source("calculate_mmd.R")
# source("seuratLiger.R")


# I have to increase the memory limit for the program to work
memory.limit(size=6e6)

tot_iterations <- 1


testset <- "cbmc"
asu <- "RNA"

so <- readRDS(paste0("datasets/output/",testset,".rds"))

if(file.exists(paste0("results/",testset,"/mmd.rds"))){
  mmd <- readRDS(paste0("results/",testset,"/mmd.rds"))
} else {
  mmd <- cal_mmd(so)
  saveRDS(mmd,file = paste0("results/",testset,"/mmd.rds"))
}


m_mmd <- apply(data.frame(mmd),1, function(x) min(x[x>0]))
# 0 -> split in groups of equal size
# 1- 999 (integer values) -> split in equal sizes, but cell type X only occurs in 1 batch
# 0.001 - 0.999 split in unequal groups, with x*100% in one group
# 1.001 - 1.999 split in equal groups, but in 1 group, the counts are subsampled by (x-1)*100%
# 
testvalues <- c(1.25,1.5, 1.75, which.min(m_mmd),
                 which.min(abs(m_mmd - quantile(m_mmd,0.25))),
                 which.min(abs(m_mmd - median(m_mmd))),
                 which.min(abs(m_mmd - quantile(m_mmd,0.75))),
                 which.max(m_mmd)
)
#
#testvalues <- c(0.5,0.1,0.25, 1.1)

nPCA <- 25

for (split_batch in testvalues){
  print(split_batch)

  nx1<- matrix(0,nrow = tot_iterations, ncol = 1)
  nx2 <- matrix(0,nrow = tot_iterations, ncol = 2)
  nx3 <- matrix(0,nrow = tot_iterations, ncol = 3)
  hvg <- list(nx1,nx1,nx1,nx1,nx1,nx1)
  kbet <- list(nx3,nx3,nx3,nx3,nx3,nx3)
  lisi <- list(nx3,nx3,nx3,nx3,nx3,nx3)
  asw <- list(nx2,nx2,nx2,nx2,nx2,nx2)
  nmi<- list(nx2,nx2,nx2,nx2,nx2,nx2)
  ari <- list(nx2,nx2,nx2,nx2,nx2,nx2)
  cms <- list(1,1,1,1,1,1)
  rm("nx1","nx2","nx3")
  
  
  filename <-  paste0("results/",testset,"/",split_batch,"_it",tot_iterations,".RData")
  
  
  if (file.exists(filename)) {
    load(filename) 
  } else {
    i_START <- 1
  }
  
  for (i in i_START:tot_iterations){
    so <- readRDS(paste0("datasets/output/",testset,".rds"))
    ######### split the dataset ########### 
    gt <- Idents(so)#so@meta.data$cell_type
    if (split_batch >1 && split_batch<2) {
      groups <- sample(c("1", "2"), 
                       size = dim(GetAssayData(object = so, slot = "counts"))[2],
                       replace = TRUE)
      names(groups) <- colnames(so)
      so <- AddMetaData(object = so, metadata = groups, col.name = "group")
      so <- AddMetaData(object = so, metadata = gt, col.name = "cell_type")
      
      so.list <- SplitObject(so, split.by = "group")
      
      md <- so.list[[2]]@meta.data
      
      group2 <- as.sparse(GetAssayData(object = so.list[[2]], slot = "counts"))
      
      group2@x <- as.numeric(sapply(group2@x,rbinom,n= 1,prob = (split_batch-1)))
      
      so2 <- CreateSeuratObject(counts = group2)
      so2 <- AddMetaData(object = so2, metadata = md$group, col.name = "group")
      so2 <- AddMetaData(object = so2, metadata = md$cell_type, col.name = "cell_type")
      
      so <- merge(so.list[[1]],so2)
      
      #so <- NormalizeData(so)
      #so <- FindVariableFeatures(so)
      #so <- ScaleData(so)
      #so <- RunPCA(so, verbose = FALSE)
      rm('so.list',"so2","group2","md","groups")
      
    } else {
      if (split_batch ==0) {
        groups <- sample(c("1", "2"), 
                         size = dim(GetAssayData(object = so, slot = "counts"))[2],
                         replace = TRUE)
        # split the object in batches -> add metadata with batch information
        names(groups) <- colnames(so)
        so <- AddMetaData(object = so, metadata = groups, col.name = "group")
        so <- AddMetaData(object = so, metadata = gt, col.name = "cell_type")
        
        rm('groups')
      } else if (split_batch<1) {
        groups <- sample(c("1", "2"),
                         size = dim(GetAssayData(object = so, slot = "counts"))[2],
                         replace = TRUE, 
                         prob = c(split_batch,1-split_batch))
        # split the object in batches -> add metadata with batch information
        names(groups) <- colnames(so)
        so <- AddMetaData(object = so, metadata = groups, col.name = "group")
        so <- AddMetaData(object = so, metadata = gt, col.name = "cell_type")
        
        rm('groups')
      } else {
        groups <- sample(c("1", "2"),
                         size = dim(GetAssayData(object = so, slot = "counts"))[2],
                         replace = TRUE)
        groups[so@meta.data$cell_type== levels(so)[split_batch]] <- 2
        
        so <- AddMetaData(object = so, metadata = groups, col.name = "group")
        so <- AddMetaData(object = so, metadata = gt, col.name = "cell_type")
        
        so.list <- SplitObject(so, split.by = "group")
        
        md <- so.list[[2]]@meta.data
        
        group2 <- as.sparse(GetAssayData(object = so.list[[2]], slot = "counts"))
        
        group2@x <- as.numeric(sapply(group2@x,rbinom,n= 1,prob = 0.5))
        
        so2 <- CreateSeuratObject(counts = group2)
        so2 <- AddMetaData(object = so2, metadata = md$group, col.name = "group")
        so2 <- AddMetaData(object = so2, metadata = md$cell_type, col.name = "cell_type")
        
        so <- merge(so.list[[1]],so2)
        rm('groups')
      }
    }
    rm('gt')
    print(so@meta.data$group[1:50])
    so@meta.data <- so@meta.data[, c("group","cell_type")]
    
    # print("liger")
    # LIGERl <- LIGER(so, batch = "group", nPCAs = nPCA)
    # so_LIGER <- LIGERl[[1]]
    # LIGERl[[1]] <- 0
    # perf<- seurat.performance(so_LIGER,"cell_type","group","seurat_clusters")
    # hvg[[6]][i]<- as.numeric(perf[1])
    # kbet[[6]][i,] <- as.vector(unlist(perf[2]))
    # lisi[[6]][i,] <- as.vector(unlist(perf[3]))
    # asw[[6]][i,] <- as.vector(unlist(perf[4]))
    # nmi[[6]][i,] <- as.vector(unlist(perf[5]))
    # ari[[6]][i,] <- as.vector(unlist(perf[6]))
    # cms[6] <- perf[7]
    # rm("so_LIGER","perf")
    # 


    print("merge")
    merge <- smerge(so, batch = "group", nPCAs = nPCA)
    so_merge <- merge[[1]] 
      merge[[1]] <- 0
    perf <- seurat.performance(so_merge,"cell_type","group","seurat_clusters")
    hvg[[1]][i]<- as.numeric(perf[1])
    kbet[[1]][i,] <- as.vector(unlist(perf[2]))
    lisi[[1]][i,] <- as.vector(unlist(perf[3]))
    asw[[1]][i,] <- as.vector(unlist(perf[4]))
    nmi[[1]][i,] <- as.vector(unlist(perf[5]))
    ari[[1]][i,] <- as.vector(unlist(perf[6]))
    cms[1] <- perf[7]
    rm("so_merge","perf")
    
    print("scMerge")
    scm <- scMerge_so(so, batch = "group", nPCAs = nPCA)
    so_scm <- scm[[1]] 
    scm[[1]] <- 0
    perf<- seurat.performance(so_scm,"cell_type","group","seurat_clusters")
    hvg[[5]][i]<- as.numeric(perf[1])
    kbet[[5]][i,] <- as.vector(unlist(perf[2]))
    lisi[[5]][i,] <- as.vector(unlist(perf[3]))
    asw[[5]][i,] <- as.vector(unlist(perf[4]))
    nmi[[5]][i,] <- as.vector(unlist(perf[5]))
    ari[[5]][i,] <- as.vector(unlist(perf[6]))
    cms[5] <- perf[7]
    rm("so_scm","perf")
    
    print("seurat")
    seuratl <- seurat(so, batch = "group", nPCAs = nPCA)
    so_seurat <- seuratl[[1]] 
      seuratl[[1]] <- 0
    perf<- seurat.performance(so_seurat,"cell_type","group","seurat_clusters")
    hvg[[2]][i]<- as.numeric(perf[1])
    kbet[[2]][i,] <- as.vector(unlist(perf[2]))
    lisi[[2]][i,] <- as.vector(unlist(perf[3]))
    asw[[2]][i,] <- as.vector(unlist(perf[4]))
    nmi[[2]][i,] <- as.vector(unlist(perf[5]))
    ari[[2]][i,] <- as.vector(unlist(perf[6]))
    cms[2] <- perf[7]
    rm("so_seurat","perf")
    
    print("harmony")
    harm <- harmony(so, batch = "group", nPCAs = nPCA, au = asu)
    so_harmony <- harm[[1]] 
      harm[[1]] <- 0
    perf<- seurat.performance(so_harmony,"cell_type","group","seurat_clusters")
    hvg[[3]][i]<- as.numeric(perf[1])
    kbet[[3]][i,] <- as.vector(unlist(perf[2]))
    lisi[[3]][i,] <- as.vector(unlist(perf[3]))
    asw[[3]][i,] <- as.vector(unlist(perf[4]))
    nmi[[3]][i,] <- as.vector(unlist(perf[5]))
    ari[[3]][i,] <- as.vector(unlist(perf[6]))
    cms[3] <- perf[7]
    rm("so_harmony","perf")
    
    print("sct")
    sct <- scTransform(so, batch = "group", nPCAs = nPCA)
    so_sct <- sct[[1]]
      sct[[1]] <- 0
    perf<- seurat.performance(so_sct,"cell_type","group","seurat_clusters")
    hvg[[4]][i]<- as.numeric(perf[1])
    kbet[[4]][i,] <- as.vector(unlist(perf[2]))
    lisi[[4]][i,] <- as.vector(unlist(perf[3]))
    asw[[4]][i,] <- as.vector(unlist(perf[4]))
    nmi[[4]][i,] <- as.vector(unlist(perf[5]))
    ari[[4]][i,] <- as.vector(unlist(perf[6]))
    cms[4] <- perf[7]
    rm("so_sct","perf")
    
    i_START <- i+1
    save(hvg, kbet, lisi, nmi, ari, asw, cms, i_START, file = filename)
    if (tot_iterations == 1){
      save(scm, merge, seuratl, harm, sct, file = paste0("results/",testset,"/",split_batch,"_it",tot_iterations,"_umap.RData"))
      rm(scm, merge, seuratl, harm, sct)
      #, LIGERl
    }
  
  }
}
