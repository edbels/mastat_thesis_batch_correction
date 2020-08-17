
rm(list = ls())

library(devtools)
#devtools::install_github('theislab/kBET',force = TRUE)
#devtools::install_github("immunogenomics/harmony", force = TRUE)
#devtools::install_github("immunogenomics/lisi")
#BiocManager::install("scMerge")
library(ggplot2)
library(SingleCellExperiment)
library(scMerge)
library(scater)
library(Seurat)
library(patchwork)
library(mclust)
library(lisi)
library(cluster)
library(NMI)
library(harmony)
library(Matrix)
library(sctransform)
library(BiocSingular)

setwd("C:/Users/edbels/Documents/GitHub/mastat_thesis")

source("seurat_metrics.R")
source("DI_function.R")

tot_iterations <- 1e2

conflicts(detail=TRUE)### choose option of split 
# -> 1: equal split
# -> 2: unequal split
# -> 3: equal split, but with all "Naive CD4 T" in batch one
split_batch <- 1

### load in the file and make a seurat object
cbmc.rna <- as.sparse(read.csv(file = "dataset_seurat_h0/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna) # deletes HUMAN-
cbmc <- CreateSeuratObject(counts = cbmc.rna)
rm("cbmc.rna") #regularly delete variables that are not used anymore, to keep memory usage low

# Normalize and find HVG, default number is 2000
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)

## select PCA's
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
ElbowPlot(cbmc, ndims = 50) # 20 to 30 PCA's seems good
nPCA <- 25

### other method to choose number of PCAs?

## because no cell type labels are present: create "ground truth" from clustering
cbmc <- FindNeighbors(cbmc, dims = 1:nPCA)
cbmc <- FindClusters(cbmc, resolution = 0.8) 
cbmc.rna.markers <- FindAllMarkers(cbmc, 
                                   max.cells.per.ident = 100,
                                   min.diff.pct = 0.3,
                                   only.pos = TRUE)
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", 
                     "CD14+ Mono", "Mouse", "B","CD8 T", "CD16+ Mono", 
                     "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse",
                     "Eryth", "Mk","Mouse", "DC", "pDCs")
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)
gt <- Idents(cbmc)
rm("cbmc.rna.markers")

### explore the data + QC
cbmc[["percent.mt"]] <- PercentageFeatureSet(cbmc, pattern = "^MT-")
head(cbmc@meta.data, 10) # Show QC metrics for the first 5 cells --> no mitochondrial data?
sum(cbmc@meta.data$percent.mt==0) # all percent.mt == 0 --> or naming genes is different?
VlnPlot(cbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#cbmc <- subset(cbmc, subset = nCount_RNA <25000 &
#                 nFeature_RNA < 4500 &
#                 percent.mt < 15) # --> keep all datapoints
rm('plot1','plot2')

nx1<- matrix(0,nrow = tot_iterations, ncol = 1)
nx2 <- matrix(0,nrow = tot_iterations, ncol = 2)
nx3 <- matrix(0,nrow = tot_iterations, ncol = 3)
hvg <- list(nx1,nx1,nx1,nx1,nx1)
kbet <- list(nx3,nx3,nx3,nx3,nx3)
lisi <- list(nx3,nx3,nx3,nx3,nx3)
asw <- list(nx2,nx2,nx2,nx2,nx2)
nmi<- list(nx2,nx2,nx2,nx2,nx2)
ari <- list(nx2,nx2,nx2,nx2,nx2)
rm("nx1","nx2","nx3")

load("cbmc_stability_1.RData")

for (i in 1:tot_iterations){
  ######### split the dataset ########### 
  if (split_batch ==1) {
    groups <- sample(c("1", "2"), 
                     size = dim(GetAssayData(object = cbmc, slot = "counts"))[2],
                     replace = TRUE)
  } else if (split_batch==2) {
    groups <- sample(c("1", "2"),
                     size = dim(GetAssayData(object = cbmc, slot = "counts"))[2],
                     replace = TRUE, 
                     prob = c(0.2,0.8))
  } else if (split_batch ==3) {
    groups <- sample(c("1", "2"),
                     size = dim(GetAssayData(object = cbmc, slot = "counts"))[2],
                     replace = TRUE)
    groups[gt=="Naive CD4 T"] <- 1
  }
  # split the object in batches -> add metadata with batch information
  names(groups) <- colnames(cbmc)
  cbmc <- AddMetaData(object = cbmc, metadata = groups, col.name = "group")
  cbmc <- AddMetaData(object = cbmc, metadata = gt, col.name = "cell_type")
  
  cbmc_merge <- smerge(cbmc, batch = "group", nPCAs = nPCA)
  perf <- seurat.performance(cbmc_merge,"cell_type","group","seurat_clusters")
  hvg[[1]][i]<- as.numeric(perf[1])
  kbet[[1]][i,] <- as.vector(unlist(perf[2]))
  lisi[[1]][i,] <- as.vector(unlist(perf[3]))
  asw[[1]][i,] <- as.vector(unlist(perf[4]))
  nmi[[1]][i,] <- as.vector(unlist(perf[5]))
  ari[[1]][i,] <- as.vector(unlist(perf[6]))
  rm("cbmc_merge","perf")
  
  cbmc_seurat <- seurat(cbmc, batch = "group", nPCAs = nPCA)
  perf<- seurat.performance(cbmc_seurat,"cell_type","group","seurat_clusters")
  hvg[[2]][i]<- as.numeric(perf[1])
  kbet[[2]][i,] <- as.vector(unlist(perf[2]))
  lisi[[2]][i,] <- as.vector(unlist(perf[3]))
  asw[[2]][i,] <- as.vector(unlist(perf[4]))
  nmi[[2]][i,] <- as.vector(unlist(perf[5]))
  ari[[2]][i,] <- as.vector(unlist(perf[6]))
  rm("cbmc_seurat","perf")
  
  cbmc_harmony <- harmony(cbmc, batch = "group", nPCAs = nPCA)
  perf<- seurat.performance(cbmc_harmony,"cell_type","group","seurat_clusters")
  hvg[[3]][i]<- as.numeric(perf[1])
  kbet[[3]][i,] <- as.vector(unlist(perf[2]))
  lisi[[3]][i,] <- as.vector(unlist(perf[3]))
  asw[[3]][i,] <- as.vector(unlist(perf[4]))
  nmi[[3]][i,] <- as.vector(unlist(perf[5]))
  ari[[3]][i,] <- as.vector(unlist(perf[6]))
  rm("cbmc_harmony","perf")
  
  cbmc_sct <- scTransform(cbmc, batch = "group", nPCAs = nPCA)
  perf<- seurat.performance(cbmc_sct,"cell_type","group","seurat_clusters")
  hvg[[4]][i]<- as.numeric(perf[1])
  kbet[[4]][i,] <- as.vector(unlist(perf[2]))
  lisi[[4]][i,] <- as.vector(unlist(perf[3]))
  asw[[4]][i,] <- as.vector(unlist(perf[4]))
  nmi[[4]][i,] <- as.vector(unlist(perf[5]))
  ari[[4]][i,] <- as.vector(unlist(perf[6]))
  rm("cbmc_sct","perf")
  
  cbmc_scm <- scMerge_so(cbmc, batch = "group", nPCAs = nPCA)
  perf<- seurat.performance(cbmc_scm,"cell_type","group","seurat_clusters")
  hvg[[5]][i]<- as.numeric(perf[1])
  kbet[[5]][i,] <- as.vector(unlist(perf[2]))
  lisi[[5]][i,] <- as.vector(unlist(perf[3]))
  asw[[5]][i,] <- as.vector(unlist(perf[4]))
  nmi[[5]][i,] <- as.vector(unlist(perf[5]))
  ari[[5]][i,] <- as.vector(unlist(perf[6]))
  rm("cbmc_scm","perf")
  
  
  save(hvg, kbet, lisi, nmi, ari, asw, file = "cbmc_stability_1.RData")
  
}


save(hvg, kbet, lisi, nmi, ari, asw, file = "cbmc_stability_1.RData")

