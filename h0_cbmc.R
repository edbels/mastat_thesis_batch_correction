
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


### choose option of split 
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
nPCAs <- 25

### other method to choose number of PCAs?

## because no cell type labels are present: create "ground truth" from clustering
cbmc <- FindNeighbors(cbmc, dims = 1:nPCAs)
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

################## No data integration method #################

# simple merge, so do nothing in this case

#clustering -> creates metadata suerat_clusters
cbmc <- FindNeighbors(cbmc, dims = 1:nPCAs)
cbmc <- FindClusters(cbmc, resolution = 0.8)

# UMAP and TSNE
cbmc <- RunUMAP(cbmc, reduction = "pca", dims = 1:nPCAs)
#cbmc <- RunTSNE(cbmc, perplexity=30, do.fast = T)

# Run the standard workflow for visualization
# for the different batches
reduc <- "umap"
p1merge <- DimPlot(cbmc, reduction = reduc, group.by = "group")
p2merge <- DimPlot(cbmc, reduction = reduc, group.by = "seurat_clusters")
p1merge + p2merge 

perf_merge <- seurat.performance(cbmc,"cell_type","group","seurat_clusters")


######### Seurat v3 ########### 
cbmc.list <- SplitObject(cbmc, split.by = "group")

# normalize data and HVG
#for (i in 1:length(cbmc.list)) {
#  cbmc.list[[i]] <- NormalizeData(cbmc.list[[i]], verbose = FALSE)
#  cbmc.list[[i]] <- FindVariableFeatures(cbmc.list[[i]], selection.method = "vst", 
#                                             nfeatures = 2000, verbose = FALSE)
#}

# merge
reference.list <- cbmc.list
cbmc.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:nPCAs)
cbmc.integrated <- IntegrateData(anchorset = cbmc.anchors, dims = 1:nPCAs)
rm('cmbc.list','reference.list')

# scale, find HVG and run PCA
cbmc.integrated <- ScaleData(cbmc.integrated, verbose = FALSE)
cbmc.integrated <- FindVariableFeatures(cbmc.integrated)
cbmc.integrated <- RunPCA(cbmc.integrated, npcs = 30, verbose = FALSE)

#clustering -> creates metadata suerat_clusters
cbmc.integrated <- FindNeighbors(cbmc.integrated, dims = 1:nPCAs)
cbmc.integrated <- FindClusters(cbmc.integrated, resolution = 0.8)

#UMAP and TSNE
cbmc.integrated <- RunUMAP(cbmc.integrated, reduction = "pca", dims = 1:nPCAs)
#cbmc.integrated <- RunTSNE(cbmc.integrated, perplexity=30, do.fast = T)

# Run the standard workflow for visualization
# for the different batches
reduc <- "umap"
p1seurat <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "group")
p2seurat <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "seurat_clusters")
p1seurat + p2seurat

perf_seurat <- seurat.performance(cbmc.integrated,"cell_type","group","seurat_clusters")
rm("cbmc.integrated")

save.image()

######### Harmony ########### 

# merge
cbmc.integrated <- RunHarmony(cbmc,"group")

#clustering -> creates metadata suerat_clusters
cbmc.integrated <- FindNeighbors(cbmc.integrated, reduction = "harmony",dims = 1:nPCAs)
cbmc.integrated <- FindClusters(cbmc.integrated, resolution = 0.8)

# tsne and umap for visualisation
cbmc.integrated <- RunUMAP(cbmc.integrated, reduction = "harmony", dims = 1:nPCAs)
#cbmc.integrated <- RunTSNE(cbmc.integrated, reduction = "harmony", perplexity=30, do.fast = T)

# Run the standard workflow for visualization
# for the different batches
reduc <- "umap"
p1harmony <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "group")
p2harmony <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "seurat_clusters")

p1harmony + p2harmony 

perf_Harmony <- seurat.performance(cbmc.integrated,"cell_type","group","seurat_clusters")
rm("cbmc.integrated")

######### SCtransform ########### 

# cbmc is actually already normalised here (standard log-normalisation -> change?)

# merge
cbmc.integrated <- SCTransform(cbmc, vars.to.regress = "group", verbose = FALSE)

# run PCA (scale, and HVG imbedded in SCTransform)
cbmc.integrated <- RunPCA(cbmc.integrated, npcs = nPCAs, verbose = FALSE)

#clustering -> creates metadata suerat_clusters
cbmc.integrated <- FindNeighbors(cbmc.integrated, reduction = "pca",dims = 1:nPCAs)
cbmc.integrated <- FindClusters(cbmc.integrated, resolution = 0.8)

# UMAP and tsne
cbmc.integrated <- RunUMAP(cbmc.integrated, reduction = "pca", dims = 1:nPCAs)
#cbmc.integrated <- RunTSNE(cbmc.integrated, reduction = "pca", perplexity=30, do.fast = T)

# Run the standard workflow for visualization
# for the different batches
reduc <- "umap"
p11sct <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "group")
p22sct <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "seurat_clusters")

p1sct + p2sct 

perf1_SCTransform <- seurat.performance(cbmc.integrated,"cell_type","group","seurat_clusters")

rm("cbmc.integrated")


######### SCMerge ########### 

# merge
# make single cell experiment
cbmc.sce <- as.SingleCellExperiment(cbmc)
# select negative controls: stably expressed genes
exprs_mat = SummarizedExperiment::assay(cbmc.sce, 'counts')
result = scSEGIndex(exprs_mat = exprs_mat)
# number of clusters in datasets (here equal, change for other datasets)
K <- length(levels(cbmc.sce@colData@listData$cell_type))
# remove other datasets before doing the SCMerge, since it is memory demanding
rm("cbmc")
# do the merging
cbmc.sce <- scMerge(
  sce_combine = cbmc.sce, 
  ctl = rownames(result),
  kmeansK = c(K,K),
  assay_name = "scMerge_unsupervised",
  batch_name ="group",
  BSPARAM = IrlbaParam(), 
  svd_k = 20)
# remake seurat object
cbmc.integrated <- as.Seurat(scMerge_unsupervised, counts = "scMerge_unsupervised")
rm("cbmc.sce")

# run PCA
cbmc.integrated <- ScaleData(cbmc.integrated, verbose = FALSE)
cbmc.integrated <- FindVariableFeatures(cbmc.integrated)
cbmc.integrated <- RunPCA(cbmc.integrated, npcs = nPCAs, verbose = FALSE)

#clustering -> creates metadata suerat_clusters
cbmc.integrated <- FindNeighbors(cbmc.integrated, reduction = "pca",dims = 1:nPCAs)
cbmc.integrated <- FindClusters(cbmc.integrated, resolution = 0.8)


# umap and tsne
cbmc.integrated <- RunUMAP(cbmc.integrated, reduction = "pca", dims = 1:nPCAs)
#cbmc.integrated <- RunTSNE(cbmc.integrated, reduction = "pca", perplexity=30, do.fast = T)

# Run the standard workflow for visualization
# for the different batches
reduc <- "umap"
p1scm <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "group")
p2scm <- DimPlot(cbmc.integrated, reduction = reduc, group.by = "seurat_clusters")
p1scm + p2scm 

perf_SCMerge <- seurat.performance(cbmc.integrated,"cell_type","group","seurat_clusters")

save.image()

