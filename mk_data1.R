#### creating a first data-set for the batch effect correction ####
# based on the following vignette: 
# https://satijalab.org/seurat/v3.2/multimodal_vignette.html
# we will only use the RNA part

library(Seurat)
library(ggplot2)
library(patchwork)

input.dir <-"~/Projects/10_Batchcor/input"
output.dir <-"~/Projects/10_Batchcor/output"

# load in the RNA UMI matrix
cbmc.rna <- as.sparse(read.csv(file = file.path(input.dir,"/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz"), sep = ",", 
                               header = TRUE, row.names = 1))
# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# creating the object with the labels
cbmc <- CreateSeuratObject(counts = cbmc.rna)

# standard log-normalization
cbmc <- NormalizeData(cbmc)

# choose ~1k variable features
cbmc <- FindVariableFeatures(cbmc)

# standard scaling (no regression)
cbmc <- ScaleData(cbmc)

# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
cbmc <- RunPCA(cbmc, verbose = FALSE)
ElbowPlot(cbmc, ndims = 50)

cbmc <- FindNeighbors(cbmc, dims = 1:25)
cbmc <- FindClusters(cbmc, resolution = 0.8)
cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")
cbmc <- RunUMAP(cbmc, dims = 1:25)

# looking at the clusters
DimPlot(cbmc, label = TRUE) + NoLegend()
table(Idents(cbmc))
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", 
                     "NK - set 1", "CD14+ Mono - set 1", "Mouse - set 1",
                     "B", "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK - set 2", 
                     "CD34+", "Multiplets", "Mouse - set 2", "Eryth", "Mk", 
                     "Mouse - set 3", "DC", "pDCs")
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)

DimPlot(cbmc, label = TRUE) + NoLegend()

# how many cells and genes do we have?
dim(cbmc)
# 20501 genes and 8617 cells

# saving the object 
saveRDS(cbmc,file = file.path(input.dir,"cbmc"))
