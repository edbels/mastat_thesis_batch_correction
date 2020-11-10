#### creating a second data-set for the batch effect correction ####
# based on the following vignettes: 
# https://satijalab.org/seurat/v3.2/hashing_vignette.html
# we will only use the singlets
# processing the data set in the same way but with sctransform
# https://satijalab.org/seurat/v3.2/sctransform_vignette.html
# trying to add biological lables to the clusters using
# https://satijalab.org/seurat/v4.0/reference_mapping.html

remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
remotes::install_github("jlmelville/uwot")
remotes::install_github("mojaveazure/seurat-disk")


library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(sctransform)
library(patchwork)

input.dir <-"~/Projects/10_Batchcor/input"
output.dir <-"~/Projects/10_Batchcor/output"

#### 1. Creating a data set of singlets from HTO data ####

# Load in the UMI matrix
pbmc.umis <- readRDS(file.path(input.dir,"pbmc_umi_mtx.rds"))

# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
pbmc.htos <- readRDS(file.path(input.dir,"pbmc_hto_mtx.rds"))

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using the
# default settings
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(pbmc.hashtag$HTO_classification.global)
Idents(pbmc.hashtag) <- "HTO_classification.global"

# getting out the singlets and creating an object with only the RNA counts
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")
pbmc.assay <- GetAssay(pbmc.singlet,"RNA")
pbmc <- CreateSeuratObject(counts = pbmc.assay@counts)

# removing some other objects
rm(pbmc.assay,pbmc.hashtag,pbmc.htos,pbmc.singlet,pbmc.umis)


#### 2. vignette to process the data with SCtransform to get our groundtruth clustering ####

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE) # default 0.8
DimPlot(pbmc, label = TRUE) + NoLegend()
table(Idents(pbmc)) # this is our ground truth for this data set
pbmc[["orig.clusters"]] <- Idents(pbmc)

#### 3. Try to add biological labels to your data using new reference mapping ####

# loading the reference data set, contains multiple samples and batches, batch corrected with Seurat procedure
reference <- LoadH5Seurat(file.path(input.dir,"pbmc_multimodal.h5seurat"))
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
dim(reference) #20729 genes, 161764 cells

# find anchors for transfer of labels  
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# doing the transfer
pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1", 
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

# explore the mapping results
p1 <- DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 <- DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2

# they have created two clustering labels
orig.clus.l1 <- table(pbmc$predicted.celltype.l1,pbmc$orig.clusters) # a macro level
orig.clus.l2 <- table(pbmc$predicted.celltype.l2,pbmc$orig.clusters) # a more detailed level

color <- colorRampPalette(c("white","green", "blue", "red","yellow"))
hm.l1 <- pheatmap::pheatmap(orig.clus.l1, 
                   display_numbers = orig.clus.l1, number_format = "%.0f",
                   color = color(100))
hm.l2 <- pheatmap::pheatmap(orig.clus.l2, 
                   display_numbers = orig.clus.l2, number_format = "%.0f",
                   color= color(100))

# we can use this to try and annotate our originial clusters
# our original clusters on original umap
p.orig <- DimPlot(pbmc, label = TRUE) + NoLegend()
# color by macro label
p.l1 <- DimPlot(pbmc,  group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p.orig + p.l1
# color by more detailed label
p.l2 <- DimPlot(pbmc,  group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p.orig + p.l2

# based in this we can relabel the cells
new.cluster.ids <- c("T - 1", "NK", "T - 2","B - 1", "CD14 Mono - 1",
                     "T - 3", "CD14 Mono - 2","T- 4","T - 5",
                     "CD16 Mono","CD14 Mono - 3","B - 2","B - 3",
                     "DC","T - 5","Eryth","pDC","Plasmablast")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc[["label.clusters"]] <- Idents(pbmc)

p.l <- DimPlot(pbmc, reduction = "umap", group.by = "label.clusters", 
        label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p.l2 <- DimPlot(pbmc, reduction = "umap", group.by = "predicted.celltype.l1", 
               label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p.l+p.l2

# saving the object 
saveRDS(pbmc,file = file.path(input.dir,"pbmc"))

