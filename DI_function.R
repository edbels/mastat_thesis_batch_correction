smerge <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  #clustering -> creates metadata suerat_clusters
  seurat_object <- FindNeighbors(seurat_object, dims = 1:nPCAs)
  seurat_object <- FindClusters(seurat_object, resolution = 0.8)
  
  # UMAP and TSNE
  seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:nPCAs)
  seurat_object <- RunTSNE(seurat_object, perplexity=30, do.fast = T)
  
  # Run the standard workflow for visualization
  # for the different batches
  p1merge <- DimPlot(seurat_object, reduction = reduc, group.by = batch)
  p2merge <- DimPlot(seurat_object, reduction = reduc, group.by = "seurat_clusters")
  p1merge + p2merge 
  
  return(seurat_object)
}

seurat <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  seurat_object.list <- SplitObject(seurat_object, split.by = batch)
  
  # normalize data and HVG
  #for (i in 1:length(seurat_object.list)) {
  #  seurat_object.list[[i]] <- NormalizeData(seurat_object.list[[i]], verbose = FALSE)
  #  seurat_object.list[[i]] <- FindVariableFeatures(seurat_object.list[[i]], selection.method = "vst", 
  #                                             nfeatures = 2000, verbose = FALSE)
  #}
  
  # merge
  reference.list <- seurat_object.list
  seurat_object.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:nPCAs)
  seurat_object.integrated <- IntegrateData(anchorset = seurat_object.anchors, dims = 1:nPCAs)
  rm('cmbc.list','reference.list')
  
  # scale, find HVG and run PCA
  seurat_object.integrated <- ScaleData(seurat_object.integrated, verbose = FALSE)
  seurat_object.integrated <- FindVariableFeatures(seurat_object.integrated)
  seurat_object.integrated <- RunPCA(seurat_object.integrated, npcs = 30, verbose = FALSE)
  
  #clustering -> creates metadata suerat_clusters
  seurat_object.integrated <- FindNeighbors(seurat_object.integrated, dims = 1:nPCAs)
  seurat_object.integrated <- FindClusters(seurat_object.integrated, resolution = 0.8)
  
  #UMAP and TSNE
  seurat_object.integrated <- RunUMAP(seurat_object.integrated, reduction = "pca", dims = 1:nPCAs)
  #seurat_object.integrated <- RunTSNE(seurat_object.integrated, perplexity=30, do.fast = T)
  
  # Run the standard workflow for visualization
  # for the different batches
  p1seurat <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = batch)
  p2seurat <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = "seurat_clusters")
  p1seurat + p2seurat
  
  return(seurat_object.integrated)
}

harmony <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  # merge
  seurat_object.integrated <- RunHarmony(seurat_object,batch)
  
  #clustering -> creates metadata suerat_clusters
  seurat_object.integrated <- FindNeighbors(seurat_object.integrated, reduction = "harmony",dims = 1:nPCAs)
  seurat_object.integrated <- FindClusters(seurat_object.integrated, resolution = 0.8)
  
  # tsne and umap for visualisation
  seurat_object.integrated <- RunUMAP(seurat_object.integrated, reduction = "harmony", dims = 1:nPCAs)
  #seurat_object.integrated <- RunTSNE(seurat_object.integrated, reduction = "harmony", perplexity=30, do.fast = T)
  
  # Run the standard workflow for visualization
  # for the different batches
  p1harmony <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = batch)
  p2harmony <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = "seurat_clusters")
  
  p1harmony + p2harmony
  return(seurat_object.integrated)
}

scTransform <- function(seurat_object, batch, nPCAs = 20, reduc = "umap") {
  # merge
  seurat_object.integrated <- SCTransform(seurat_object, vars.to.regress = batch, verbose = FALSE)
  
  # run PCA (scale, and HVG imbedded in SCTransform)
  seurat_object.integrated <- RunPCA(seurat_object.integrated, npcs = nPCAs, verbose = FALSE)
  
  #clustering -> creates metadata suerat_clusters
  seurat_object.integrated <- FindNeighbors(seurat_object.integrated, reduction = "pca",dims = 1:nPCAs)
  seurat_object.integrated <- FindClusters(seurat_object.integrated, resolution = 0.8)
  
  # UMAP and tsne
  seurat_object.integrated <- RunUMAP(seurat_object.integrated, reduction = "pca", dims = 1:nPCAs)
  #seurat_object.integrated <- RunTSNE(seurat_object.integrated, reduction = "pca", perplexity=30, do.fast = T)
  
  # Run the standard workflow for visualization
  # for the different batches
  p1sct <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = batch)
  p2sct <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = "seurat_clusters")
  p1sct + p2sct 
  
  
  return(seurat_object.integrated)
}

scMerge_so <- function(seurat_object, batch, nPCAs = 20, reduc = "umap") {
  # merge
  # make single cell experiment
  seurat_object.sce <- as.SingleCellExperiment(seurat_object)
  # select negative controls: stably expressed genes
  exprs_mat = SummarizedExperiment::assay(seurat_object.sce, 'counts')
  result = scSEGIndex(exprs_mat = exprs_mat)
  # number of clusters in datasets (here equal, change for other datasets)
  K <- length(levels(seurat_object.sce@colData@listData$cell_type))
  # remove other datasets before doing the SCMerge, since it is memory demanding
  rm("seurat_object")
  # do the merging
  seurat_object.sce <- scMerge(
    sce_combine = seurat_object.sce, 
    ctl = rownames(result),
    kmeansK = c(K,K),
    assay_name = "scMerge_unsupervised",
    batch_name ="group",
    BSPARAM = IrlbaParam(), 
    svd_k = 20)
  # remake seurat object
  seurat_object.integrated <- as.Seurat(seurat_object.sce, counts = "scMerge_unsupervised")
  rm("seurat_object.sce")
  
  # run PCA
  seurat_object.integrated <- ScaleData(seurat_object.integrated, verbose = FALSE)
  seurat_object.integrated <- FindVariableFeatures(seurat_object.integrated)
  seurat_object.integrated <- RunPCA(seurat_object.integrated, npcs = nPCAs, verbose = FALSE)
  
  #clustering -> creates metadata suerat_clusters
  seurat_object.integrated <- FindNeighbors(seurat_object.integrated, reduction = "pca",dims = 1:nPCAs)
  seurat_object.integrated <- FindClusters(seurat_object.integrated, resolution = 0.8)
  
  
  # umap and tsne
  seurat_object.integrated <- RunUMAP(seurat_object.integrated, reduction = "pca", dims = 1:nPCAs)
  #seurat_object.integrated <- RunTSNE(seurat_object.integrated, reduction = "pca", perplexity=30, do.fast = T)
  
  # Run the standard workflow for visualization
  # for the different batches
  p1scm <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = batch)
  p2scm <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = "seurat_clusters")
  p1scm + p2scm 
  
  
  return(seurat_object.integrated)
}