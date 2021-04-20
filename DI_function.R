# below are all the functions to perform data integration. 

smerge <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  #clustering -> creates metadata suerat_clusters
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object,split.by = batch)
  seurat_object <- RunPCA(seurat_object, verbose = FALSE)
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
  
  return(list(seurat_object,p1merge,p2merge))
}

seurat <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object,split.by = batch)
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
  seurat_object.integrated <- ScaleData(seurat_object.integrated,split.by = batch, verbose = FALSE)
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
  
  return(list(seurat_object.integrated,p1seurat,p2seurat))
}

harmony <- function(seurat_object, batch, nPCAs = 30, reduc = "umap", au ) {
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object,split.by = batch)
  seurat_object <- RunPCA(seurat_object, verbose = FALSE)
  # merge
  seurat_object.integrated <- RunHarmony(seurat_object,batch, assay.use = au)
  
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
  
  return(list(seurat_object.integrated,p1harmony,p2harmony))
}

scTransform <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  # merge
  seurat_object.integrated <- SCTransform(seurat_object, vars.to.regress = batch, verbose = FALSE)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object,split.by = batch)
  
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
  
  
  return(list(seurat_object.integrated,p1sct,p2sct))
}

scMerge_so <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  seurat_object.list <- SplitObject(seurat_object, split.by = batch)
  cs <- rep(1,length(seurat_object.list))
  # Find number of cells types in each batch
  for (i in 1:length(seurat_object.list)) {
    cs[i] <- nlevels(as.factor(seurat_object.list[[i]]@meta.data$cell_type))
  }
  
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object,split.by = batch)
  # merge
  # make single cell experiment
  seurat_object.sce <- as.SingleCellExperiment(seurat_object)
  seurat_object.sce  = runPCA(seurat_object.sce, exprs_values = "logcounts")
  # select negative controls: stably expressed genes
  exprs_mat = SummarizedExperiment::assay(seurat_object.sce, 'counts')
  result = scSEGIndex(exprs_mat = exprs_mat)
  rm("exprs_mat")
  # number of clusters in datasets (here equal, change for other datasets)
  # remove other datasets before doing the SCMerge, since it is memory demanding
  rm("seurat_object")
  # do the merging
  seurat_object.sce <- scMerge(
    sce_combine = seurat_object.sce, 
    ctl = rownames(result),
    kmeansK = cs,
    assay_name = "scMerge_supervised",
    batch_name = batch,
    cell_type = seurat_object.sce$cell_type)
  # remake seurat object
  seurat_object.integrated <- as.Seurat(seurat_object.sce, counts = "counts")# "scMerge_supervised")
  print(1)
  rm("seurat_object.sce")
  
  # run PCA
  #seurat_object.integrated <- NormalizeData(seurat_object.integrated)
  seurat_object.integrated <- ScaleData(seurat_object.integrated, verbose = FALSE)
  seurat_object.integrated <- FindVariableFeatures(seurat_object.integrated)
  seurat_object.integrated <- RunPCA(seurat_object.integrated, npcs = nPCAs, verbose = FALSE)
  nPCAs = max(dim(seurat_object.integrated@reductions$pca@feature.loadings)[2],nPCAs)
  
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
  
  
  return(list(seurat_object.integrated,p1scm,p2scm))
}

LIGERs <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  seurat_object.list <- SplitObject(seurat_object, split.by = batch)
  data.liger <- seuratToLiger(seurat_object, combined.seurat = TRUE, meta.var = "group")
  data.liger <- normalize(data.liger)
  data.liger <- selectGenes(data.liger, var.thresh = 0.1, do.plot = T)
  data.liger <- scaleNotCenter(data.liger)
  # merge
  data.liger <- optimizeALS(data.liger, k = nPCAs, lamda = 5) 
  data.liger <- quantile_norm(data.liger)
  
  ligerToSeurat(
    data.liger
  )
  
  #clustering -> creates metadata suerat_clusters
  seurat_object.integrated <- FindNeighbors(seurat_object.integrated, reduction = "iNMF",dims = 1:nPCAs)
  seurat_object.integrated <- FindClusters(seurat_object.integrated, resolution = 0.8)
  
  # tsne and umap for visualisation
  seurat_object.integrated <- RunUMAP(seurat_object.integrated, reduction = "iNMF", dims = 1:nPCAs)
  #seurat_object.integrated <- RunTSNE(seurat_object.integrated, reduction = "iNMF", perplexity=30, do.fast = T)
  
  # Run the standard workflow for visualization
  # for the different batches
  p1liger <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = batch)
  p2liger <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = "seurat_clusters")
  
  p1liger + p2liger
  return(list(seurat_object.integrated,p1liger,p2liger))
}

LIGER <- function(seurat_object, batch, nPCAs = 30, reduc = "umap") {
  # merge
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object,split.by = batch,  do.center=FALSE)
  
  
  seurat_object.integrated <- RunOptimizeALS(seurat_object,k = nPCAs, lambda = 5, split.by = batch)
  seurat_object.integrated <- RunQuantileNorm(seurat_object.integrated, split.by = batch)
  
  
  #clustering -> creates metadata suerat_clusters
  seurat_object.integrated <- FindNeighbors(seurat_object.integrated, reduction = "iNMF",dims = 1:nPCAs)
  seurat_object.integrated <- FindClusters(seurat_object.integrated, resolution = 0.8)
  
  # tsne and umap for visualisation
  seurat_object.integrated <- RunUMAP(seurat_object.integrated, reduction = "iNMF", dims = 1:nPCAs)
  #seurat_object.integrated <- RunTSNE(seurat_object.integrated, reduction = "iNMF", perplexity=30, do.fast = T)
  
  # Run the standard workflow for visualization
  # for the different batches
  p1liger <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = batch)
  p2liger <- DimPlot(seurat_object.integrated, reduction = reduc, group.by = "seurat_clusters")
  
  p1liger + p2liger
  return(list(seurat_object.integrated,p1liger,p2liger))
}