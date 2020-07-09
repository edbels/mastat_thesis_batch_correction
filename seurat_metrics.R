################ function description ##########################

#### returns HVG conservation, kBET (k = 10% and 20%), LISI, ASW, NMI and ARI 
#### it takes 4 arguments: - a seurat object
####                       - the name of the column in the metadata that contains the ground truth
####                       - the name of the column in the metadata that contains the batch labels
####                       - the name of the column in the metadata that contains cluster labels 
####                                              obtained after data integration

#### remarks:
#### UMAP required and used in case of kBET, LISI and aSW


seurat.performance <- function(seurat_object,cell_type,batch,clustering){
  ######### percentage of features in the separate lists that is also in the mean list ########
  nHVG <- min( 2000, 
               floor(0.5*seurat_object@assays$RNA@counts@Dim[1]))  # number of HVGs to take. Standard is 2000, but if
                                                                  # if less than 4000 genes are present, we keep 50%
  seurat_object.list <- SplitObject(seurat_object, split.by = batch)
  # Find most variable features in the different subgroups
  for (i in 1:length(seurat_object.list)) {
    seurat_object.list[[i]] <- FindVariableFeatures(seurat_object.list[[i]],
                                                    selection.method = "vst", 
                                                  nfeatures = nHVG,
                                                  verbose = FALSE)
  }
  seurat_object <- FindVariableFeatures(seurat_object,
                                        selection.method = "vst",
                                        nfeatures = nHVG)
  # compare most variable features
  hvg_group = VariableFeatures(seurat_object.list[[1]]) #union of all HVG in the different batches
  for (i in 1:length(seurat_object.list)) {
    hvg_group = union(hvg_group,VariableFeatures(seurat_object.list[[i]]))
  }                #this gives the union of HVGs from the different batches
  
  hvg_perc <- length(intersect(hvg_group,VariableFeatures(seurat_object)))/ # this gives the percentage of 
                                                                                #HVGs conserved
    min(length(hvg_group),length(VariableFeatures(seurat_object)))
  rm("seurat_object.list") # to keep memory usage low
  
  ############# kbet test #############
  kbet <- c(0,0,0) #kBET test done for k= 5%, 10% and 20% (higher than 20% always gives very low rejection rates )
  
  df <- data.frame(Embeddings(seurat_object[["umap"]])) 
  if (dim(df)[1]>1e4) subset_size <- 1e4 else subset_size <- dim(df)[1] # for larger datasets, take a subset (1e4 
                                                                        # cells) to keep calculation time reasonable
  subset_id <- sample.int(n = length(seurat_object@meta.data[,batch]),
                          size = subset_size, replace=FALSE)
  df <- df[subset_id,]
  batch.estimate <- kBET::kBET(df, batch = seurat_object@meta.data[subset_id,batch],
                               k0 = round(0.05*dim(df)[1]))
  kbet[1] <- mean(batch.estimate$stats$kBET.observed )
  batch.estimate <- kBET::kBET(df, batch = seurat_object@meta.data[subset_id,batch],
                               k0 = round(0.1*dim(df)[1]))
  kbet[2] <- mean(batch.estimate$stats$kBET.observed )
  batch.estimate <- kBET::kBET(df, batch = seurat_object@meta.data[subset_id,batch],
                               k0 = round(0.2*dim(df)[1]))
  kbet[3] <- mean(batch.estimate$stats$kBET.observed )
  
  #### LISI test ##########
  df <- data.frame(Embeddings(seurat_object[["umap"]]))
  md <- data.frame(seurat_object@meta.data[,c(batch,clustering)])
  lisi_res <- lisi::compute_lisi(df, md, colnames(md),perplexity = 30)
  head(lisi_res)
  boxplot(as.matrix(lisi_res))
  mean(as.matrix(lisi_res))
  ilisi <- median(as.matrix(lisi_res[,1]))
  clisi <- median(as.matrix(lisi_res[,2]))
  # calculate the expected iLISI
  probabilities <- table(seurat_object@meta.data$group)/length(seurat_object@meta.data$group)
  nlisi <- 1/ sum(probabilities^2) # expected value of iLISI with perfect data integration
  lisi <- c(ilisi, clisi, nlisi)
  
  ##### asw test ######
  svalues <- silhouette(as.numeric(seurat_object@meta.data[,batch]),
                        dist = dist(df,method = "euclidean"))
  asw_batch <- mean(svalues[,3])
  svalues <- silhouette(as.numeric(as.factor(seurat_object@meta.data[,clustering])),
                        dist = dist(df,method="euclidean"))
  asw_cluster <- mean(svalues[,3])
  asw <- c(asw_batch, asw_cluster)
  
  ##### NMI test ######
  cell_labels <- data.frame(names(Idents(seurat_object)),
                            seurat_object@meta.data[,cell_type])
  cluster_labels <- data.frame(names(Idents(seurat_object)),
                               seurat_object@meta.data[,clustering])
  batch_labels <- data.frame(names(Idents(seurat_object)),
                             seurat_object@meta.data[,batch])
  nmi_batch <- NMI(batch_labels, cluster_labels)
  nmi_cluster <- NMI(cell_labels, cluster_labels)
  nmi <- c(nmi_batch, nmi_cluster)
  
  #### ARI #######
  ARI_batch <- adjustedRandIndex(seurat_object@meta.data[,batch],
                                 seurat_object@meta.data[,clustering])
  ARI_cluster <- adjustedRandIndex(seurat_object@meta.data[,cell_type],
                                   seurat_object@meta.data[,clustering])
  ari <- c(ARI_batch, ARI_cluster)
  
  return_list <- list("hvg_perc" = hvg_perc, "kbet"=kbet, "lisi" = lisi, "asw" = asw, "NMI" = nmi, "ARI" = ari)
  return (return_list)
  
}
