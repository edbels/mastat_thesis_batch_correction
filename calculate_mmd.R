library(kernlab)

cal_mmd <- function(seurat_object){

  seurat_object.list <- SplitObject(seurat_object, split.by = "cell_type")
  
  output <- matrix(, nrow = length(seurat_object.list), ncol = length(seurat_object.list))
  
  for (i in c(1:length(seurat_object.list))){
    for (j in c(i:length(seurat_object.list))){
      if (i == j){
        output[i,i] <- 0
      } else { 
        x <- t(as.matrix(GetAssayData(object = seurat_object.list[[i]], slot = "counts")))
        y <- t(as.matrix(GetAssayData(object = seurat_object.list[[j]], slot = "counts")))
        
        mmdstat <- kernlab::kmmd(x,y)
        
        output[i,j] <- mmdstat@mmdstats[1]
        output[j,i] <- mmdstat@mmdstats[1]
      }
    }
  }
  
  return (output)
}
