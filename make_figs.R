# make boxplots over iterations of performance metrics

rm(list = ls())


### choose option of split 
# -> 1: equal split
# -> 2: unequal split
# -> 3: equal split, but with all "Naive CD4 T" in batch one
setwd("C:/Users/edbels/Documents/GitHub/mastat_thesis")

testset <- "Pbmc"
tot_iterations <- 1

i_START <- 1

mmd <- readRDS(paste0("results/",testset,"/mmd.rds"))
m_mmd <- apply(data.frame(mmd),1, function(x) min(x[x>0]))
testvalues <- c(0, 0.1,0.25,0.5, 1.1,1.25,1.5, which.min(m_mmd),
                which.min(abs(m_mmd - quantile(m_mmd,0.25))),
                which.min(abs(m_mmd - median(m_mmd))),
                which.min(abs(m_mmd - quantile(m_mmd,0.75))),
                which.max(m_mmd))

split_batch <- 0.5

fn <- paste0("results/",testset,"/",split_batch,"_it",tot_iterations,"_umap.RData")

load(fn)


# scm
  jpeg(paste0("results/",testset,"/scm_",split_batch,".jpeg"))
  scm[[2]]+scm[[3]]
  dev.off()
#merge  
  jpeg(paste0("results/",testset,"/merge_",split_batch,".jpeg"))
  merge[[2]]+merge[[3]]
  dev.off()
#seuratl  
  jpeg(paste0("results/",testset,"/seurat_",split_batch,".jpeg"))
  seuratl[[2]]+seuratl[[3]]
  dev.off()  
#harm 
  jpeg(paste0("results/",testset,"/harm_",split_batch,".jpeg"))
  harm[[2]]+harm[[3]]
  dev.off()
#sct  
  jpeg(paste0("results/",testset,"/sct_",split_batch,".jpeg"))
  sct[[2]]+sct[[3]]
  dev.off()
#LIGERl 
  jpeg(paste0("results/",testset,"/LIGER_",split_batch,".jpeg"))
  LIGERl[[2]]+LIGERl[[3]]
  dev.off()


  