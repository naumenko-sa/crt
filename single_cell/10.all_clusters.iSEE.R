# Rscript 10.all_clusters.iSEE.R
library(iSEE)
library(Seurat)
# Viewing clusters in iSEE
if (file.exists("data/sse.clusters.RDS")){
    sse <- readRDS("data/sse.clusters.RDS")
}else{
    seurat <- readRDS("data/seurat.clusters.RDS")
    sse <- as.SingleCellExperiment(seurat)
    saveRDS(sse, "data/sse.clusters.RDS")
}
iSEE(sse)
