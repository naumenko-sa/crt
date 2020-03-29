# Rscript 10.all_clusters.iSEE.R
library(iSEE)
library(Seurat)
library(org.Mm.eg.db)
library(singleCellTK)
# Viewing clusters in iSEE
if (file.exists("data/pfc/sce.clusters.RDS")){
    sce <- readRDS("data/pfc/sce.clusters.RDS")
}else{
    seurat <- readRDS("data/pfc/seurat.clusters.RDS")
    sce <- as.SingleCellExperiment(seurat)
    sce <- convertGeneIDs(inSCE = sce,
                          inSymbol = "ENSEMBL",
                          outSymbol = "SYMBOL",
                          database = "org.Mm.eg.db")
    saveRDS(sce, "data/pfc/sce.clusters.RDS")
}
iSEE(sce)
