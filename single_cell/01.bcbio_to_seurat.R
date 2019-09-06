# cluster: use interact job with 20G RAM
# conda activate r
# which R
# Rscript 00_bcbio_to_seurat.R
# conda deactivate
# Read10X expects to see matrix.mtx.gz barcodes.tsv.gz features.tsv.gz (genes)
library(Seurat)
counts <- Read10X(data.dir = "data", gene.column = 1)
seurat_object <- CreateSeuratObject(counts = counts, min.features = 100)
saveRDS(seurat_object, "seurat.bcbio.RDS")
