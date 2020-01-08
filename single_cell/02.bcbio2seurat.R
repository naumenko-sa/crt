# O2: use interactive job with 20G RAM
# sometimes Rscript is not working but works from R
# conda activate r
# which R
# Rscript 00_create_seurat_object.R
# conda deactivate
# should have 3 files from bcbio
library(R.utils)

file.rename("tagcounts.mtx", "matrix.mtx")
file.rename("tagcounts.mtx.rownames", "features.tsv")
file.rename("tagcounts.mtx.colnames", "barcodes.tsv")

gzip("matrix.mtx")
gzip("features.tsv")
gzip("barcodes.tsv")

library(Seurat)
counts <- Read10X(data.dir = ".", gene.column = 1)
seurat_object <- CreateSeuratObject(counts = counts, min.features = 100)
saveRDS(seurat_object, "seurat.bcbio.RDS")
