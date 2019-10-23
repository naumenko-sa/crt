# O2: use interact job with 20G RAM
# sometimes Rscript is not working but works from R
# conda activate r
# which R
# Rscript 00_create_seurat_object.R
# conda deactivate
library(Seurat)
counts <- Read10X(data.dir = "data", gene.column = 1)
seurat_object <- CreateSeuratObject(counts = counts, min.features = 100)
saveRDS(seurat_object, "seurat.bcbio.RDS")
