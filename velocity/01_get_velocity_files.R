library(BUSpaRse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(AnnotationHub)

ah <- AnnotationHub()
query(ah, pattern = c("Ensembl", "97", "Mus musculus", "EnsDb"))
# Get mouse Ensembl 97 annotation
edb <- ah[["AH73905"]]
# note L = read lengths for transcript read, L=61 for indrops3/61
get_velocity_files(edb, L=61, Genome=BSgenome.Mmusculus.UCSC.mm10, out_path = "./veloindex", isoform_action = "separate")

