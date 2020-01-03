# from Lorena Pantano's original template: hbc_kressler_sRNAseq_anxiety_brain_mouse_hbc02185

# run:
# cd to project/final/project
# conda activate r
# Rscript 01.bcbio2sse.R

library(SingleCellExperiment)
library(Matrix)
library(AnnotationHub)
library(tidyverse)

species = "Mus musculus"
counts = readMM(file.path("tagcounts.mtx"))
dupcounts = readMM(file.path("tagcounts-dupes.mtx"))
rownames = read.csv(file.path("tagcounts.mtx.rownames"), header = F)[["V1"]]
rownames = as.character(rownames)
colnames = read.csv(file.path("tagcounts.mtx.colnames"), header = F)[["V1"]]
colnames = make.names(as.character(colnames))
reads = read.csv(file.path("cb-histogram.txt"), header = F, sep="\t", row.names = 1)
rownames(reads) = make.names(rownames(reads))

counts =  as(counts, "dgCMatrix")
rownames(counts) = rownames
colnames(counts) = colnames
metadata = read.csv(file.path("tagcounts.mtx.metadata"))
rownames(metadata) = colnames
metadata[["nUMI"]] = colSums(counts)
metadata[["nGenes"]] = colSums(counts>0)
metadata[["log10GenesPerUMI"]] = log10(metadata$nGene) / log10(metadata$nUMI)
metadata[["nReads"]] = reads[colnames,]
metadata[["saturation_rate"]] = 1-(colSums(counts)/colSums(dupcounts))
metadata[["dupReads"]] = colSums(dupcounts)
metadata[["dupMeanReads"]] = colMeans(dupcounts)

# check if file is empty and skip if the case
# annotation was download from ensembl biomart to match the version GRCh38.92
# AnnotationHub can be used.

## Load the annotation resource.
ah <- AnnotationHub()
ahDb <- query(ah, pattern= c(species, "EnsDb") )
ahEdb <- ahDb[[rev(names(ahDb))[1]]] # last one is chosen
rows = genes(ahEdb) %>%
    as.data.frame() %>%
    janitor::clean_names() %>%
    dplyr::select(gene_id,
                  gene_name,
                  description,
                  biotype = gene_biotype,
                  entrezid,
                  chrom = seqnames) %>%
    group_by(gene_id, biotype, chrom, description) %>%
    summarise(gene_name = paste(unique(gene_name), collapse = ","),
              entrezid = paste(unique(entrezid), collapse = ",")) %>%
    mutate(gene_name=ifelse(gene_name=="", gene_id, gene_name)) %>%
    as.data.frame()

# mit
rrna = rows %>% dplyr::filter(chrom == "MT") %>% .[["gene_id"]] %>% intersect(., rownames)

metadata[["mtUMI"]] = colSums(counts[rrna,], na.rm = T)
metadata[["mtUMI"]][is.na(metadata[["mtUMI"]])] = 0
metadata[["mitoRatio"]] = metadata$mtUMI/metadata$nUMI

se = SingleCellExperiment(assays=list(raw = counts), colData = metadata)
saveRDS(se, "se.RDS")
