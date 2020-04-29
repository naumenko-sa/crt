# create SummarizedExperiment object (se) from bcbio output
# usage: 
# Rscript bcbio2se.R /bcbio/result/final/project
# output: data/bcbio.se.RDS

library(tidyverse)
library(tximport)
library(SummarizedExperiment)
library(janitor)
library(DESeq2)

args <- commandArgs(trailingOnly = T)
project_dir <- args[1]

metadata <- read_csv(file.path(project_dir, "metadata.csv"))
colnames(metadata)[1] <- "sample"

# metadata$genotype <- as.factor(metadata$genotype)

metrics <- read_tsv(file.path(project_dir, 
			      "multiqc", 
			      "multiqc_data", 
			      "multiqc_bcbio_metrics.txt")) %>% 
	   clean_names(case = "snake")

sample_dirs <- file.path(project_dir, "..", metadata$sample)
salmon_files <- file.path(sample_dirs, "salmon", "quant.sf")
names(salmon_files) <- metadata$sample

transcripts2genes_file <- file.path(project_dir, "tx2gene.csv")
transcripts2genes <- read_csv(transcripts2genes_file, col_names = c("ensembl_transcript_id", "ensembl_gene_id"))

txi_salmon <- tximport(salmon_files, type = "salmon", 
                       tx2gene = transcripts2genes,
                       countsFromAbundance = "lengthScaledTPM")

col_data <- metadata %>% column_to_rownames(var = "sample")
col_data$sample <- rownames(col_data)

raw_counts <- round(data.frame(txi_salmon$counts, check.names = FALSE), 0) %>% as.matrix()
colnames(raw_counts) <- rownames(col_data)


se_metadata <- list(metrics = metrics,
                 countsFromAbundance = txi_salmon$countsFromAbundance)

vst <- vst(raw_counts)

dir.create("data")

se <- SummarizedExperiment(assays = list(raw = raw_counts,
                                         tpm = txi_salmon$abundance,
                                         length = txi_salmon$length,
                                         vst = vst),
                           colData = col_data,
                           metadata = se_metadata)
saveRDS(se, "data/bcbio.se.RDS")
