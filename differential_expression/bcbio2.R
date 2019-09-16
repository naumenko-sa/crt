# https://github.com/vbarrera/bcbioLite/blob/master/R/readBCBIO.R

# Rscript bcbio2se.R /bcbio/result/final/project

library(tidyverse)

project_dir <- "/home/sergey/cluster/2019-08-29_9way_bcbio/final/2019-08-22_project"

project_dir <- args[1]

metadata <- read_csv(file.path(project_dir, "metadata.csv"))

metrics <- rio::import(file.path(projectDir,"multiqc","multiqc_data","multiqc_bcbio_metrics.txt"))
metrics <- janitor::clean_names(metrics,case = "small_camel")
sampleDirs = file.path(projectDir,"..", rownames(sampleMetadata))
tx2genes_file <- file.path(projectDir,"tx2gene.csv")
salmon_files = file.path(sampleDirs, "salmon", "quant.sf")
sf_files = salmon_files
names(sf_files) = rownames(sampleMetadata)
tx2gene = read.table(tx2genes_file, sep=",", row.names=NULL, header=FALSE)
txi.salmon = tximport::tximport(sf_files, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")
rawCounts = round(data.frame(txi.salmon$counts, check.names=FALSE),0) %>% as.matrix()
colData <- sampleMetadata %>% as.data.frame() %>% .[colnames(rawCounts),
                                                        , drop = FALSE] %>% tibble::rownames_to_column() %>% dplyr::mutate_all(as.factor) %>%
mutate_all(droplevels) %>% tibble::column_to_rownames() %>% as("data.frame")

metadata <- list(metrics = metrics, countsFromAbundance = txi.salmon$countsFromAbundance)
se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = rawCounts,
                                                                   tpm = txi.salmon$abundance,
                                                                   length = txi.salmon$length),
                                                     colData = colData,
                                                     metadata = metadata)
    invisible(se)