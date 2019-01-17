###############################################################################
# Bulk RNA-seq project from Michael Tropak : creatinine regulation in mouse
###############################################################################
init <- function(){
    library(edgeR)
    library(GO.db)
    library(genefilter)
    source("~/crt/crt.utils.R")
    setwd("~/Desktop/work/creatinine/")
}

# modified for mouse from crt.gonorazky.naumenko.2018.R
mds_mouse <- function(refresh_files = F){
    
    counts <- read.feature_counts_dir(update = F)
    
    sample_names <- colnames(counts)
    sample_names <- strsplit2(sample_names,"_",fixed=T)[,1]
    sample_labels <- sample_names
    colnames(counts) <- sample_names
    n_samples <- length(sample_names)
    
    group <- factor(c(rep(1,n_samples/2), rep(2,n_samples/2)))    
    v_colors <- c(rep("green",5),rep("red",5))
    
    if (file.exists("~/cre/data/ensembl_w_description.mouse.csv")){
        protein_coding_genes <- read.csv("~/cre/data/ensembl_w_description.mouse.csv", 
                                         stringsAsFactors = F)
        counts <- counts[row.names(counts) %in% protein_coding_genes$ensembl_gene_id,]
    }else{
        print("Please provide protein_coding_genes.list with ENS_GENE_ID")
    }
    
    y <- DGEList(counts = counts, group = group, remove.zeros = T)
    
    png("mds.png", res = 300, width = 2000, height = 2000)
    mds <- plotMDS(y, labels = sample_labels)
    plot(mds,
         col = v_colors,
         pch=19,
         xlab = "MDS dimension 1", 
         ylab = "MDS dimension 2")
    
    legend("topright",
           title="Colors",
           c("Controls",
             "Treatment"
           ),
           fill=c("green",
                  "red"))
    
    dev.off()
    
    png("mds.no_sample6_7.labels.png", res = 300, width = 2000, height = 2000)
    plotMDS(y, labels = sample_labels)
    dev.off()
}

differential_expression <- function(){
    ensembl_w_description.mouse <- read.csv("~/cre/data/ensembl_w_description.mouse.csv",stringsAsFactors = F)
    #counts = read.feature_counts_dir()
    counts <- read.table("raw_counts.txt", stringsAsFactors = F)
    counts <- counts[row.names(counts) %in% ensembl_w_description.mouse$ensembl_gene_id,]
    
    counts <- counts[,c("rh30ctrl7_srh30ctrl7",
                       "rh30ctrl8_srh30ctrl8",
                       "rh30treat7_srh30treat7",
                       "rh30treat8_srh30treat8")]
    
    #group = factor(c(rep(1,ncol(counts))))
    
    sample_names <- colnames(counts)
    sample_names <- strsplit2(sample_names, "_", fixed = T)[,1]
    sample_labels <- sample_names
    colnames(counts) <- sample_names
    n_samples <- length(sample_names)
    
    #group = factor(c(rep(1,n_samples/2), rep(2,n_samples/2)))
    #group = factor(c(1,1,2,2,1,1,2,2))
    group <- factor(c(1,1,2,2))
    
    max_genes <- nrow(counts)
    
    gene_lengths <- read.csv("~/crt/data/gene_lengths.mouse.csv", stringsAsFactors = F)
    gene_lengths <- gene_lengths[gene_lengths$ensembl_gene_id %in% row.names(counts),]
    gene_lengths <- gene_lengths[order(gene_lengths$ensembl_gene_id),]
    #merging by ensembl_gene_id
    gene_lengths <- merge(gene_lengths,
                         ensembl_w_description.mouse,
                         by.x = "ensembl_gene_id",
                         by.y = "ensembl_gene_id", all.x = T, all.y = F)
    #row.names(gene_lengths)=gene_lengths$ENSEMBL_GENE_ID
    #gene_lengths$ENSEMBL_GENE_ID = NULL
    
    y <- DGEList(counts <- counts,
              group = group,
              genes = gene_lengths,
              remove.zeros = T)
    
    rpkm.counts <- rpkm(y)
    
    plotMDS(y)
    #filter - 1 or 0.5
    filter <-  0.5
    keep <- rowSums(cpm(y) > filter) >= n_samples/2
    y <- y[keep, ,keep.lib.sizes = F]
    
    y$samples$lib.size <- colSums(y$counts)
    
    
    #normalization for RNA composition (2.7.3)
    y=calcNormFactors(y)
    
    #nc=cpm(y,normalized.lib.sizes=F)
    #write.table(nc,"filtered.normalized_counts.txt",col.names=NA)
    
    plotMDS(y,las=1)
    
    design <- model.matrix(~group)
    
    y <- estimateDisp(y,design)
    
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit)
    
    #o=order(lrt$table$PValue)
    #cpm(y)[o[1:10],]
    #write.table(cpm(y)[o[1:12566],],"allgenes.cpm")
    
    de_results <- topTags(lrt, n=max_genes, sort.by="PValue", p.value=1, adjust.method = "fdr")
    
    de_results$table$ensembl_gene_id <-  NULL
    de_results$table <- merge(de_results$table, 
                              rpkm.counts, 
                              by.x = "row.names", 
                              by.y = "row.names",
                              all.x = T, all.y = F)
    colnames(de_results$table)[1] = "Ensembl_gene_id"
    write.csv(de_results, "de.csv", quote = T, row.names = F)
}

# combining counts for different tissues from the article 
# https://www.nature.com/articles/sdata2017185#data-records
# and Tropak's samples
combine_counts = function()
{
    counts <- read.feature_counts_dir(update = F)
    mouse_counts <- read.delim("mouse_counts.txt", stringsAsFactors = F)
    counts <- merge(counts, mouse_counts, by.x = "row.names", by.y = "row.names",
                   all.x = T, all.y = F)
    row.names(counts) <- counts$Row.names
    counts$Row.names <- NULL
    
    if (file.exists("~/cre/data/ensembl_w_description.mouse.csv"))
    {
        protein_coding_genes <- read.csv("~/cre/data/ensembl_w_description.mouse.csv", 
                                         stringsAsFactors = F)
        counts <- counts[row.names(counts) %in% protein_coding_genes$ensembl_gene_id,]
    }else{
        print("Please provide protein_coding_genes.list with ENS_GENE_ID")
    }
    
    counts <- na.omit(counts)
    
    n_samples <- ncol(counts)
    
    group <- factor(rep(1,n_samples))    
    
    y <- DGEList(counts = counts, group = group, remove.zeros = T)
    
    png("mds.tissues.png",res=300, width = 2000, height = 2000)
    plotMDS(y, labels = colnames(counts))
    dev.off()
}
