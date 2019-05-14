# MDS for protein coding genes
# run: qsub ~/crt.mds.pbs -v refresh=TRUE
# or Rscript ~/crt/crt.mds.R TRUE
mds_work <- function(update = F, top_genes = 500, mds_samples){
    #update <- T
    counts <- feature_counts2rpkm_dir(update = update)
    print("Subsetting protein coding genes ...")
    counts <- counts %>% filter(ensembl_gene_id %in% protein_coding_genes$ensembl_gene_id)
    counts$ensembl_gene_id <- NULL
    counts$external_gene_name <- NULL
    counts$gene_description <- NULL
    
    sample_names <- tibble(sample_name = colnames(counts))
    samples <- read_csv("mds_samples.csv")
    
    print("Removing zeroes ...")
    group <- factor(c(rep(1,ncol(counts))))
    y <- DGEList(counts = counts, group = group, remove.zeros = T)
    
    print("Plotting ...")
    png("mds.png", res=300, width = 2000, height = 2000)
    print(paste0("Top genes: ", top_genes))
    mds <- plotMDS(y, top = top_genes)
    
    v_colors <- left_join(sample_names, samples, by = "sample_name") %>% select(color) %>% unlist(use.names = F)
    plot(mds,
         col = v_colors,
         pch = 19,
         xlab = "MDS dimension 1", 
         ylab = "MDS dimension 2")
    dev.off()
    
    png("mds.legend.png", res=300, width=2000, height=2000)
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("topleft",
           title = "Tissue",
           c("GTEx blood",
             "NebNext blood",
             "RNADirect blood",
             "Primary fibroblasts",
             "GTEx transformed fibroblasts",
             "RNADirect fibroblasts",
             "Transdifferentiated myotubes",
             "Muscle",
             "GTEx muscle",
             "RNADirect muscle"
           ),
           fill = c("cornflowerblue",
                    "cyan",
                    "navyblue",
                    "orange",
                    "yellow",
                    "darkviolet",
                    "red",
                    "chartreuse",
                    "darkgreen",
                    "darkkhaki"))
    
    dev.off()
    v_labels <- left_join(sample_names, samples, by="sample_name") %>% 
        select(sample_label) %>% unlist(use.names = F)
    
    png("mds.labels.png", res = 300, width = 2000, height = 2000)
    plotMDS(y, labels = v_labels, cex = 0.2)
    dev.off()
}

###############################################################################
args <- commandArgs(trailingOnly = T)
if (length(args) == 0 || args[1] == "--help"){
    cat("Usage: Rscript function_name function_args\n")
    cat("Available functions:\n")
    cat("mds_work [refresh_files=TRUE|F] [top_genes=500]\n")
}else{
    cat(paste0("Running function: ", args[1],"\n"))
    source("~/crt/scripts/crt.utils.R")
    fcn_name <- args[1]
    if (fcn_name == "mds_work"){
        fcn <- get(args[1])
        fcn(update = as.logical(args[2]), top_genes = as.numeric(args[3]))    
    }
}
###############################################################################