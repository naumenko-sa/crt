# MDS for protein coding genes
# run: qsub ~/crt.mds.pbs -v refresh=TRUE
# or Rscript ~/crt/crt.mds.R TRUE
mds_work <- function(update = F){
    update <- T
    counts <- read_feature_counts_dir(update = update)
    
    sample_names <- tibble(sample_name = colnames(counts))
    
    samples <- read_csv("mds_samples.csv")
    
    print("Subsetting protein coding genes ...")
    #a file with ENSEMBL IDs
    if (file.exists("protein_coding_genes.list")){
        protein_coding_genes <- read.csv("protein_coding_genes.list", sep="", stringsAsFactors=F)
        counts <- counts[row.names(counts) %in% protein_coding_genes$ENS_GENE_ID,]
    }else{
        print("Please provide protein_coding_genes.list with ENS_GENE_ID")
    }
    
    print("Removing zeroes ...")
    group <- factor(c(rep(1,ncol(counts))))
    y <- DGEList(counts = counts, group = group,remove.zeros = T)
    
    print("Plotting ...")
    png("mds.png", res=300, width=2000, height=2000)
    mds <- plotMDS(y)
    
    v_colors <- left_join(sample_names, samples, by="sample_name") %>% select(color) %>% unlist(use.names = F)
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
    v_labels <- left_join(sample_names, samples, by="sample_name") %>% select(sample_label) %>% unlist(use.names = F)
    png("mds.labels.png", res = 300, width = 2000, height = 2000)
    plotMDS(y, labels = v_labels)
    dev.off()
}

source("~/crt/crt.utils.R")
args = commandArgs(trailingOnly = T)
print(args[1])
mds_work(update = as.logical(args[1]))

