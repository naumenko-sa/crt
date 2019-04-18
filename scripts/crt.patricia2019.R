init <- function(){
    setwd("~/Desktop/work/patricia")
    library(tidyverse)
    library(edgeR)
    ensembl_w_description <- read.delim2("~/cre/data/ensembl_w_description.txt", stringsAsFactors = F)
    source("~/crt/scripts/crt.utils.R")
}

plot_mds <- function(){
    # ENSG00000187474.4 row was corrupted - 26 columns
    counts <- read_csv("raw_counts_joined.csv")
    # remove suffix
    counts$Ensembl_gene_id <- str_replace(counts$Ensembl_gene_id,"\\.\\d+","")

    sample_names <- tibble(sample_name = colnames(counts)) %>% tail(-1)

    samples <- read_csv("samples.csv")
    legend <- read_csv("legend.csv")

    print("Subsetting protein coding genes ...")
    #a file with ENSEMBL IDs
    if (file.exists("protein_coding_genes.list")){
        protein_coding_genes <- read.csv("protein_coding_genes.list", sep="", stringsAsFactors=F)
        counts <- inner_join(counts, protein_coding_genes, by = c("Ensembl_gene_id" = "ENS_GENE_ID"))
    }else{
        print("Please provide protein_coding_genes.list with ENS_GENE_ID")
    }

    print("Removing zeroes ...")
    counts$Ensembl_gene_id <- NULL
    group <- factor(c(rep(1, ncol(counts))))

    #counts <- drop_na(counts)
    y <- DGEList(counts = counts, group = group, remove.zeros = T)
    
    top_genes = 1000
    print("Plotting ...")
    png("mds.png", res=300, width=2000, height=2000)
    print(paste0("Top genes: ", top_genes))
    mds <- plotMDS(y, top = top_genes)

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
           title = "Cell type", 
           legend$sample_type, fill = legend$color)

    dev.off()
    #v_labels <- left_join(sample_names, samples, by="sample_name") %>% select(sample_label) %>% unlist(use.names = F)
    v_labels = sample_names$sample_name
    png("mds.labels.png", res = 300, width = 2000, height = 2000)
    plotMDS(y, labels = v_labels, cex = 0.8)
    dev.off()
}

rpkm_table <- function(){
    counts <- read_csv("raw_counts.csv")
    # remove suffix
    counts$Ensembl_gene_id <- str_replace(counts$Ensembl_gene_id,"\\.\\d+","")
    
    protein_coding_genes <- read_csv("protein_coding_genes.list")
    counts <- inner_join(counts, protein_coding_genes, by = c("Ensembl_gene_id" = "ENS_GENE_ID")) 
    counts <- column_to_rownames(counts, var = "Ensembl_gene_id")
    
    sample_names <- tibble(sample_name = colnames(counts))
    
    n_samples <- nrow(sample_names)
    
    group <- factor(c(rep(1, n_samples)))
    filter <- 0.5
    
    gene_lengths <- read.delim("~/crt/data/gene_lengths.txt", stringsAsFactors = F)
    gene_lengths <- gene_lengths[gene_lengths$ENSEMBL_GENE_ID %in% row.names(counts),]
    gene_lengths <- gene_lengths[order(gene_lengths$ENSEMBL_GENE_ID),]
    gene_lengths <- merge(gene_lengths, ensembl_w_description, by.x = "ENSEMBL_GENE_ID",
                          by.y = "row.names", all.x = T, all.y = F)
    #row.names(gene_lengths)=gene_lengths$ENSEMBL_GENE_ID
    #gene_lengths$ENSEMBL_GENE_ID = NULL
    
    counts <- counts[row.names(counts) %in% gene_lengths$ENSEMBL_GENE_ID,]
    
    y <- DGEList(counts = counts,
                 group = group,
                 genes = gene_lengths,
                 remove.zeros = T)
    
    rpkm_counts <- as_tibble(rownames_to_column(as.data.frame(rpkm(y))))
    write_excel_csv(rpkm_counts, "rpkms.csv")
}

#gsea_input_tsv <- "gsea/EV_MUT.GM.no_EG2/result/MUELLER_METHYLATED_IN_GLIOBLASTOMA.xls"
#result_file <-  "ev_mut.gs.no_eg2.mueller.csv"
add_rpkm_to_gsea <- function(gsea_input_tsv, result_file){
    rpkm_counts <- read_csv("rpkms.csv")
    gsea_result <- read_tsv(gsea_input_tsv) %>% 
        select(PROBE, `DESCRIPTION<br>(from dataset)`, `CORE ENRICHMENT`)
    colnames(gsea_result) <- c("gene","ensembl_gene_id", "core_enrichment")
    
    gsea_result <- left_join(gsea_result, rpkm_counts, by = c("ensembl_gene_id" = "rowname"))
    write_excel_csv(gsea_result, result_file)
}

differential_expression <- function()
{
    counts <- read_csv("raw_counts.csv")
    # remove suffix
    counts$Ensembl_gene_id <- str_replace(counts$Ensembl_gene_id,"\\.\\d+","")
    sample_names <- tibble(sample_name = colnames(counts)) %>% tail(-1)
    
    protein_coding_genes <- read_csv("protein_coding_genes.list")
    counts <- inner_join(counts, protein_coding_genes, by = c("Ensembl_gene_id" = "ENS_GENE_ID")) %>% 
        select(Ensembl_gene_id, GLINS1, HCNSM_K27M, KN2, KN3)
    
    # contrast 1
    # select(Ensembl_gene_id, HCGM_WT, WG2, WG3, HCGM_K27M, KG2, KG3)
    
    # contrast 2
    # select(Ensembl_gene_id, WG2, WG3, KG2, KG3)
    
    # contrast 3
    # select(Ensembl_gene_id, HCGM_EV, EG2, EG3, HCGM_K27M, KG2, KG3)
    
    # contrast 4
    # select(Ensembl_gene_id, HCGM_EV, EG3, HCGM_K27M, KG2, KG3)
    
    # contrast 5
    # select(Ensembl_gene_id, HCGM_K27M, KG2, KG3, HCNSM_K27M, KN2, KN3)
    
    # contrast 6
    # select(Ensembl_gene_id, KG2, KG3, KN2, KN3)
    
    # contrast 7
    # select(Ensembl_gene_id, GLINS1, HCNSM_K27M, KN2, KN3)
    
    counts <- column_to_rownames(counts, var = "Ensembl_gene_id")
    samples <- colnames(counts)
    n_samples <- length(samples)
    
    #group <- factor(c(rep(1, n_samples/2), rep(2, n_samples/2)))
    #group <- factor(c(1,1,2,2,2))
    group <- factor(c(1,2,2,2))
    filter <- 1
    
    gene_lengths <- read.delim("~/crt/data/gene_lengths.txt", stringsAsFactors = F)
    gene_lengths <- gene_lengths[gene_lengths$ENSEMBL_GENE_ID %in% row.names(counts),]
    gene_lengths <- gene_lengths[order(gene_lengths$ENSEMBL_GENE_ID),]
    gene_lengths <- merge(gene_lengths, ensembl_w_description, by.x = "ENSEMBL_GENE_ID",
                          by.y = "row.names", all.x = T, all.y = F)
    #row.names(gene_lengths)=gene_lengths$ENSEMBL_GENE_ID
    #gene_lengths$ENSEMBL_GENE_ID = NULL
    
    counts <- counts[row.names(counts) %in% gene_lengths$ENSEMBL_GENE_ID,]
    
    y <- DGEList(counts = counts,
                 group = group,
                 genes = gene_lengths,
                 remove.zeros = T)
    
    rpkm.counts <- rpkm(y)
    
    max_genes <- nrow(counts)
    
    logcpm <- cpm(counts, prior.count = 1, log = T)
    t_cpm <- cpm(counts, prior.count = 1, log = F)
    logcpm <- logcpm[,samples]
    t_cpm <- t_cpm[,samples]
    
    plotMDS(y)
    
    keep <- rowSums(cpm(y)>filter) >= n_samples/2
    y <- y[keep,, keep.lib.sizes = F]
    
    #necessary for goana
    #idfound = y$genes$genes %in% mappedRkeys(org.Hs.egENSEMBL)
    #y = y[idfound,]
    
    #egENSEMBL=toTable(org.Hs.egENSEMBL)
    #m = match (y$genes$genes,egENSEMBL$ensembl_id)
    #y$genes$EntrezGene = egENSEMBL$gene_id[m]
    #egSYMBOL = toTable(org.Hs.egSYMBOL)
    #m = match (y$genes$EntrezGene,egSYMBOL$gene_id)
    #y$genes$Symbol = egSYMBOL$symbol[m]
    
    #remove duplications - just 1 gene in this dataset
    #first order by counts to remove duplicated names with 0 counts
    #o = order(rowSums(y$counts),decreasing = T)
    #y = y[o,]
    #d = duplicated(y$genes$Symbol)
    #dy = y[d,]$genes
    #y = y[!d,]
    #nrow(y)
    
    #y$samples$lib.size = colSums(y$counts)
    #rownames(y$counts) = y$genes$EntrezGene 
    #rownames(y$genes) = y$genes$EntrezGene
    #y$genes$EntrezGene = NULL
    
    #normalization for RNA composition (2.7.3)
    y <- calcNormFactors(y)
    
    #nc=cpm(y,normalized.lib.sizes=F)
    #write.table(nc,"filtered.normalized_counts.txt",col.names=NA)
    
    png(paste0(prefix,".mds.png"),res = 300,width=2000,height=2000)
    plotMDS(y,las=1)
    dev.off()
    
    design <- model.matrix(~group)
    
    y <- estimateDisp(y,design)
    
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit)
    
    prefix <- "2019-03-14"
    efilename <- paste0(prefix, ".csv")
    de_results <- topTags(lrt, n = max_genes, sort.by = "PValue", p.value = 1, adjust.method = "fdr")
    
    de_results$table$ensembl_gene_id <- NULL
    de_results$table$external_gene_name <- NULL
    de_results$table$Gene_description <- NULL
    de_results$table$ENSEMBL_GENE_ID <- NULL
    
    de_results$table <- merge(de_results$table, rpkm.counts, by.x = "row.names", by.y = "row.names", 
                              all.x = T, all.y = F)
    colnames(de_results$table)[1] = "ensembl_gene_id"
    write.csv(de_results$table, efilename, quote = T, row.names = F)
    
    de_results = read.csv(efilename, stringsAsFactors = F)
    de_results = merge(de_results, ensembl_w_description,
                       by.x = "ensembl_gene_id", by.y = "row.names", all.x = T)
    de_results <- de_results[c("ensembl_gene_id", "external_gene_name", "Gene_description", "Length",
                               "logFC", "logCPM", "LR", "PValue", "FDR", samples)]
    de_results <- de_results[de_results$FDR < 0.05,]
    de_results <- de_results[abs(de_results$logFC)>=2,]
    de_results <- de_results[order(de_results$logFC),]
    write.csv(de_results, efilename, quote = T, row.names = F)
    
    prepare_file_4gsea(counts, samples, prefix)
    
    #plot_heatmap_separate (counts,samples,de_results,prefix)
    #plot_heatmap_separate (counts,samples,de_results,paste0(prefix,".top50genes"),50)
}

merge_with_nishani2019 <- function(){
    nishani2019 <- read_delim("LGK_counts.txt", delim = " ")

    counts <- read_csv("raw_counts.csv")
    tmp <- as_tibble(str_split(counts$Ensembl_gene_id, "\\.", n=2, simplify = T))
    counts$Ensembl_gene_id <- tmp$V1
    
    counts <- left_join(counts, nishani2019, by = c("Ensembl_gene_id" = "ensembl_gene_id"))
    
    write_excel_csv(counts, "raw_counts_joined.csv")
    
}