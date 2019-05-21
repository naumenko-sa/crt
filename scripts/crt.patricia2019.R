init <- function(){
    setwd("~/Desktop/work/patricia")
    library(tidyverse)
    library(edgeR)
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
        protein_coding_genes <- read_csv("protein_coding_genes.list.csv")
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
    counts$ensembl_gene_id <- str_replace(counts$ensembl_gene_id,"\\.\\d+","")
    sample_names <- tibble(sample_name = colnames(counts)) %>% tail(-1)
    
    counts <- inner_join(counts, protein_coding_genes, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
        dplyr::select(ensembl_gene_id, `1130_BD-B175`, `1158_AC-A79`, `1256_TS-T11`, `1275_BK-B225`, 
                      `1365_SD-S169`, `1388_MJ-M219`, `1400_SN-S169`, `1481_TJ-T11`, `2180_CB-C365`, 
                      `2296_GA-G207`, `6034_SA-G207`, `6087_PD-B317G`, `6100_DD-D360`, 
                      `S07_4-1-M`, `S11_8-1-M`, `S12_9-1-M`, `S14_10-1-M`, 
                      `S18_3-1-M_IV-21_9776_RNAseq_polya`, `S19_3-3-M_IV-24_RNAseq_polya`, 
                      `S20_12-1-Mpv`, `S21_12-1-Mvl`, `S22_13-1-M`, `S43_14-1-M`, 
                      `S44_14-2-M`, `S46_17-1-M`, 
                      `S54_33-1-M`, `S56_31-1-M`, `S57_32-1-M`, 
                      `S63_34-1-M`, `S64_35-1-M`, `S65_36-1-M`, 
                      `S67_38-1-M`,
                      `S68_39-1-M`, `S69_39-2-M`, `S70_40-1-M`, `S71_40-2-M`, 
                      `S72_28-1-M`, `S73_18-1-M`, `S77_5-1-M`, `S78_6-1-M`, `S91_43-1-M`, `S92_44-1-M`)
        
         # dplyr::select(Ensembl_gene_id, HCNSM_EV, EN2, EN3, HCNSM_K27M, KN2, KN3)

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
    
    # contrast 8
    # select(Ensembl_gene_id, HCNSM_EV, EN2, EN3, HCNSM_K27M, KN2, KN3)
    
    counts <- column_to_rownames(counts, var = "ensembl_gene_id")
    samples <- colnames(counts)
    n_samples <- length(samples)
    
    #group <- factor(c(rep(1, n_samples/2), rep(2, n_samples/2)))
    group <- factor(c(rep(1,13),rep(2,29)))
    #group <- factor(c(1,1,2,2,2))
    #group <- factor(c(1,2,2,2))
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
    
    keep <- rowSums(cpm(y)>filter) >= n_samples/10
    y <- y[keep,, keep.lib.sizes = F]
    
    # necessary for goana - it uses ENTREZ gene ids
    idfound <- y$genes$ENSEMBL_GENE_ID %in% mappedRkeys(org.Hs.egENSEMBL)
    y <- y[idfound, ]
    egENSEMBL <- toTable(org.Hs.egENSEMBL)
    m <- match (y$genes$ENSEMBL_GENE_ID, egENSEMBL$ensembl_id)
    y$genes$ENTREZID <- egENSEMBL$gene_id[m]
    egSYMBOL <- toTable(org.Hs.egSYMBOL)
    m <- match(y$genes$ENTREZID,egSYMBOL$gene_id)
    y$genes$Symbol <- egSYMBOL$symbol[m]
    
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
    
    y <- estimateDisp(y, design)
    
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit)
    
    prefix <- "contrast_mh_"
    efilename <- paste0(prefix, ".csv")
    de_results <- topTags(lrt, n = max_genes, sort.by = "PValue", p.value = 1, adjust.method = "fdr")$table
    de_results$ENSEMBL_GENE_ID <- NULL
    de_results$ensembl_gene_id <- NULL
    de_results$ENTREZID <- NULL
    de_results$external_gene_name <- NULL
    de_results$gene_description <- NULL
    
    rpkm.counts <- read_csv("rpkms.muscle.csv")
    de_results <- rownames_to_column(as.data.frame(de_results), var = "ensembl_gene_id")
    rpkm.counts$external_gene_name <- NULL
    rpkm.counts$gene_description <- NULL
    
    de_results <- left_join(de_results, rpkm.counts, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
                    left_join(ensembl_w_description, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
                    dplyr::select(ensembl_gene_id, external_gene_name, gene_description, Length,
                               logFC, logCPM, LR, PValue, FDR, one_of(samples)) %>% 
                    filter(FDR < 0.05, abs(logFC)>=2) %>% 
                    arrange(logFC)
    write_excel_csv(de_results, efilename)
    
    prepare_file_4gsea(counts, samples, prefix)
    
    go_analysis(lrt, "contrast8_")
    kegg_analysis(lrt, "contrast8_")
    
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