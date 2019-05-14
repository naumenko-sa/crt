init <- function(){
    setwd("~/Desktop/work/wilson/muscle")
    source("~/crt/scripts/crt.utils.R")    
}

plot_mds <- function(file.png,counts,filter=0){

    n_samples <- ncol(counts)
    group <- factor(rep(1, ncol(counts)))
    
    y <- DGEList(counts=counts,
              group=group,
              remove.zeros = T)

    keep <- rowSums(cpm(y) > filter) >= n_samples/2
    y <- y[keep,, keep.lib.sizes = F]
        
    png(file.png, res = 300, width=2000, height=2000)
    plotMDS(y, las=1)
    dev.off()    
}

mds_plots <- function(){
    # all samples
    counts <- read_feature_counts_dir()
    counts <- counts[row.names(counts) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    colnames(counts) <- strsplit2(colnames(counts),"_",fixed=T)[,1]
    plot_mds("mds.all_samples.png",counts)

# removing outliers
counts$WL1101 <- NULL
counts$WL1723 <- NULL
counts$WL1724 <- NULL
plot_mds("mds.no_outliers.png",counts)

# filter 0.1
plot_mds("mds.filter_0.1.png",counts,0.1)
plot_mds("mds.filter_1.png",counts,1)
plot_mds("mds.filter_5.png",counts,5)
plot_mds("mds.filter_less10.png",counts,10)
plot_mds("mds.filter_50.png",counts,50)

setwd("muscle")
counts <- read_feature_counts_dir()
counts <- counts[row.names(counts) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
#colnames(counts) <- strsplit2(colnames(counts),"_",fixed=T)[,1]
plot_mds("mds.muscle.png",counts)
#remove S100 which is fibroblast
counts$S100_26.1.M <- NULL
counts$S49_26.1.M <- NULL
plot_mds("mds.muscle.no_outliers.png",counts)

setwd("../fibro/")
counts <- read_feature_counts_dir()
counts <- counts[row.names(counts) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
#colnames(counts) <- strsplit2(colnames(counts),"_",fixed=T)[,1]
plot_mds("mds.fibro.png",counts)

# only protein coding transcripts
#protein_coding_transcripts = genes_transcripts[genes_transcripts$Ensembl_gene_id %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
#tpm = tpm[row.names(tpm) %in% protein_coding_transcripts$Ensembl_transcript_id,]
#plot_mds("mds.no_outliers.protein_coding_genes.png",tpm)

#tpm$mean = rowMeans(tpm)
#tpm = merge(tpm,genes_transcripts,by.x="row.names",by.y="Ensembl_transcript_id",all.x=T,all.y=F)
#row.names(tpm) = tpm$Row.names
#tpm$Row.names = NULL

#tpm = tpm [order(tpm$mean,decreasing = T),]

#tpm100 = head(tpm,100)
#write.csv(tpm100,"tpm100.csv",quote = F)
}

ranked_expression <- function(){
    rpkms <- feature_counts2rpkm_dir()
    muscular_genes <- get_genes_in_panels()
    rpkms <- rpkms %>% filter(external_gene_name %in% muscular_genes)
    write_excel_csv(rpkms, "expression.in_muscle.csv")
}

diff_expression <- function(){
    prefix="2019-01-28"
    counts <- read_feature_counts_dir()
    counts <- counts[row.names(counts) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    colnames(counts) <- strsplit2(colnames(counts),"_",fixed=T)[,1]
    
    # removing outliers
    counts$WL1101 <- NULL
    counts$WL1723 <- NULL
    counts$WL1724 <- NULL
    
    samples <- colnames(counts)
    n_samples <- length(samples)
    
    group <- factor(c(rep(1, n_samples/2), rep(2, n_samples/2)))
    filter <- 0.5
    
    gene_lengths <- read.delim("~/crt/gene_lengths.txt", stringsAsFactors = F)
    gene_lengths <- gene_lengths[gene_lengths$ENSEMBL_GENE_ID %in% row.names(counts),]
    gene_lengths <- gene_lengths[order(gene_lengths$ENSEMBL_GENE_ID),]
    gene_lengths <- merge(gene_lengths,ensembl_w_description,by.x="ENSEMBL_GENE_ID",by.y="row.names",all.x=T,all.y=F)
    #row.names(gene_lengths)=gene_lengths$ENSEMBL_GENE_ID
    #gene_lengths$ENSEMBL_GENE_ID = NULL
    
    y <- DGEList(counts = counts,
                 group = group,
                 genes = gene_lengths,
                 remove.zeros = T)
    
    rpkm.counts <- rpkm(y)
    
    max_genes <- nrow(counts)
    
    logcpm <- cpm(counts,prior.count=1,log=T)
    t_cpm <- cpm(counts,prior.count=1,log=F)
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
    y=calcNormFactors(y)
    
    #nc=cpm(y,normalized.lib.sizes=F)
    #write.table(nc,"filtered.normalized_counts.txt",col.names=NA)
    
    png(paste0(prefix,".mds.png"),res = 300,width=2000,height=2000)
    plotMDS(y,las=1)
    dev.off()
    
    design=model.matrix(~group)
    
    y=estimateDisp(y,design)
    
    fit <- glmFit(y,design)
    lrt=glmLRT(fit)
    
    efilename=paste0(prefix,".csv")
    de_results = topTags(lrt,n=max_genes,sort.by="PValue",p.value=1,adjust.method = "fdr")
    
    de_results$table$ENSEMBL_GENE_ID = NULL
    de_results$table = merge(de_results$table,rpkm.counts,by.x="row.names",by.y="row.names",all.x=T,all.y=F)
    write.csv(de_results,efilename,quote=T,row.names=F)
    
    de_results = read.csv(efilename, sep="", stringsAsFactors=FALSE)
    s_rownames = row.names(de_results)
    #setnames(de_results,"genes","ensembl_gene_id")
    #de_results = lrt$table
    
    gene_descriptions = read.delim2(paste0("~/cre/ensembl_w_description.txt"), stringsAsFactors=FALSE)
    
    de_results = merge(de_results,gene_descriptions,by.x="genes",by.y="ensembl_gene_id",all.x=T)
    #de_results = rename(de_results,c("Row.names"="ensembl_gene_id"))
    de_results = merge(de_results,counts,by.x = "genes", by.y="row.names",all.x=T)
    #rownames(de_results) = s_rownames
    
    top_genes_cpm = logcpm[de_results$genes,]
    colnames(top_genes_cpm)=paste0(colnames(top_genes_cpm),".log2cpm")
    
    de_results = merge(de_results,top_genes_cpm,by.x = "genes", by.y="row.names",all.x=T)
    de_results = de_results[order(abs(de_results$logFC),decreasing = T),]
    
    colnames(de_results)[1]="Ensembl_gene_id"
    colnames(de_results)[2]="Gene_name"
    de_results$external_gene_name = NULL
    result_file=paste0(prefix,".txt")
    write.table(de_results,result_file,quote=T,row.names=F)
    
    prepare_file_4gsea(counts,samples,prefix)
    
    #plot_heatmap_separate (counts,samples,de_results,prefix)
    plot_heatmap_separate (counts,samples,de_results,paste0(prefix,".top50genes"),50)
    
}
