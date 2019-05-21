###############################################################################
# Functions for RNA-seq article Gonorazky.Naumenko.et_al.Dowling.2019
# expression outlier detection in ~/crt/crt.expression.outliers.R
# 1. Clustering, MDS, PCA plots
# 2. Heatmaps
# 3. Expression
# 4. Coverage
# 5. Splicing
###############################################################################
install <- function(){
    source("http://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
}

init <- function(){
    library(edgeR)
    library(RColorBrewer)
    library(org.Hs.eg.db)
    library(plyr)
    library(VennDiagram)    
    library(genefilter)
    library(reshape2)
    library(tidyverse)
    source("~/crt/scripts/crt.utils.R")
    source("~/bioscripts/gene_panels/genes.R")
    setwd("~/Desktop/work")
    #gene_lengths = read.delim("~/Desktop/project_muscular/reference/gene_lengths.txt", stringsAsFactors=F, row.names=1)
}

###############################################################################
# 1. clustering, MDS, PCA plots
###############################################################################
fig2A_mds <- function()
{
    plot_mds("rpkms.protein_coding.csv", "mds_samples.csv")    
    png("mds.legend.png", res = 300, width = 2000, height = 2000)
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
    plot_mds("rpkms.protein_coding.csv", "mds_samples.csv")
}

# fig 1B - expression profiles in muscular samples vs sample age
# input = TableS1.Samples
fig2B.mds_plot_colored_by_muscle_age <- function(){
    setwd("~/Desktop/work/expression/muscle")
  
    muscle_samples = read.csv("TableS1.Samples.csv", stringsAsFactors=F)
    muscle_samples$Sample_name = gsub("-",".",muscle_samples$Sample_name)
    muscle_samples$Bioinf_sample_id = paste0(muscle_samples$Bioinf_sample_id,"_",muscle_samples$Sample_name)
    muscle_samples = muscle_samples[muscle_samples$Tissue_type == "Muscle",]
  
    refresh_files=T
  
    print("Reading counts ...")
    counts = read.feature_counts_dir(update=refresh_files)
    #group = factor(c(rep(1,ncol(counts))))
  
    #https://gist.github.com/jdblischak/11384914
    cpm_log <- cpm(counts, log = TRUE)
    median_log2_cpm <- apply(cpm_log, 1, median)
    hist(median_log2_cpm)
    expr_cutoff <- -1
    abline(v = expr_cutoff, col = "red", lwd = 3)
    sum(median_log2_cpm > expr_cutoff)
    data_clean <- counts[median_log2_cpm > expr_cutoff, ]
    
    
    cpm_log <- cpm(data_clean, log = TRUE)
    heatmap(cor(cpm_log))
    
    
    png("fig1B.pca1_2.png",res=300,width=2000,height=2000)
    pca <- prcomp(t(cpm_log), scale. = TRUE)
    plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
    text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
    summary(pca)
    dev.off()
    
    test = as.data.frame(pca$x[,1])
    test = cbind(test, muscle_samples$Age)
    
    png("fig1B.pca1_vs_muscle_age.png",res=300,width=2000,height=2000)
    plot(test$`pca$x[, 1]` ~ test$`muscle_samples$Age`)
    model <- lm(test$`pca$x[, 1]` ~ test$`muscle_samples$Age`)
    abline(model, col = "red")
    dev.off()
    summary(model)
    coef(model)
    
    sample_names = colnames(counts)
  
    i=1
  
    group = factor(rep(1,length(sample_names)))    
  
    v_colors = c()
    for (i in 1:nrow(muscle_samples))
    {
        sample_id = muscle_samples$Bioinf_sample_id[i]
        #v_colors[i] = muscle_samples$Color[i]
        age = muscle_samples$Age[i]
        
        #iage = as.integer(age)
        #if (iage == 0){
        #    clr = "cornflowerblue"
        #}else (age == 2){
        #    clr = "yellow"
        #}else if (iage <=10){
        #    clr = "orange"
        #}else if (iage <=20){
        #    clr = "red"
        #}else if (iage <=35){
        #    clr = "chartreuse"
        #}else{
        #    clr = "darkgreen"
        #}
    
        #v_colors[i] = clr
        print(paste(age,clr))
    }
  
    print("Subsetting protein coding genes ...")
    #a file with ENSEMBL IDs
    if (file.exists("~/Desktop/reference_tables/protein_coding_genes.list"))
    {
        protein_coding_genes <- read.csv("~/Desktop/reference_tables/protein_coding_genes.list", sep="", stringsAsFactors=FALSE)
        counts = counts[row.names(counts) %in% protein_coding_genes$ENS_GENE_ID,]
    }else{
        print("Please provide protein_coding_genes.list with ENS_GENE_ID")
    }
  
    print("Removing zeroes ...")
  
    y=DGEList(counts=counts,group=group,remove.zeros = T)
  
    
    
    png("mds.muscle_age.png",res=300,width=2000,height=2000)
  
    print("Plotting ...")
    
    set.seed(1)
    par(xpd=F)
  
    mds = plotMDS(y,
                labels=sample_names
    )
    plot(mds,
       col = v_colors,
       pch=19,
       xlab = "MDS dimension 1", 
       ylab = "MDS dimension 2")
  
    par(xpd=T)
    legend(0,0,
         title="Tissue age",
         c("10 days",
           "<=2 y",
           "2-10y",
           "10-20y",
           "20-35y",
           ">35y"
         ),
         fill=c("cornflowerblue",
                "yellow",
                "orange",
                "red",
                "chartreuse",
                "darkgreen"))
  
    dev.off()
  
    png("mds.muscle_age.labels.png",res=300,width=2000,height=2000)
    plotMDS(y, labels=sample_names)
    dev.off()
  
}

###############################################################################
# 2. Heatmaps
###############################################################################


###############################################################################
# 3. Expression
###############################################################################
tableS5.get_1rpkm_genes <- function(rpkms,samples,use_sample_names=T){
    rpkms = rpkms[row.names(rpkms) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    
    if (use_sample_names == T)
    {
        colnames(rpkms) = gsub("S.._","",colnames(rpkms))
        rpkms = subset(rpkms,select = samples)
    }
    
    rpkms$mean = rowMeans(rpkms)
    rpkms = rpkms[rpkms$mean > 1,]
    return(rpkms)
}

# intersect row names from two dataframes
tableS5.intersect = function(tissue1,tissue2){
    tissue1.rpkms = get(paste0('rpkms.',tissue1))
    tissue2.rpkms = get(paste0('rpkms.',tissue2))
    print(paste0(tissue1,':',nrow(tissue1.rpkms)))
    print(paste0(tissue2,':',nrow(tissue2.rpkms)))
    print(paste0(tissue1,' vs ',tissue2,':',length(intersect(rownames(tissue1.rpkms),rownames(tissue2.rpkms)))))
}

# how many genes are covered at 10X in blood from the muscular panel?
# using output of bam.coverage.bamstats05.sh
fig1C.genes_covered_at_10x_in_panel = function()
{
    setwd("~/Desktop/work/coverage")
    files = list.files(".","\\.coverage$")
    #samples = unlist(read.table("samples.txt", stringsAsFactors=F))
    
    coverage = read.delim(files[1],header=T,stringsAsFactors = F)
    coverage = coverage[,c("gene","avg")]
    colnames(coverage)[2]=files[1]
    
    for (file in tail(files,-1))
    {
        sample_coverage = read.delim(file,header=T,stringsAsFactors = F)
        sample_coverage = sample_coverage[,c("gene","avg")]
        colnames(sample_coverage)[2]=file
        coverage = cbind(coverage,sample_coverage[2])
    }
    row.names(coverage) = coverage$gene
    coverage$gene=NULL
    
    coverage$mean = rowMeans(coverage)
    
    coverage = coverage[coverage$mean >= 10,]
    
    print(paste0(
        nrow(coverage),
        " genes out of ",
        length(gene_list),
        " covered at 10X"
    ))
}

# how many genes expressed >1RPKM in blood in muscular gene panels?
fig1C.genes_expressed_at_1rpkm_in_panel = function()
{
    # prepare gtex blood table
    rpkms.gtex_blood = read.csv("rpkms.gtex_blood.txt", sep="", stringsAsFactors = F)
    rpkms.gtex_blood$external_gene_name = NULL
    rpkms.gtex_blood = tableS5.get_1rpkm_genes(rpkms.gtex_blood,NULL,F)
    
    gene_list = c()
    for (gene_panel_name in panel_list)
    {
        gene_panel = get(gene_panel_name)
        gene_list = unique(c(gene_list,gene_panel))
    }
    
    mart = init_mart()
    
    ensembl_ids = get_ensemble_gene_ids_by_gene_names(gene_list)
    
    rpkms.gtex_blood = rpkms.gtex_blood[rownames(rpkms.gtex_blood) %in% ensembl_ids$ensembl_gene_id,]
    
    print(paste0(
            nrow(rpkms.gtex_blood),
            " genes out of ",
            length(gene_list),
            " expressed at 1RPKM"
            ))
}

# among genes expressed in muscle at 1RPKM, what is % expressed in myo at > 1RPKM
# based on trios (muscle + myo + fibro) and GTEx blood controls
# fig1C - genes expressed at 1RPKM in 4 tissues
# data from tableS5
supplementary_table_3.genes_expressed_at_1rpkm = function()
{
    # prepare muscle table
    samples.muscle = paste0(trios,".M")
    rpkms.muscle = read.csv("rpkms.muscle.txt", sep="", stringsAsFactors=F)
    
    # average two muscle tissues in 12.1
    rpkms.muscle$"12.1.M" = (rpkms.muscle$S20_12.1.Mpv + rpkms.muscle$S21_12.1.Mvl)/2
    rpkms.muscle$S20_12.1.Mpv = NULL
    rpkms.muscle$S21_12.1.Mvl = NULL
    
    rpkms.muscle = tableS5.get_1rpkm_genes(rpkms.muscle,samples.muscle)
    
    # prepare myo table
    samples.myo = paste0(trios,".Myo")
    rpkms.myo = read.csv("rpkms.myo.txt", sep="", stringsAsFactors=F)
    rpkms.myo = tableS5.get_1rpkm_genes(rpkms.myo,samples.myo)
    
    # prepare fibro table
    samples.fibro = paste0(trios,".F")
    rpkms.fibro = read.csv("rpkms.fibro.txt", sep="", stringsAsFactors = F)
    rpkms.fibro = tableS5.get_1rpkm_genes(rpkms.fibro,samples.fibro)
    
    #rpkms.fibro$external_gene_name=NULL
    #rpkms.fibro = tableS5.get_1rpkm_genes(rpkms.fibro,NULL,F)
    
    # prepare gtex blood table
    rpkms.gtex_blood = read.csv("rpkms.gtex_blood.txt", sep="", stringsAsFactors = F)
    rpkms.gtex_blood$external_gene_name = NULL
    rpkms.gtex_blood = tableS5.get_1rpkm_genes(rpkms.gtex_blood,NULL,F)
    
    #prepare gtex fibro table
    rpkms.gtex_fibro = read.csv("rpkms.gtex_fibro.txt", sep="", stringsAsFactors = F)
    rpkms.gtex_fibro$external_gene_name = NULL
    rpkms.gtex_fibro = tableS5.get_1rpkm_genes(rpkms.gtex_fibro,NULL,F)
    
    print("genes expressed in at 1RPKM")
    tissues = c("muscle","myo","fibro","gtex_blood") 
    tissue_labels = c("Muscle","Transdifferentiated \nmyotubes","Primary fibroblasts","Gtex blood")

    png("fig1C.genes_expressed_at_1_rpkm.png",res=300,width=2000,height=2000)
    grid.draw(venn.diagram(list(rownames(rpkms.muscle),
                      rownames(rpkms.myo),
                      rownames(rpkms.fibro),
                      rownames(rpkms.gtex_blood)),NULL,
                    category.names = tissue_labels,
                    print.mode = c("percent","raw"),
                    sigdigs = 2,
                    margin = 0.1))
    dev.off()
    
    # just checking
    tissue_pairs = t(combn(tissues,2))
    for (i in 1:nrow(tissue_pairs))
    {
        tableS5.intersect(tissue_pairs[i,1],tissue_pairs[i,2])
    }
}

# expression of SMN1 and SMN2 isoforms
expression_smn <- function(){
    setwd("~/Desktop/work/SMN/")

    read_gene_kallisto_dir("SMN1")
    read_gene_kallisto_dir("SMN2")
    
    S93.smn1 <- get_isoform_expression_kallisto("S93_45-1-F.tsv", "SMN1", "S93_45-1-F.tsv")
    S93.smn2 <- get_isoform_expression_kallisto("S93_45-1-F.kallisto.tsv","SMN2","S93")
    
    S94.smn1 <- get_isoform_expression_kallisto("S94_46-1-F.kallisto.tsv","SMN1","S94")
    S94.smn2 <- get_isoform_expression_kallisto("S94_46-1-F.kallisto.tsv","SMN2","S94")

    smn2 <- genes_transcripts[genes_transcripts$external_gene_name=="SMN2",]
    
    S93 <- read.delim("S93_45-1-F.kallisto.tsv",header=T,stringsAsFactors = F)
    S94 <- read.delim("S94_46-1-F.kallisto.tsv",header=T,stringsAsFactors = F)
    
    S93.smn1 <- S93[S93$target_id %in% smn1$Ensembl_transcript_id,]
    S93.smn2 <- S93[S93$target_id %in% smn2$Ensembl_transcript_id,]
    
    S94.smn1 <- S94[S94$target_id %in% smn1$Ensembl_transcript_id,]
    S94.smn2 <- S94[S94$target_id %in% smn2$Ensembl_transcript_id,]
    
    S93.smn1 <- S93.smn1[,c("target_id","tpm")]
    colnames (S93.smn1) = c("ensembl_transcript_id","S93")
    S93.smn2 = S93.smn2[,c("target_id","tpm")]
    colnames (S93.smn2) = c("ensembl_transcript_id","S93")
    
    S94.smn1 = S94.smn1[,c("target_id","tpm")]
    colnames(S94.smn1) = c("ensembl_transcript_id","S94")
    S94.smn2 = S94.smn2[,c("target_id","tpm")]
    colnames(S94.smn2) = c("ensembl_transcript_id","S94")
    
    smn1 = cbind(S93.smn1,S94.smn1$S94)
    colnames(smn1)[3]="S94"
    smn2 = cbind(S93.smn2,S94.smn2$S94)
    colnames(smn2)[3]="S94"
    
    write.csv(smn1,"smn1.csv",row.names = F)
    write.csv(smn2,"smn2.csv",row.names = F)
}

dexpression <- function()
{
    #counts <- read_csv
    
    counts = counts[row.names(counts) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    
    #counts = counts[,c("S84_GOR1.1.F.Fiech","S85_GOR1.2.F.Fieja1","S86_GOR1.3.F.Fiejo2","S87_GOR1.4.F.Avoel",
    #                   "S08_5.1.F","S09_6.1.F","S10_7.1.F","S25_15.1.F","S26_15.2.F",
    #                   "S27_15.3.F","S28_17.1.F","S29_20.1.F","S30_27.1.F","S31_16.1.F",
    #                   "S33_18.1.F","S35_21.1.F","S36_22.1.F","S37_12.1.F","S38_24.1.F",
    #                   "S39_25.1.F","S40_28.1.F","S41_29.1.F","S42_30.1.F","S59_14.1.F",
    #                   "S60_14.2.F","S61_3.2.F","S62_9.1.F","S75_26.1.F","S76_41.1.F")]
    
    counts = counts[,c("S84_GOR1.1.F.Fiech","S85_GOR1.2.F.Fieja1","S86_GOR1.3.F.Fiejo2","S87_GOR1.4.F.Avoel")]
    
    #attach(counts)
    #counts$fibro_average = (S08_5.1.F+S09_6.1.F+S10_7.1.F+S25_15.1.F+S26_15.2.F+
    #    +S27_15.3.F + S28_17.1.F + S29_20.1.F + S30_27.1.F + S31_16.1.F + 
    #    S33_18.1.F +  S35_21.1.F + S36_22.1.F + S37_12.1.F + S38_24.1.F + 
    #    S39_25.1.F + S40_28.1.F + S41_29.1.F + S42_30.1.F + S59_14.1.F + 
    #        S60_14.2.F + S61_3.2.F + S62_9.1.F + S75_26.1.F + S76_41.1.F)/25
    
    #counts = counts[,c("S84_GOR1.1.F.Fiech","S85_GOR1.2.F.Fieja1","S86_GOR1.3.F.Fiejo2","S87_GOR1.4.F.Avoel","fibro_average")]
    
    samples = colnames(counts)
    n_samples = length(samples)
    #group=factor(c(rep(1,n_samples/2),rep(2,n_samples/2)))
    #group = factor(c(rep(1,2),rep(2,25)))
    group = factor(c(rep(1,1),rep(2,3)))
    filter=0.5 
    prefix = "gor1"
    
    #for 2 x 2 experiment
    #group = factor(c(1,1,2,2))
    
    y=DGEList(counts=counts,group=group,genes=row.names(counts),remove.zeros = T)
    
    max_genes = nrow(counts)
    
    logcpm = cpm(counts,prior.count=1,log=T)
    t_cpm = cpm(counts,prior.count=1,log=F)
    logcpm = logcpm[,samples]
    t_cpm = t_cpm[,samples]
    
    plotMDS(y)
    
    keep=rowSums(cpm(y)>filter) >= n_samples/2
    y=y[keep,,keep.lib.sizes=F]
    
    #necessary for goana
    idfound = y$genes$genes %in% mappedRkeys(org.Hs.egENSEMBL)
    y = y[idfound,]
    
    egENSEMBL=toTable(org.Hs.egENSEMBL)
    m = match (y$genes$genes,egENSEMBL$ensembl_id)
    y$genes$EntrezGene = egENSEMBL$gene_id[m]
    egSYMBOL = toTable(org.Hs.egSYMBOL)
    m = match (y$genes$EntrezGene,egSYMBOL$gene_id)
    y$genes$Symbol = egSYMBOL$symbol[m]
    
    #remove duplications - just 1 gene in this dataset
    #first order by counts to remove duplicated names with 0 counts
    o = order(rowSums(y$counts),decreasing = T)
    y = y[o,]
    d = duplicated(y$genes$Symbol)
    dy = y[d,]$genes
    y = y[!d,]
    nrow(y)
    
    y$samples$lib.size = colSums(y$counts)
    rownames(y$counts) = y$genes$EntrezGene 
    rownames(y$genes) = y$genes$EntrezGene
    y$genes$EntrezGene = NULL
    
    #normalization for RNA composition (2.7.3)
    y=calcNormFactors(y)
    
    #nc=cpm(y,normalized.lib.sizes=F)
    #write.table(nc,"filtered.normalized_counts.txt",col.names=NA)
    
    png(paste0(prefix,".mds.png"),res = 300,width=2000,height=2000)
    plotMDS(y,las=1)
    dev.off()
    
    design=model.matrix(~group)
    
    y=estimateDisp(y,design)
    
    fit=glmFit(y,design)
    lrt=glmLRT(fit)
    
    efilename=paste0(prefix,".de_genes.txt")
    de_results = topTags(lrt,p.value=0.05,n=max_genes,sort.by="logFC")
    write.table(de_results,efilename,quote=F,row.names=F)
    
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

# table S10 - what outliers in muscle we can detect with myotubes
expression.outliers.trio_analysis <- function(){
    trios = paste0(c("12.1","14.1","14.2","17.1","18.1","26.1","28.1","5.1","6.1","9.1"),".Myo")
    
    outliers.muscle = read.csv("rpkms.muscle.expression.outliers.txt")
    # we have two muscle samples in case12
    outliers.muscle$Sample = gsub("Mpv","M",outliers.muscle$Sample)
    outliers.muscle$Sample = gsub("Mvl","M",outliers.muscle$Sample)
    outliers.muscle$Sample = gsub("S.._","",outliers.muscle$Sample)
    outliers.muscle$Sample = gsub("X","",outliers.muscle$Sample)
    outliers.muscle$Sample = gsub("M","Myo",outliers.muscle$Sample)
    outliers.muscle = outliers.muscle[outliers.muscle$Sample %in% trios,]
    
    outliers.muscle.number = nrow(unique(subset(outliers.muscle,select=c("Sample","Gene","Regulation"))))
    
    print(paste0(outliers.muscle.number," expression outliers detected in 10 muscle samples"))
    
    outliers.myotubes = read.csv("rpkms.myo.expression.outliers.txt")
    outliers.myotubes$Sample = gsub("S.._","",outliers.myotubes$Sample)
    outliers.myotubes = outliers.myotubes[outliers.myotubes$Sample %in% trios,]
    outliers.myotubes.number = nrow(unique(subset(outliers.myotubes,select=c("Sample","Gene","Regulation"))))
    print(paste0(outliers.myotubes.number," expression outliers detected in 10 myotubes samples"))
    
    test = merge(outliers.muscle, outliers.myotubes, by.x = c("Sample","Gene","Regulation"), by.y = c("Sample","Gene","Regulation"))
    print(paste0(nrow(test)," outlier muscle genes detected in 10 myotubes"))
    
    outliers.fibroblasts = read.csv("rpkms.fibro.expression.outliers.txt")
    outliers.fibroblasts$Sample = gsub("S.._","",outliers.fibroblasts$Sample)
    outliers.fibroblasts$Sample = gsub("F","Myo",outliers.fibroblasts$Sample)
    outliers.fibroblasts = outliers.fibroblasts[outliers.fibroblasts$Sample %in% trios,]
    test = merge(outliers.muscle, outliers.fibroblasts, by.x = c("Sample","Gene","Regulation"), by.y = c("Sample","Gene","Regulation"))
    print(paste0(nrow(test)," outlier muscle genes detected in 10 fibroblasts"))
    
    #myo vs fibroblast
    test = merge(outliers.myotubes, outliers.fibroblasts, by.x = c("Sample","Gene","Regulation"), by.y = c("Sample","Gene","Regulation"))
    print(paste0(nrow(test)," outlier myotubes genes detected in 10 fibroblasts"))
}

expression.tissue_comparison.tableS12 = function()
{
    rpkms.muscle = read.csv("rpkms.muscle.txt", sep="", stringsAsFactors=F)
    rpkms.muscle = rpkms.muscle[rpkms.muscle$external_gene_name %in% protein_coding_genes$Gene_name,]
    
    rpkms.gtex_blood = read.csv("rpkms.gtex_blood.txt", sep="", stringsAsFactors=F)
    rpkms.gtex_blood = rpkms.gtex_blood[rpkms.gtex_blood$external_gene_name %in% protein_coding_genes$Gene_name,]
    
    rpkms.muscle$external_gene_name = NULL
    rpkms.gtex_blood$external_gene_name = NULL
    
    rpkms.muscle$avg = rowMeans(rpkms.muscle)
    rpkms.gtex_blood$avg = rowMeans(rpkms.gtex_blood)
    
    rpkms.muscle = rpkms.muscle[rpkms.muscle$avg > 1,]
    
    rpkms.gtex_blood = rpkms.gtex_blood[rpkms.gtex_blood$avg > 1,]
    
    test = rpkms.gtex_blood[row.names(rpkms.gtex_blood) %in% row.names(rpkms.muscle),]
}

expression4gene_in_a_gene_panel = function(sample,gene_panel,counts)
{
    #sample = "S12_9.1.M"
    #gene_panel = congenital_myopathy
    for (gene in gene_panel)
    {
        expression4gene(gene,sample,counts)
    }
}

expression4gene = function(gene,sample,counts)
{
    #gene="LAMA2"
    #sample="S12_9.1.M"
    print(paste0(gene," ",sample))
    
    #https://datascienceplus.com/building-barplots-with-error-bars/
    gene_expression = subset(counts,external_gene_name == gene)     
    gene_expression$external_gene_name=NULL
    v_muscle=as.numeric(gene_expression[1,])
    
    se = sd(v_muscle)/sqrt(length(v_muscle))
    
    muscle_mean = mean(v_muscle)
    means = c(muscle_mean,gene_expression[[sample]])
    ses = c(se,0)
    
    #significance
    ttest  = t.test(v_muscle,mu=gene_expression[[sample]])
    
    if ((muscle_mean > 0) && (ttest$p.value < 0.01) && 
        ((muscle_mean > 1.5*gene_expression[[sample]]) || (1.5*muscle_mean < gene_expression[[sample]])))
    {
        if (gene_expression[[sample]]<muscle_mean)
        {
            file_name = paste0("_",sample,"_",gene,".png")
        }
        else
        {
            file_name = paste0(sample,"_",gene,".png")
        }
        png(file_name) 
        par(mar=c(2,3,2,0)+0.1)
        plotTop = max(mean(v_muscle)+6*se, gene_expression[[sample]]+6*se)
        barCenters = barplot(
            height = c(mean(v_muscle),gene_expression[[sample]]),
            names.arg = c("Mean muscle",sample),
            beside = F,
            las = 2,
            ylim = c(0,plotTop),
            xaxt = "n",
            axes=T,
            main = paste0(sample,",",gene," expression, RPKM"),width=c(1,1))
        
        text (x=barCenters,y=par("usr")[3]-1,
              labels = c("Mean muscle",sample),xpd=T)
        
        
        segments (barCenters,means-ses*2,barCenters,means+ses*2,lwd=1.5)
        arrows (barCenters,means-ses*2,barCenters,means+ses*2,lwd=1.5,angle=90,
                code=3,length=0.05)
        
        text(x = barCenters, y = means, label = round(means,2), 
             pos = 2, col = "red")
        
        dev.off()
    }
}

#calculate expression variability for genes downregulated in family gor1
gor1.expression_variability = function()
{
    gor1_panel = c("C2","CACNB2","CD70","CKB","CTGF","CXCL14", "EBF2", "EPB41L3", "F2R", "FABP3", 
                   "FAT3", "HES1","HSPB7", "IGFBP7", "ITGA8", "KCNE4", "MAMDC2", "NDUFA4L2","PCDHGB6", "PCSK9",
                   "PI16", "PNRC2", "PPP1R14A", "PRLR", "PRUNE2", "RHOB", "SCUBE3", "SDK2", "SPINT2", "STYK1", 
                   "TPD52L1")
    rpkm.file = "rpkms.fibro.txt"
    tissue = "fibro"
    expression_variability(gor1_panel,rpkm.file,tissue)
}

expression.variability.tableS12 = function()
{
    panel = c("DMD","NEB","ACTA1","LMNA","RYR1","MTM1","KCNIP4","FKRP","CAV3","LAMA2")
    for (tissue in c("fibro","myo","muscle"))
    {
        rpkm.file = paste0("rpkms.",tissue,".txt")
        expression_variability(panel,rpkm.file,tissue)
    }
}

# calculate a df expression variability in a gene_panel (list of genes) for 
# file like rpkms.muscle.txt - produced by crt.utils.read.coverage2counts_dir
expression.variability.get_cov = function(gene_panel,rpkm.file,tissue,debug=F)
{
    if (debug == T)
    {
        rpkm.file = "rpkms.muscle.txt"
        tissue = "Muscle"
        genes = get_genes_in_panels()
    }
    d
    rpkms = read.table(rpkm.file)
    
    agene = gene_panel[1]
    if (grepl("ENSG",agene,fixed=T)){
        rpkms = rpkms[rownames(rpkms) %in% genes,]
    }else{
        rpkms = rpkms[rpkms$external_gene_name %in% genes,]
        row.names(rpkms)=rpkms$external_gene_name
    }
    
    rpkms$external_gene_name = NULL
    
    result.table = data.frame(row.names = row.names(rpkms))
    result.table$mean = rowMeans(rpkms)
    result.table$sd = rowSds(rpkms)
    result.table$cov = with(result.table,100*sd/mean)
    
    result.table$mean = NULL
    result.table$sd = NULL
    
    colnames(result.table)[1] = tissue
    result.table = result.table[order(rownames(result.table)),,drop=F]
    
    return(result.table)
}

# boxplot of coefficient of variation (cov) for different tissue types
expression.variation.boxplot = function(genes,png.file,title)
{
    variability.muscle = expression.variability.get_cov(genes,"rpkms.muscle.txt","Muscle")
    variability.myotubes = expression.variability.get_cov(genes,"rpkms.myo.txt","Myotubes")
    variability.fibro = expression.variability.get_cov(genes,"rpkms.fibro.txt","Fibroblast")
    
    variability.gtex.muscle = expression.variability.get_cov(genes,"rpkms.50gtex.txt","Gtex muscle")
    variability.gtex.fibro = expression.variability.get_cov(genes,"rpkms.gtex_fibro.txt","Gtex fibroblast")
    variability.gtex.blood = expression.variability.get_cov(genes,"rpkms.gtex_blood.txt","Gtex_blood")
    
    variability = cbind(variability.muscle,variability.myotubes,variability.fibro,
                        variability.gtex.muscle,variability.gtex.fibro,variability.gtex.blood)
    
    png(png.file,res=300,width=3000,height=2000)
    boxplot(variability,main=title)
    dev.off()
}

expression.variation = function()
{
    genes = get_genes_in_panels()
    expression.variation.boxplot(genes,"expression.cov.panels.png",
                                 "Expression variation in muscular genes")
    
    genes = read.table("genes1rpkm_in_muscle.trios.txt",header=T,stringsAsFactors = F)$ensembl_gene_id
    expression.variation.boxplot(genes,"expression.cov.1rpkm_genes.png",
                                 "Expression variation in 1rpkm genes")
}

expression_in_a_sample <- function(){
    gene <- "NEB"
    rpkms <- read_csv("rpkms.csv") %>% filter (external_gene_name == gene)
    rpkms$ensembl_gene_id <- NULL
    rpkms$external_gene_name <- NULL
    rpkms$gene_description <- NULL
    rowMeans(rpkms)
    sd(rpkms)
}

###############################################################################
# 4. Coverage
###############################################################################
# plots coverage for every gene for all samples having sample.coverage in the current dir - output of
# bam.coverage.bamstats05.sh -v bam=file.bam,bed=~/crt/data/muscular_genes_exons.bed
# I tried to plot a panel 5x3 - too much information, some panels are too wide (many genes)
supplementary_figure4.get_coverage = function()
{
    files = list.files(".","\\.coverage$")
    #samples = unlist(read.table("samples.txt", stringsAsFactors=F))
    
    coverage = read.delim(files[1],header=T,stringsAsFactors = F)
    coverage = coverage[,c("gene","avg")]
    colnames(coverage)[2]=files[1]
    
    for (file in tail(files,-1))
    {
        sample_coverage = read.delim(file,header=T,stringsAsFactors = F)
        sample_coverage = sample_coverage[,c("gene","avg")]
        colnames(sample_coverage)[2]=file
        coverage = cbind(coverage,sample_coverage[2])
    }
    row.names(coverage) = coverage$gene
    coverage$gene=NULL
    return(coverage)
}
supplementary_figure4.plot_coverage_bins = function (coverage)
{
    n_genes = nrow(coverage)
    
    bin1 = c()
    bin2 = c()
    bin3 = c()
    bin4 = c()
    bin5 = c()
    
    for (i in 1:n_genes)
    {
        mean_expression = rowMeans(coverage[i,])
        gene = row.names(coverage)[i]
        #print(paste0(gene," ",mean_expression))
        if (mean_expression >= 10000){
            bin1 = unique(sort(c(bin1,gene)))
        }else if (mean_expression >= 1000){
            bin2 = unique(sort(c(bin2,gene)))
        }else if (mean_expression >= 100){
            bin3 = unique(sort(c(bin3,gene)))
        }else if (mean_expression >= 10){
            bin4 = unique(sort(c(bin4,gene)))
        }else{
            bin5 = unique(sort(c(bin5,gene)))
        }
    }
    
    for (bin in (c("bin1","bin2","bin3","bin4","bin5")))
    {
    
        #png(paste0(bin,".png"),res=500,width=length(get(bin))*150,height=1500)
        png(paste0(bin,".png"),res=500,width=5000,height=2000)
        boxplot(t(coverage[get(bin),]),
                las=2,
                cex.axis=0.8)
        dev.off()
    }
}

supplementary_figure4.coverage.panel = function()
{
    setwd("~/Desktop/work/RNAseq_muscle_coverage/muscle/")
    coverage.muscle = supplementary_figure4.get_coverage()   
    supplementary_figure4.plot_coverage_bins(coverage.muscle)
    
    setwd("~/Desktop/work/RNAseq_muscle_coverage/myo/")
    coverage.myo = supplementary_figure4.get_coverage()   
    supplementary_figure4.plot_coverage_bins(coverage.myo)
    
    setwd("~/Desktop/work/RNAseq_muscle_coverage/fibro/")
    coverage.fibro = supplementary_figure4.get_coverage()
    supplementary_figure4.plot_coverage_bins(coverage.fibro)
}

expression_panels = function()
{
    ensembl_w_description <- read.delim2("~/cre/ensembl_w_description.txt", row.names=1, stringsAsFactors=FALSE)
    counts = merge(counts,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(counts)=counts$Row.names
    counts$Row.names=NULL
    
    #rewrite
    #counts=subset(counts,Muscle5 !=0 & Fibroblast5 !=0 & Myotubes5!=0 & X10.1.M.bam!=0 &
    #                X11.1.K.bam != 0 & X9.1.Myo.bam !=0 & X8.1.M.bam != 0)
    counts$Gene_description = NULL
    
    counts = merge(counts,gene_lengths,by.x="row.names",by.y="row.names")
    row.names(counts)=counts$Row.names
    counts$Row.names=NULL
    
    write.table(counts,"9-1-Myo.comparison.txt",quote=F,sep = "\t")
    
    sample5 <- read.csv("9-1-Myo.comparison.txt", row.names=1, sep="", stringsAsFactors=F)
    
    x=sample5
    x$external_gene_name=NULL
    x$Length = NULL
    
    group = factor(rep(1,nrow(samples)))
    
    y=DGEList(counts=x,group=group,genes=row.names(x),remove.zeros = F)
    plotMDS(y)
    #generate counts in the other way
    rpkms = rpkm(y,sample5$Length)
    
    rpkms = merge(rpkms,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(rpkms)=rpkms$Row.names
    rpkms$Row.names=NULL
    rpkms$Gene_description=NULL
    
    write.table(rpkms,"9-1-Myo.rpkms.txt",quote=F,sep = "\t")
    
    rpkms <- read.delim("~/Desktop/project_muscular/counts/muscular_filtered/9-1-Myo.rpkms.txt", stringsAsFactors=FALSE)
    
    #MuscleGeneRPKM <- read.csv("~/Desktop/project_muscular/counts/muscular_filtered/MuscleGeneRPKM.txt", sep="", stringsAsFactors = F)
    
    plot_all_panels(rpkms)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,500,1000,2000,3000,4000)
    plot_panel(congenital_myopathy, rpkms, MuscleGeneRPKM, "congenital_myopaties",breaks)
    
    breaks = c(0,5,10,50,100,200)
    plot_panel(congenital_myastenic_syndromes, rpkms, MuscleGeneRPKM, "congenital_myastenic_syndromes",breaks)
    
    breaks = c(0,5,10,50)
    plot_panel(channelopathies, rpkms, MuscleGeneRPKM, "channelopathies",breaks)
    
    breaks = c(0,5,10,50,100,500,1000,1500)
    plot_panel(vacuolar_and_others, rpkms, MuscleGeneRPKM, "vacoular_and_others",breaks)
    
    breaks = c(0,5,10,50,100,500,1000,2000,5000,5200)
    plot_panel(distal_myopathies, rpkms, MuscleGeneRPKM, "distal_myopathies",breaks)
    
    breaks = c(0,5,10,50,100,500,1000,2000,5000,10000,15000)
    plot_panel(congenital_muscular_dystrophies, rpkms, MuscleGeneRPKM, "congenital_muscular_dystrophies",breaks)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,500,1000,1400)
    #2000,5000,5200)
    plot_panel(limb_girdle, rpkms, MuscleGeneRPKM, "Limb_girdle_dystrophies",breaks)
  
    breaks = c(0,5,10,50,100,500,1000)
    plot_panel(muscular_dystrophies, rpkms, MuscleGeneRPKM, "Muscular_dystrophies",breaks)
}

fibroblast8 = function()
{
  setwd("~/Desktop/project_muscular/counts/mh_unfiltered/")
  samples = 
  mh = read.delim("mh.txt", stringsAsFactors=F, row.names=1)
  fibroblast8 = read.delim("fibroblast8.counts", stringsAsFactors=F, row.names=1)
  
  counts = merge(mh,fibroblast8, by.x='row.names',by.y='row.names')
  row.names(counts) = counts$Row.names
  counts$Row.names = NULL
  
  congenital_muscular_dystrophies=c("LAMA2", "COL6A1", "COL6A2", "COL6A3", 
                                    "SPEN1", "FHL1", "ITGA7", "DNM2","TCAP", "LMNA", "FKTN", 
                                    "POMT1", "POMT2", "FKRP", "POMGNT1", "ISPD", "GTDC2", "B3GNT1", 
                                    "POMGNT1", "GMPPB", "LARGE", "DPM1", "DPM2", "ALG13", "B3GALNT2", 
                                    "TMEM5",  "POMK", "CHKB", "ACTA1", "TRAPPC11")
  
  
  panel.counts = counts[counts$symbol %in% gene_panel,]
  
  
  row.names(all.rpkm) = all.rpkm$external_gene_name
  all.rpkm$external_gene_name=NULL
  
  library(pheatmap)
  library(RColorBrewer)
  png(paste0(title,".png"),res=100,width=500)
  pheatmap(all.rpkm,treeheight_row=0,treeheight_col=0,cellwidth = 40,
           display_number =T,cluster_rows=T, cluster_cols=T,
           main=paste0(title," RPKM"),
           breaks=breaks,
           colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaks)-1))
  dev.off()
  
}

# use case: to plot exon coverage for muscular gene panel 
# generate coverage tables first with ~/bioscripts/bam.gene_coverage.sh
plot_exon_coverage = function(sample1,sample2,gene)
{
    barwidth = 5
    #gene="DMD"
    #sample1="9-1-M"
    #sample2="9-1-Myo"
    
    s1 = read.delim(paste0(sample1,".",gene,".coverage"), header=T, stringsAsFactors=F)
    s2 = read.delim(paste0(sample2,".",gene,".coverage"), header=T, stringsAsFactors=F)
    n_exons = nrow(s1)
    
    strand = 1
    if (n_exons > 1 && s1[1,2] > s1[2,2])
    {
         strand = -1
    }
    
    width=barwidth*n_exons*3+100
    png(paste0(gene,".png"),width=width)
    if (strand == 1)
    {
        exon_numbers = seq(1,n_exons)
    }else{
        exon_numbers = seq(n_exons,1)
    }
    df = cbind(s1[5],s2[5])
    colnames(df)=c(sample1,sample2)
    barplot(t(df),beside=T,col=c("red","blue"),names.arg = exon_numbers,width=5,xlim=c(0,3*barwidth*n_exons),
            las=2,main=gene,legend=c(sample1,sample2))
    dev.off()
}

plot_exon_coverage_main = function()
{
    setwd("/home/sergey/Desktop/project_RNAseq_diagnostics/Sample12_9-1-M_Muscle10_DMDcase/exon_coverage") 
    gene_list = read.table("gene.list", stringsAsFactors=F)
    for (gene in gene_list[[1]])
    {
        #gene="ACVR1"
        print(gene)
        plot_exon_coverage("9-1-M","9-1-Myo",gene)
    }
}

count_rpkm_for_exons = function()
{
    setwd("reference/")
    gtex.exon_reference = read.delim("gtex.exon_reference.txt", stringsAsFactors=F)
    gtex.exon_reference$length = gtex.exon_reference$stop - gtex.exon_reference$start + 1
}

###################################################################################################
# 5. Splicing 
###################################################################################################
# Supplementary_Table_9
splicing_panels_frequency = function(){
    setwd("~/Desktop/work/splicing_flank2/")
    files <- list.files(".","*rare_junctions.txt")
    events <- splicing.read_novel_splice_events(files[1])
    for (file in tail(files, -1))
    {
        events_buf <- splicing.read_novel_splice_events(file)
        events <- bind_rowd(events, events_buf)
    }
    
    events$norm_read_count <- as.numeric(events$norm_read_count)
    events$read_count <- as.numeric(events$read_count)
    
    filtered <- events[events$norm_read_count >= 0.5,]
    
    filtered <- filtered[filtered$read_count >= 30,]
    
    frequencies = as.data.frame(table(filtered$pos))
    colnames(frequencies) = c("pos","frequency")
    
    filtered = merge(filtered,frequencies,by.x = "pos", by.y = "pos",all.x = T)
    filtered = filtered[filtered$frequency <= 7,]
    
    #filtered = filtered[filtered$dup == F,]
    
    filtered$tissue = ''
    filtered$tissue = apply(filtered[,c("sample","tissue")],1,decode_tissue_type)
    
    filtered$Omim = NULL
    
    write.csv(filtered,"Supplemental_table_10.Splicing_all_genes.csv",quote=T, row.names = F)
    
    panel_genes = get_genes_in_panels()
    events = events[events$gene %in% panel_genes,]
    
    frequencies = as.data.frame(table(events$pos))
    colnames(frequencies) = c("pos","frequency")
    events = merge(events,frequencies,by.x = "pos", by.y = "pos",all.x = T)
    
    events$tissue = ''
    events$tissue = apply(events[,c("sample","tissue")],1,decode_tissue_type)
    
    #events$dup = c(duplicated(events$pos,fromLast=T) | duplicated(events$pos))
    
    eoutliers <- read.csv("~/Desktop/work/expression/outliers_panels/outliers.txt", stringsAsFactors=F)
    eoutliers = subset(eoutliers,select=c("Sample","Gene","Regulation","Abs_FC_cohort","Abs_FC_GTex"))
    eoutliers$Sample = gsub("[.]","-",eoutliers$Sample)
    eoutliers = unique(eoutliers)
    
    res = merge(events,eoutliers,by.x = c("sample","gene"),by.y = c("Sample","Gene"),all.x=T,all.y=F)
    
    res$Omim = NULL
    
    write.csv(res,"TableS9.Splicing.Panels.Frequency.csv",quote=T,row.names = F)
}

TableS10.S11.S12.Splicing.Statistics = function()
{
    setwd("~/Desktop/work/splicing_flank2/")
    junctions <- read.csv("TableS9.Splicing.Panels.Frequency.csv")

    junctions.muscle = junctions[junctions$tissue =='Muscle',]
    junctions.muscle.pos = junctions.muscle$pos
    junctions.unique = unique(junctions.muscle[,c("gene","pos")])
    
    print(paste0(nrow(junctions.unique)," unique junctions discovered in muscle samples"))
    
    junctions.muscle.frequency = count(junctions.muscle.pos)
    junctions.muscle.frequency = merge(junctions.muscle.frequency,junctions.unique,by.x="x",by.y="pos",all.x=T,all.y=F)
    colnames(junctions.muscle.frequency) = c("Junction","Frequency","Gene")
    
    junctions.muscle.type = unique(junctions.muscle[,c("pos","annotation")])
    colnames(junctions.muscle.type) = c("pos","Annotation")
    junctions.muscle.frequency = junctions.muscle.frequency[,c("Gene","Junction","Frequency")]
    junctions.muscle.frequency = merge(junctions.muscle.frequency,junctions.muscle.type,by.x="Junction",by.y="pos",all.x=T,all.y=F)
    junctions.muscle.frequency = junctions.muscle.frequency[order(-junctions.muscle.frequency$Frequency),]
    write.csv(junctions.muscle.frequency,"TableS10.Splicing.Statistics.By_junction.csv",quote=T,row.names =F)
    
    frequency.by_sample = count(junctions.muscle$sample)
    colnames(frequency.by_sample)=c("Sample","Novel_junctions")
    frequency.by_sample = frequency.by_sample[order(-frequency.by_sample$Novel_junctions),]
    write.csv(frequency.by_sample,"TableS11.Splicing.Novel_junctions.By_Sample.csv",quote=T,row.names =F)
    
    frequency.by_type = count(junctions.muscle$annotation)
    colnames(frequency.by_type)=c("Even_type","Frequency")
    frequency.by_type = frequency.by_type[order(-frequency.by_type$Frequency),]
    write.csv(frequency.by_type,"TableS12.Splicing.Novel_junctions.By_type.csv",quote=T,row.names =F)
}

TableS13 = function()
{
    samples <- read.csv("TableS1.Samples.csv")
    samples = samples[samples$Trio == 'y',]
    samples$full_name = paste0(samples$Bioinf_sample_id,"_",samples$Sample_name)
    
    junctions.trios = junctions[junctions$sample %in% samples$full_name,]
    
    junctions.muscle = junctions.trios[junctions.trios$tissue == "Muscle",]
    frequencies.muscle = as.data.frame(table(junctions.muscle$pos))
    colnames(frequencies.muscle) = c("pos","frequency.muscle")
    junctions.trios = merge(junctions.muscle,frequencies.muscle,by.x = "pos", by.y = "pos",all.x = T)
    
    
    junctions.myo = junctions.trios[junctions.trios$tissue == "Myotubes",]
    frequencies.myo = as.data.frame(table(junctions.myo$pos))
    colnames(frequencies.myo) = c("pos","frequency.myo")
    junctions.trios = merge(junctions.trios,frequencies.myo,by.x = "pos", by.y = "pos",all.x = T)
    
    res$unique_in_myo = ifelse(res$frequency > res$frequency.myo, 'N', 'Y')
}

TableS15.expression.1rpkm <- function(rpkms.file){
    genes=get_genes_in_panels()
    rpkms.muscle <- read.csv(rpkms.file, sep="", stringsAsFactors=FALSE)
    rpkms.muscle = rpkms.muscle[rpkms.muscle$external_gene_name %in% genes,]
    row.names(rpkms.muscle) = rpkms.muscle$external_gene_name
    rpkms.muscle$external_gene_name = NULL
    rpkms.muscle$average = rowMeans(rpkms.muscle)
    print(paste0("Genes at >=1 RPKM:",length(row.names(rpkms.muscle[rpkms.muscle$average>=1,]))))
    print(paste(sort(row.names(rpkms.muscle[rpkms.muscle$average>=1,])),collapse=","))
    
    print(paste0("Genes at <1 RPKM:",length(row.names(rpkms.muscle[rpkms.muscle$average<1,]))))
    print(paste(sort(row.names(rpkms.muscle[rpkms.muscle$average<1,])),collapse=","))
}

###############################################################################
args <- commandArgs(trailingOnly = T)
if (length(args) == 0 || args[1] == "--help"){
    cat("Usage: Rscript function_name function_args\n")
    cat("Available functions:\n")
    cat("splicing_panels_frequency\n")
}else{
    cat(paste0("Running function: ", args[1],"\n"))
    init()
    fcn <- get(args[1])
    fcn(tail(args,-1))
}
###############################################################################