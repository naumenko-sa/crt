###############################################################################
# Expression outlier detection based on Z-scores, detects control outliers well
###############################################################################
expression_outliers_init <- function(){
    library(genefilter)
    source("~/crt/scripts/crt.utils.R")
    source("~/bioscripts/gene_panels/genes.R")
}

###############################################################################
# Supplemental tables for Gonorazky.Naumenko.2019 RNA-seq article
###############################################################################
# returns 4 tables:
# - for muscle genes
# - for omim genes
# - for mitochondrial genes in case40
# - for all genes
get_expression_outliers_muscle <- function(){
    setwd("~/Desktop/work/expression")
    
    muscular_genes <- read_csv("~/bioscripts/gene_panels/gonorazky2019_muscular_genes.csv")
    get_expression_outliers_zscore("rpkms.muscle.csv",
                               "rpkms.50gtex.csv",
                               muscular_genes,
                              "Supplemental_table_9.Expression_outliers_in_muscular_panels.csv")
    
    omim.genes <- read.csv("~/cre/data/omim.genes.csv")
    get_expression_outliers_zscore("rpkms.muscle.txt",
                               "rpkms.50gtex.txt",
                               omim.genes$ensembl_gene_id,
                               "Supplemental_table_10.Expression_outliers_in_OMIM.csv")
    
    mitocarta.genes <- read.csv("~/cre/data/mitocarta.ensembl_ids.txt")
    get_expression_outliers_zscore("rpkms.muscle.txt",
                               "rpkms.50gtex.txt",
                               mitocarta.genes$ensembl_gene_id,
                               "Supplemental_table_11.Expression_outliers.case40.csv")
    
    
    case40.outliers <- read.csv("Supplemental_table_11.Expression_outliers.case40.csv",
                                stringsAsFactors = F)
    case40.outliers <- subset(case40.outliers, sample == "S70_40.1.M" | sample == "S71_40.2.M")
    write.csv(case40.outliers, "Supplemental_table_11.Expression_outliers.case40.csv", 
              row.names = F)
}

# comparing fibroblasts to fibroblasts because GTex fibroblasts have different expression profile
get_expression_outliers_fibro <- function(){
    setwd("~/Desktop/work/expression/feature_counts/fibro/")
    #read_feature_counts_dir()
    omim_genes <- read.csv("~/cre/data/omim.genes.csv")
    get_expression_outliers_zscore("rpkms.fibro.txt",
                               "rpkms.fibro.txt",
                               omim_genes$Ensembl_gene_id,
                               "Supplemental_table_X.Expression_outliers_in_fibroblasts.csv")
}

get_expression_outliers_myotubes <- function(){
    setwd("~/Desktop/work/expression/feature_counts/myotubes/")
    #read_feature_counts_dir()
    muscular_genes <- read.csv("~/cre/data/muscular_genes.csv", stringsAsFactors = F)
    get_expression_outliers_zscore("rpkms.myo.txt",
                                   "rpkms.myo.txt",
                                   muscular_genes$ensembl_gene_id,
                                   "Supplemental_table_X.Expression_outliers_in_myotubes.csv")
}

# expression outliers detection using Z-scores:
# mean GTEx rpkm (sample rpkm might be lower) >= 0.1 (1 rpkm misses some cases)
# abs(z-score) >= 1.5
# input:
# - rpkms.patients.filename - table of RPKM values for all genes and all patients, output of ~/crt/crt.utils.R/read.rpkm_counts_dir
# - rpkms.gtex.filename - table of RPKM values for all genes for gtex controls, output of ~/crt/crt.utils.R/read.rpkm_counts_dir
# - gene_panel - vector of ENSEMBL_IDs of genes to filter (muscular genes, omim, mitochondrial, protein_coding)
# - output.file.name = name of csv file to write expression outliers table
get_expression_outliers_zscore <- function(rpkms_patients_filename,
                                      rpkms_gtex_filename,
                                      gene_panel,
                                      output_file_name = "expression_outliers.csv")
{
    # debug
    # setwd("~/Desktop/work/expression")
    # rpkms_patients_filename = "rpkms.muscle.csv"
    # rpkms_gtex_filename = "rpkms.50gtex.csv"
    expression_outliers_init()
    
    rpkms_gtex <- read_csv(rpkms_gtex_filename) %>%
        semi_join(protein_coding_genes, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
        mutate(gtex_mean = rowMeans(select(., starts_with("GTEX"))),
               gtex_sd = rowSds(select(., starts_with("GTEX")))) %>% 
        select(ensembl_gene_id, external_gene_name, gene_description, gtex_mean, gtex_sd)
    
    rpkms_patients <- read_csv(rpkms_patients_filename) %>%
        semi_join(protein_coding_genes, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
        select(-external_gene_name, -gene_description)
        
    expression_stats <- rpkms_patients %>% mutate(cohort_mean = rowMeans(select(., -ensembl_gene_id)),
               cohort_sd = rowSds(select(., -ensembl_gene_id))) %>% 
        select(ensembl_gene_id, cohort_mean, cohort_sd) %>% 
        left_join(rpkms_gtex, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
        filter(gtex_mean >= 0.1)
    
    rpkms_patients <- rpkms_patients %>% gather(sample, expression, -ensembl_gene_id) %>% 
        inner_join(expression_stats, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
        mutate(zscore_gtex = (expression-gtex_mean)/gtex_sd,
               zscore_cohort = (expression-cohort_mean)/cohort_sd,
               fold_change = ifelse(expression>=gtex_mean, expression/gtex_mean, gtex_mean/expression),
               regulation = ifelse(zscore_gtex > 0, "UP", "DOWN"),
               comment = "") %>% 
        filter(abs(zscore_gtex) >= 1.5) %>% 
        select(one_of("comment", "sample", "ensembl_gene_id", "external_gene_name", "expression",
                                                 "gtex_mean", "cohort_mean", "fold_change", "regulation",
                                                 "zscore_gtex", "zscore_cohort","gtex_sd","cohort_sd")) %>% 
        semi_join(gene_panel, by = c("ensembl_gene_id" = "ensembl_gene_id"))
    write_excel_csv(rpkms_patients, output_file_name)
}

###################################################################################################
### development and testing
###################################################################################################

###################################################################################################
### OUTRIDER 
###################################################################################################
expression.outliers.outrider.install = function()
{
    # using outrider: https://github.com/gagneurlab/OUTRIDER
    # https://github.com/gagneurlab/OUTRIDER/issues/8
    # https://bioconductor.org/packages/devel/bioc/html/OUTRIDER.html
    # https://bioconductor.org/packages/devel/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf
    
    install.packages("devtools")
    source("https://bioconductor.org/biocLite.R")
    biocLite("BiocInstaller")
    
    #also works to update
    devtools::install_github("gagneurlab/OUTRIDER",dependencies=TRUE)
}

Supplemental_table9.expression.outliers.OUTRIDER = function()
{
    library("OUTRIDER")
    setwd("~/Desktop/work")
    patients_counts = read.table("muscle.raw_counts.txt")
    patients_counts = patients_counts[row.names(patients_counts) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    
    gtex_counts = read.table("muscle.gtex.raw_counts.txt")
    gtex_counts = gtex_counts[row.names(gtex_counts) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    
    ensembl_ids_names = read.csv("~/cre/data/genes.transcripts.csv")
    ensembl_ids_names$Ensembl_transcript_id=NULL
    ensembl_ids_names = unique(ensembl_ids_names)
    ensembl_ids_names = ensembl_ids_names[ensembl_ids_names$Ensembl_gene_id %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    
    genes = get_genes_in_panels()
    omim.genes = read.csv("~/cre/data/omim.genes.csv")
    
    outliers.panels = data.frame()
    outliers.omim = data.frame()
    
    for (sample in colnames(patient_counts))
    {
        #debug:        sample="X10.1.M"
        #debug: 
        sample="S12_9.1.M"
        print(sample)
        patient_count = subset(patients_counts,select=sample)
        counts = cbind(gtex_counts, patient_count)
        
        ods <- OutriderDataSet(countData = counts)
        
        # 1 option) changing to minCounts (only filter genes with less then 1 read over all samples)
        ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
        
        # 2 option) Use annotation to filter by FPKM values (basepair column are needed for that, as stated in the error)
        #library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        #library(org.Hs.eg.db)
        #txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
        #map <- select(org.Hs.eg.db, keys=keys(txdb, keytype = "GENEID"),
        #              keytype="ENTREZID", columns=c("SYMBOL"))
        #ods <- filterExpression(ods, filterGenes=TRUE, mapping=map, gtf=txdb)
        
        ods <- OUTRIDER(ods)
        
        res = results(ods, all=T)   
        
        #res.pvalue = res[order(abs(l2fc)), p_rank:=1:.N, by=sampleID]
        #res.pvalue = res.pvalue[p_rank <= 10 & sampleID == 'patient_counts$S12_9.1.M']
        
        #rank by abc(lfc)
        #res.lfc = res[res$sampleID=="X10.1.M",]
        #res.lfc = merge(res.lfc,ensembl_ids_names,by.x="geneID",by.y="Ensembl_gene_id",all.x=T)
        #res.lfc = res.lfc[order(abs(l2fc),decreasing=T)]
        #write.csv(res.lfc,"10-1-M.table.csv",row.names = F)
        #res.pvalue = res.pvalue[p_rank <= 10 & sampleID == 'patient_counts$S12_9.1.M']
        
        # Rank by Z score
        res.zscore = res[order(abs(zScore), decreasing=TRUE), z_rank:=1:.N, by=sampleID]
        res.zscore = res.zscore[sampleID == sample,]
        
        res.zscore = merge(res.zscore,ensembl_ids_names,by.x='geneID',by.y='Ensembl_gene_id',all.x=T,all.y=F)
        
        res.zscore.mgenes = res.zscore[res.zscore$external_gene_name %in% genes,]
        res.zscore.mgenes = res.zscore.mgenes[abs(res.zscore.mgenes$zScore)>=2]
        res.zscore.mgenes = res.zscore.mgenes[order(abs(res.zscore.mgenes$zScore),decreasing = T)]
        
        outliers.panels = rbind(outliers.panels,res.zscore.mgenes)
        
        res.zscore.omim = res.zscore[res.zscore$geneID %in% omim.genes$Ensembl_gene_id,]
        res.zscore.omim = res.zscore.omim[abs(res.zscore.omim$zScore)>=2]
        res.zscore.omim = res.zscore.omim[order(abs(res.zscore.omim$zScore),decreasing = T)]
        
        outliers.omim = rbind(outliers.omim,res.zscore.omim)
    }
    
    outliers.panels = outliers.panels[,c("sampleID","external_gene_name","geneID","zScore","l2fc","rawcounts","normcounts","mu","disp","meanCorrected")]
    write.csv(outliers.panels,"outliers.panels.csv",row.names=F)
    outliers.omim = outliers.omim [,c("sampleID","external_gene_name","geneID","zScore","l2fc","rawcounts","normcounts","mu","disp","meanCorrected")]
    write.csv(outliers.omim,"outliers.omim.csv",row.names=F)
    
}

expression.outliers.outrider.test = function()
{
    pres = res[res$sampleID=="patient_counts$S12_9.1.M",]
    pres = pres[pres$zScore < -2,]
    
    View(res)
    
    dim(res)
    
    head(res,10)
    
    ods = OutriderDataSet(countData=counts)
    
    plotAberrantPerSample(ods,padjCutoff=0.05)
    
    
    URL <- paste0("https://media.nature.com/original/nature-assets/",
                  "ncomms/2017/170612/ncomms15824/extref/ncomms15824-s1.txt")
    ctsTable <- read.table(URL, sep="\t")
    # create OutriderDataSet object
    ods <- OutriderDataSet(countData=ctsTable)
    
    # get annotation
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    map <- select(org.Hs.eg.db, keys=keys(txdb, keytype = "GENEID"),
                  keytype="ENTREZID", columns=c("SYMBOL"))
    
    
    ods <- filterExpression(ods, txdb, mapping=map,
                            filterGenes=FALSE, savefpkm=TRUE)  
    
    plotFPKM(ods)
    
    ods <- ods[mcols(ods)$passedFilter,]
    
    ods <- estimateSizeFactors(ods)
    ods <- autoCorrect(ods, q=13)
    
    ods <- plotCountCorHeatmap(ods, normalized=FALSE, nCluster=4)
    
    ods = fit(ods)
    plotDispEsts(ods)
    
    ods <- computePvalues(ods, alternative="two.sided", method="BY")
    
    ods <- computeZscores(ods)
    
    res = results(ods,all=T )   
    
    
    gene="ENSG00000198947"  
    plotVolcano(ods, gene, basePlot=TRUE)
    plotExpressionRank(ods, gene, basePlot = F)
    plotQQ(ods,res[1,])
}
###################################################################################################
###  Naive t-test
###################################################################################################
# output expression outliers as a table eoutliers.txt
# this test is too sensitive to apply for all protein coding genes
# decided to compare with GTEx rpkms not with cohort,
# fold change reported vs cohort and gtex
Supplemental_Table9.Expression.Outliers.Panels.Naive.Ttest = function(file.rpkms, for_panels = T)
{
    # DEBUG
    for_panels = T
    setwd("~/Desktop/work/expression")
    file.rpkms = "rpkms.muscle.txt"
    # DEBUG
    
    gtex.rpkms = read.table("rpkms.50gtex.txt")  
    output = gsub("txt","expression.outliers.txt",file.rpkms)
    cat("Sample,Gene_panel_name,Gene,Regulation,Abs_FC_cohort,Abs_FC_GTex,Pvalue",file=output,append=F,sep="\n")
    
    counts = read.table(file.rpkms)
    
    for (wrong_gene in c("ENSG00000261832","ENSG00000264813","ENSG00000258529","ENSG00000273170"))
    {
        counts = counts[setdiff(rownames(counts),wrong_gene),]
    }
    
    counts = counts[counts$external_gene_name %in% protein_coding_genes$Gene_name,]
    samples = head(colnames(counts),-1)
    
    if (for_panels)
    {
        gene_list = c()
        for (gene_panel_name in panel_list)
        {
            gene_panel = get(gene_panel_name)
            gene_list = unique(c(gene_list,gene_panel))
        }
        gene_panel_name = "muscular_genes"
        fc_threshold = 2
    }else{
        gene_panel_name = "protein_coding_genes"
        gene_list = protein_coding_genes$Gene_name
        fc_threshold = 50
    }
    
    for (sample in samples)
    {
        for (gene in gene_list)
        {
            expression4gene.table(gene,sample,counts,gene_panel_name,output,fc_threshold,gtex.rpkms)
        }
    }
}

expression4gene.table = function(gene,sample,counts,gene_panel_name,output,fc_threshold,gtex.rpkms)
{
    #gene = "DMD"
    #sample = "S12_9.1.M"
    
    gene_expression = counts[counts$external_gene_name == gene,]
    gtex_gene_expression = gtex.rpkms[gtex.rpkms$external_gene_name == gene,]
    
    if (nrow(gene_expression)>1)
    {
        gene_expression = head(gene_expression,1)
    }
    
    if (nrow(gtex_gene_expression)>1)
    {
        gtex_gene_expression = head(gtex_gene_expression,1)
    }
    
    if (nrow(gene_expression)>0)
    {
        gene_expression$external_gene_name = NULL
        gtex_gene_expression$external_gene_name = NULL
        
        v_cohort = as.numeric(gene_expression[1,])
        v_gtex = as.numeric(gtex_gene_expression[1,])
        
        v_cohort.mean = mean(v_cohort)
        v_gtex.mean = mean(v_gtex)
        
        muscle_mean.gtex = as.numeric(gtex_rpkm[gtex_rpkm$gene_name %in% c(gene),]$GTEX)
        
        expression_sample = gene_expression[[sample]]
        #significance
        ttest = t.test(v_gtex,mu=expression_sample)
        
        if ((v_cohort.mean > 0) && (expression_sample > 0)){
            fold_change.cohort = (max(v_cohort.mean, expression_sample)/min(v_cohort.mean,expression_sample))
            fold_change.gtex = (max(v_gtex.mean, expression_sample)/min(v_gtex.mean,expression_sample))
            
            if (expression_sample > v_gtex.mean){
                regulation = "UP"
            }else{
                regulation = "DOWN"
            }
            
            #if ((fold_change.cohort > 1.5 || fold_change.gtex > 1.5) && (ttest$p.value < 0.01)){
            if (fold_change.cohort > fc_threshold && (ttest$p.value < 0.01)){
                cat(paste(sample,gene_panel_name,gene,regulation,fold_change.cohort,fold_change.gtex,ttest$p.value,sep = ","),file = output,append=T,sep="\n")   
            }
        }
    }
}
