################################################################################################################
# Current production method based just on Z-scores, detects control outliers well,
# used in Gonorazky.Naumenko.2018 RNA-seq article
################################################################################################################

# expression outliers detection using Z-scores:
# rpkm > 0.1 filter (possible 1 rpkm misses some cases)
# abs(z-score) > 1.5
# returns 3 tables:
# - for muscle genes
# - for omim genes
# - for mitochondrial genes in case40
expression.outliers.zscore = function()
{
    setwd("~/Desktop/work/expression")
    rpkms.50gtex = read.table("rpkms.50gtex.txt",stringsAsFactors = F)
    rpkms.patients = read.table("rpkms.muscle.txt",stringsAsFactors = F)
    
    rpkms.50gtex = rpkms.50gtex[row.names(rpkms.50gtex) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    gene_names = subset(rpkms.50gtex,select=c("external_gene_name"))
    rpkms.50gtex$external_gene_name = NULL
    
    rpkms.patients = rpkms.patients[row.names(rpkms.patients) %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
    rpkms.patients$external_gene_name = NULL
    
    expression.stats = data.frame(row.names = row.names(rpkms.50gtex))
    expression.stats$gtex_mean = rowMeans(rpkms.50gtex)
    expression.stats$gtex_sd = rowSds(rpkms.50gtex)
    expression.stats$cohort_mean = rowMeans(rpkms.patients)
    expression.stats$cohort_sd = rowSds(rpkms.patients)
    
    expression.stats = expression.stats[expression.stats$gtex_mean >=0.1,]
    
    rpkms.patients = rpkms.patients[row.names(rpkms.patients) %in% row.names(expression.stats),]
    
    expression.outliers = (rpkms.patients - expression.stats$gtex_mean)/expression.stats$gtex_sd
    
    measure_vars = colnames(expression.outliers)
    expression.outliers$ensembl_gene_id = row.names(expression.outliers)
    
    expression.outliers = melt(expression.outliers,
                               measure.vars = measure_vars,
                               variable.name = "sample",value.name = "zscore_gtex")
    expression.outliers = expression.outliers[abs(expression.outliers$zscore)>=1.5,]
    expression.outliers$gene = gene_names[expression.outliers$ensembl_gene_id,1]
    
    zscore.cohort = (rpkms.patients - expression.stats$cohort_mean)/expression.stats$cohort_sd
    measure_vars = colnames(zscore.cohort)
    zscore.cohort$ensembl_gene_id = row.names(zscore.cohort)
    zscore.cohort = melt(zscore.cohort,measure.vars = measure_vars,variable.name = "sample",value.name = "zscore_cohort")
    
    expression.outliers = merge(expression.outliers,zscore.cohort,
                                by.x=c("ensembl_gene_id","sample"),
                                by.y=c("ensembl_gene_id","sample"),
                                all.x=T,all.Y=F)
    
    expression.outliers$gtex_mean_rpkm = expression.stats[expression.outliers$ensembl_gene_id,"gtex_mean"]
    expression.outliers$gtex_sd = expression.stats[expression.outliers$ensembl_gene_id,"gtex_sd"]
    expression.outliers$cohort_mean_rpkm = expression.stats[expression.outliers$ensembl_gene_id,"cohort_mean"]
    expression.outliers$cohort_sd = expression.stats[expression.outliers$ensembl_gene_id,"cohort_sd"]
    
    measure_vars = colnames(rpkms.patients)
    rpkms.patients$ensembl_gene_id = row.names(rpkms.patients)
    rpkms.patients = melt(rpkms.patients,measure.vars = measure_vars,
                          variable.name="sample",
                          value.name="expression_rpkm")
    
    expression.outliers = merge(expression.outliers,rpkms.patients,by.x=c("ensembl_gene_id","sample"),
                                by.y=c("ensembl_gene_id","sample"),all.x=T,all.y=F)
    
    expression.outliers$fold_change = ifelse(expression.outliers$expression_rpkm>=expression.outliers$gtex_mean_rpkm,
                                             expression.outliers$expression_rpkm/expression.outliers$gtex_mean_rpkm,
                                             expression.outliers$gtex_mean_rpkm/expression.outliers$expression_rpkm)
    expression.outliers$regulation = ifelse(expression.outliers$zscore_gtex>0,"UP","DOWN")
    expression.outliers$comment = ""
    
    expression.outliers = expression.outliers[,c("comment","sample","ensembl_gene_id","gene","expression_rpkm",
                                                 "gtex_mean_rpkm","cohort_mean_rpkm","fold_change","regulation",
                                                 "zscore_gtex","zscore_cohort","gtex_sd","cohort_sd")]
    
    muscular_genes = get_genes_in_panels()
    
    expression.outliers.panels = expression.outliers[expression.outliers$gene %in% muscular_genes,]
    write.csv(expression.outliers.panels,
              "Supplemental_table_9.Expression_outliers_in_muscular_panels.csv",row.names = F)
    
    omim.genes = read.csv("~/cre/data/omim.genes.csv")
    expression.outliers.omim = expression.outliers[expression.outliers$ensembl_gene_id %in% omim.genes$Ensembl_gene_id,]
    write.csv(expression.outliers.omim,
              "Supplemental_table_10.Expression_outliers_in_OMIM.csv",row.names = F)
    
    mitocarta.genes = read.csv("~/cre/data/mitocarta.ensembl_ids.txt")
    expression.outliers.case40 = subset(expression.outliers,sample=="S70_40.1.M" | sample == "S71_40.2.M")
    expression.outliers.case40 = expression.outliers.case40[expression.outliers.case40$ensembl_gene_id %in% mitocarta.genes$ensembl_gene_id,]
    write.csv(expression.outliers.omim,"Supplemental_table_11.Expression_outliers.case40.csv",row.names = F)
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