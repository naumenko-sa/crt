###################################################################################################
###  Very basic allelic imbalance analysis in the spirit of Cummings2017
###  for more sophisticated see @Pejminister and @tuuliel
###################################################################################################
imbalance.init()
{
    library("genefilter")
    source("~/crt/crt.utils.R")
    
}

# calculate median allelic imbalance ratio for protein coding genes in a sample
# input: 
# - sample = sample_name, i.e. 10-1-M, het file should be sample.het.csv in the wd, output of ~/crt/crt.imbalance.get_hets.sh
# protein_coding_genes.bed is from the global namespace, recycled each time
# imbalance statistics ratio=min(ref,alt)/max(ref,alt), ideal het ratio = 0.5, ideal hom ratio=0
# output:
# - sample.ai.csv, format: gene,sample (gene_name,imbalance %)
imbalance.variants2imbalance = function(sample)
{
    # debug:
    # setwd("~/Desktop/work/imbalance/")
    # sample = "10-1-M"
    infile = paste0(sample,".het.csv")
    outfile = paste0(sample,".ai.csv")
    sample.het = read.csv(infile, stringsAsFactors=F)
    sample.het$ratio = with(sample.het,pmin(ref,alt)/pmax(ref,alt))
    
    protein_coding_genes.bed$allele_ratio = NULL
    
    for (i in (1:nrow(protein_coding_genes.bed)))
    {
        gene_chrom = protein_coding_genes.bed[i,"chrom"]
        gene_start = protein_coding_genes.bed[i,"start"]
        gene_end = protein_coding_genes.bed[i,"end"]
        gene_variants = subset(sample.het, chrom == gene_chrom & pos >= gene_start & pos <= gene_end)
        
        if (nrow(gene_variants)>=5){
            protein_coding_genes.bed[i,"allele_ratio"] = round(median(gene_variants$ratio),2)
        }else{
            protein_coding_genes.bed[i,"allele_ratio"] = NA
        }
    }
    #protein_coding_genes.bed = protein_coding_genes.bed [!is.na(protein_coding_genes.bed$ai),]
    protein_coding_genes.bed = protein_coding_genes.bed[,c("ensembl_gene_id","gene","allele_ratio")]
    colnames(protein_coding_genes.bed) = c("ensembl_gene_id","gene",sample)
    write.csv(protein_coding_genes.bed,outfile,row.names = F)
}

imbalance.gtex_reference()
{
    setwd("~/Desktop/work/imbalance/gtex")
    
    files = list.files(".","*ai.csv")
    imbalance = read.csv(files[1],stringsAsFactors = F)
    genes_names = imbalance
    genes_names$GTEX.111CU.2026=NULL
    
    
    row.names(imbalance)=imbalance$ensembl_gene_id
    imbalance$gene=NULL
    imbalance$ensembl_gene_id=NULL
    
    for (file in tail(files,-1))
    {
        print(file)
        imbalance_buf = read.csv(file, stringsAsFactors = F)
        row.names(imbalance_buf)=imbalance_buf$ensembl_gene_id
        
        imbalance_buf$gene=NULL
        imbalance_buf$ensembl_gene_id=NULL
        imbalance = cbind(imbalance,imbalance_buf)
    }
    imbalance.stats = data.frame(row.names = row.names(imbalance))
    imbalance.stats$mean = rowMeans(imbalance,na.rm = T)
    imbalance.stats$sd = rowSds(imbalance,na.rm = T)
    imbalance.stats = imbalance.stats[!is.nan(imbalance.stats$mean),]
    
    imbalance.stats = merge(imbalance.stats,genes_names,by.x="row.names",by.y="ensembl_gene_id",all.x=T,all.y=F)
    colnames(imbalance.stats)[1]="ensembl_gene_id"
    
    write.csv(imbalance.stats,"gtex.imbalance.csv",row.names=F)
    
}

imbalance.zscore = function(sample)
{
    setwd("~/Desktop/work/imbalance/")
    gtex.imbalance = read.csv("gtex.imbalance.csv",stringsAsFactors = F)
    gtex.imbalance$gene=NULL
    
    filename = paste0(sample,".ai.csv")
    buf = read.csv(filename, stringsAsFactors = F,header = T)
    colnames(buf)[3]='imbalance'
    buf = buf[!is.na(buf$imbalance),]
    
    buf = merge(buf,gtex.imbalance,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.x=T,all.y=F)
    buf$zscore = (buf$imbalance - buf$mean)/buf$sd
    buf$sample = sample
    
    buf=buf[,c("sample","ensembl_gene_id","gene","zscore","imbalance","mean","sd")]  
    write.csv(buf,paste0(sample,".imbalance.csv"),row.names=F)   
}

supplementary_table.Allelic_imbalance.z-scores.cvs()
{
    samples = read.table("samples.txt",stringsAsFactors = F)
    for (sample in samples$V1)
    {
        print(sample)
        imbalance.zscore(sample)
    }
    supplemental_table = read.csv(paste0(samples$V1[1],".imbalance.csv"),stringsAsFactors = F)
    for (sample in tail(samples$V1,-1)){
        buf = read.csv(paste0(sample,".imbalance.csv"),stringsAsFactors = F)
        supplemental_table = rbind(supplemental_table,buf)
    }
    write.csv(supplemental_table,"Supplemental_table21.Allelic_imbalance.Z-scores.csv",row.names = F)
}

