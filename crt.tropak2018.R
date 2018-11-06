init = function()
{
    library("edgeR")
    library("genefilter")
    source("~/crt/crt.utils.R")
    setwd("~/Desktop/work/creatinine/")
}

fig1.mds = function(refresh_files = F)
{
    refresh_files=F
    print("Reading counts ...")
    counts = read.feature_counts_dir(update=refresh_files)
    #group = factor(c(rep(1,ncol(counts))))
    
    sample_names = colnames(counts)
    sample_names = strsplit2(sample_names,"_",fixed=T)[,1]
    sample_labels = sample_names
    colnames(counts) = sample_names
    n_samples = length(sample_names)
    
    group = factor(c(rep(1,n_samples/2), rep(2,n_samples/2)))    
    v_colors = c(rep("green",5),rep("red",5))
    
    print("Subsetting protein coding genes ...")
    #a file with ENSEMBL IDs
    # get mouse protein coding genes
    if (file.exists("~/Desktop/reference_tables/protein_coding_genes.list"))
    {
        protein_coding_genes <- read.csv("~/Desktop/reference_tables/protein_coding_genes.list", sep="", stringsAsFactors=FALSE)
        counts = counts[row.names(counts) %in% protein_coding_genes$ENS_GENE_ID,]
    }else{
        print("Please provide protein_coding_genes.list with ENS_GENE_ID")
    }
    
    print("Removing zeroes ...")
    
    y = DGEList(counts=counts,group=group,remove.zeros = T)
    png("mds.png",res=300,width=2000,height=2000)
    
    print("Plotting ...")
    
    mds = plotMDS(y,labels=sample_labels)
    
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
    
    png("mds.no_sample6_7.labels.png",res=300,width=2000,height=2000)
    plotMDS(y, labels=sample_labels)
    dev.off()
}


differential_expression = function()
{
    counts = read.feature_counts_dir()
    #group = factor(c(rep(1,ncol(counts))))
    
    sample_names = colnames(counts)
    sample_names = strsplit2(sample_names,"_",fixed=T)[,1]
    sample_labels = sample_names
    colnames(counts) = sample_names
    n_samples = length(sample_names)
    
    group = factor(c(rep(1,n_samples/2), rep(2,n_samples/2)))
    
    x=counts
    max_genes = nrow(counts)
    
    y=DGEList(counts=x,group=group,genes=row.names(x),remove.zeros = T)
    rpkm.counts = rpkm(y)
    
    plotMDS(y)
    #filter - 1 or 0.5
    filter=0.5
    keep=rowSums(cpm(y)>filter) >= n_samples/2
    y=y[keep,,keep.lib.sizes=F]
    
    y$samples$lib.size = colSums(y$counts)
    
    
    
    #normalization for RNA composition (2.7.3)
    y=calcNormFactors(y)
    
    #nc=cpm(y,normalized.lib.sizes=F)
    #write.table(nc,"filtered.normalized_counts.txt",col.names=NA)
    
    plotMDS(y,las=1)
    
    design=model.matrix(~group)
    
    y=estimateDisp(y,design)
    
    fit=glmFit(y,design)
    lrt=glmLRT(fit)
    
    #o=order(lrt$table$PValue)
    #cpm(y)[o[1:10],]
    #write.table(cpm(y)[o[1:12566],],"allgenes.cpm")
    
    de_results = topTags(lrt,n=max_genes,sort.by="PValue",p.value=1,adjust.method = "fdr")
    
    write.csv(topTags(lrt,p.value=0.05,n=10000),"de.csv",row.names = F)
    
    de_results = read.csv("de.csv", stringsAsFactors=FALSE)
    s_rownames = row.names(de_results)
    #setnames(de_results,"genes","ensembl_gene_id")
    #de_results = lrt$table
    
    gene_descriptions = read.csv("~/cre/data/ensembl_w_description.mouse.csv", stringsAsFactors=FALSE)
    
    de_results = merge(de_results,gene_descriptions,by.x="genes",by.y="ensembl_gene_id",all.x=T)
    #de_results = rename(de_results,c("Row.names"="ensembl_gene_id"))
    de_results = merge(de_results,x,by.x = "genes", by.y="row.names",all.x=T)
    de_results = de_results[order(de_results$PValue),]
    
    
    
    write.csv(de_results,"result.csv",row.names = F)
    
    plot_heatmap_separate(x,sample_names,de_results,"tropak")    
    #return(de_results)
    
    go_analysis(lrt,prefix)
    
    #kegg_analysis(lrt,prefix)
    
}