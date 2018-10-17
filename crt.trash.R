junk=function()
{
    png("Fig2.total_counts_after_filtration.png",width=1000)
    op <- par(mar = c(10,4,4,2) + 0.1)
    barplot(all.samples.total_counts,las=2,main="Total counts after filtration")
    par(op)
    dev.off()
    
    png("Fig3.Genes_with_reads_after_filtration.png",width=1000)
    op <- par(mar = c(10,4,4,2) + 0.1)
    barplot(all.samples.non_zero,las=2,main="Genes with reads after filtration")
    par(op)
    dev.off()
    
    bars = list()
    for (col_name in colnames(all.samples.data))
    {
        #col_name = 'X1258.AC.A79'
        col_data = all.samples.data[,col_name]
        col_data = log2(col_data[col_data!=0])
        bars = c(bars, list(col_data))
        #boxplot(col_data,las=2)
        #assign(paste0("bars$",col_name),col_data)
    }
    
    png("Fig4.Log2_counts_for_covered_genes.png",width=1000)
    op <- par(mar = c(10,4,4,2) + 0.1)
    boxplot(bars,names = colnames(all.samples.data),las=2,main="Log2_counts_for_covered_genes")
    par(op)
    dev.off()
    
    
    norm_counts = cbind(symbol,cpm(all_counts[,c(2,3,4,5,6,7,9,10,11,12,13,14,15)]))
    write.table(norm_counts,"norm_counts_all.txt",col.names=NA,quote=F)  
    
    gene_locus = read.delim("genes_locus.txt",stringsAsFactors = F,header=T)
    gene_panel = read.delim("genes_muscular.txt",stringsAsFactors = F,header=T)
    
    muscular.locus = muscular.samples[muscular.samples$symbol %in% unlist(gene_locus),]
    mh.locus = mh.samples[mh.samples$symbol %in% unlist(gene_locus),]
    muscular.locus$id=NULL
    muscular.locus$symbol=NULL
    mh.locus$id=NULL
    mh.locus$symbol=NULL
    
    png("Fig5.Log2_counts_for_locus_mh_samples.png",width=1000)
    boxplot(log2(mh.locus+1),main="Log2 counts for locus mh samples")
    dev.off()
    
    png("Fig6.Log2_counts_for_locus_muscular_samples.png",width=1000)
    boxplot(log2(muscular.locus+1))
    dev.off()
    
    muscular.panel = muscular.samples[muscular.samples$symbol %in% unlist(gene_panel),]
    mh.panel = mh.samples[mh.samples$symbol %in% unlist(gene_panel),]
    muscular.panel$id=NULL
    muscular.panel$symbol=NULL
    mh.panel$id=NULL
    mh.panel$symbol=NULL
    
    png("Fig7.Log2_counts_for_panel_mh_samples.png",width=1000)
    boxplot(log2(mh.panel+1),main="Log2 counts for locus mh samples")
    dev.off()
    
    png("Fig8.Log2_counts_for_panel_muscular_samples.png",width=1000)
    boxplot(log2(muscular.panel+1))
    dev.off()
    
    work_counts=muscular_genes
    work_counts = all_counts
    
    #exploration
    row.names(work_counts) = work_counts$symbol
    work_counts$id = NULL
    work_counts$symbol = NULL
    total_counts = colSums(work_counts)
    barplot(total_counts)
    log_counts = log2(work_counts)
    png("all_genes_counts.png")
    par(mar=c(10,3,1,1))
    boxplot(log_counts,las=2)
    dev.off()
    
    attach(work_counts)
    
    x=muscular_genes[c("muscle1","X1130.BD.B175")]
    
    x=muscular_genes[c(2,3,4,5,6,7,9,10,11,12,13,14,15)]
    group=factor(c(1,2,3,4,5,6,7,8,9,10,11,12,13))
    
    group=factor(c(1,2))
    
    y=DGEList(counts=x,group=group)
    
    
    y=calcNormFactors(y)
    
    plotMDS(y)
    bcv=0.2
    #from A.thaliana experiment
    dispersion=0.04
    
    design=model.matrix(~group)
    fit=glmFit(y,design,dispersion)
    lrt=glmLRT(fit,coef=2)
    
    write.table(merge(work_counts,lrt$table,by="row.names"),"muscle1_vs_1130_panel.txt",col.names=NA,quote=F)  
    
}

expression_unfiltered = function()
{
    annotated_combined <- read.delim("~/Desktop/project_muscular/counts/muscular_unfiltered/annotated_combined.counts", row.names=1, stringsAsFactors=FALSE)
    all_genes_lengths <- read.delim("~/Desktop/project_muscular/counts/muscular_unfiltered/all_genes_lengths", row.names=1, stringsAsFactors=FALSE)
    
    samples = c("myotubes5","fibroblast5")  
    counts = annotated_combined[,samples]
    
    counts = merge(counts,all_genes_lengths,by.x="row.names",by.y="row.names")   
    row.names(counts)=counts$Row.names
    counts$Row.names = NULL
    
    x=counts[samples]
    group = factor(c(1,1))
    y=DGEList(counts=x,group=group,remove.zeros = F)
    plotMDS(y)
    #generate counts in the other way
    rpkms = rpkm(y,counts$Length)
    
    rpkms = merge(rpkms,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(rpkms)=rpkms$Row.names
    rpkms$Row.names=NULL
    
    panel.rpkm = rpkms[rpkms$external_gene_name %in% congenital_myopathy,]
    row.names(panel.rpkm) = panel.rpkm$external_gene_name
    panel.rpkm$external_gene_name=NULL
    panel.rpkm$Gene_description=NULL
    
    png("congenital_myopaties.png",res=100,width=1000)
    pheatmap(panel.rpkm,treeheight_row=0,treeheight_col=0,cellwidth = 40,
             display_number =T,cluster_rows=T, cluster_cols=T,
             main="Congenital myopaties, RPKM/2, unfiltered",breaks = c(0,10,50,100,500,1000,5000,10000,15000),
             colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8))
    dev.off()
}

#2017-02-22: expression of DGKE gene for Hernan and Mathieu
expression_rpkm_blood12 = function()
{
    setwd("~/Desktop/project_muscular/9_Blood1_2/")
    blood1 = load_rpkm_counts("blood1.counts_for_rpkm.txt")
    blood2 = load_rpkm_counts("blood2.counts_for_rpkm.txt")
    
    blood2$external_gene_name = NULL
    
    rpkm.counts = merge_row_names(blood1,blood2)
    
    write.table(rpkm.counts,"blood1_2.rpkms.txt",quote=F,sep = "\t")
    #sometimes contains duplicate entries, delete them
    rpkm.counts = read.delim("blood1_2.rpkms.txt", row.names=1, stringsAsFactors=F)
    
    breaks = c(0,5,10,50,100,500,1000)
    
    dgke = c("DGKE")
    
    plot_panel(dgke,rpkm.counts,"blood1_2.DGKE.png","Blood1_2 expr(rpkm), DGKE gene",breaks)
}

expression_rpkm_muscle2 = function()
{
    setwd("~/Desktop/project_muscular/1_Family_V_chr19_Muscle2/expression/")
    s62_AF_S5 = load_rpkm_counts("62_AF_S5.rpkm")
    write.table(s62_AF_S5,"62_AF_S5.rpkms.txt",quote=F,sep = "\t")
    
    s1258_AC_A79 = load_rpkm_counts("1258-AC-A79.rpkm")
    s1275_BK_B225 = load_rpkm_counts("1275-BK-B225.rpkm")
    s1388_MJ_M219 = load_rpkm_counts("1388-MJ-M219.rpkm")
    
    s1258_AC_A79$external_gene_name = NULL
    s1275_BK_B225$external_gene_name = NULL
    s1388_MJ_M219$external_gene_name = NULL
    
    s62_AF_S5 = merge_row_names(s62_AF_S5,s1258_AC_A79)
    s62_AF_S5 = merge_row_names(s62_AF_S5,s1275_BK_B225)
    s62_AF_S5 = merge_row_names(s62_AF_S5,s1388_MJ_M219)
    
    write.table(s62_AF_S5,"62_AF_S5.rpkms.controls.txt",quote=F,sep = "\t")
    #sometimes contains duplicate entries, delete them
    s62_AF_S5 = read.delim("62_AF_S5.rpkms.controls.txt", row.names=1, stringsAsFactors=F)
    
    plot_all_panels(s62_AF_S5,gtex_rpkm)
    
    #plot linkage region panels
    rpkms = s62_AF_S5
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,350)
    plot_panel(linkage_region1, rpkms, gtex_rpkm, "1_linkage_region.png","Linkage region part 1 RPKM",breaks)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,260)
    plot_panel(linkage_region2, rpkms, gtex_rpkm, "2_linkage_region.png","Linkage region part 2 RPKM",breaks)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,260)
    plot_panel(linkage_region3, rpkms, gtex_rpkm, "3_linkage_region.png","Linkage region part 3 RPKM",breaks)
}

sample_11_1_K = function()
{
    glomerular_diseases = c("NPHS1","NPHS2","PLCE1","CD2AP","LAMB2","ACTN4","TRPC6","WT1", "LMX1B", "SMARCAL1", "COQ2", 
                            "PDSS2", "MTTL1", "SCARB2", "FN1", "COL4A5", "COL4A6", "COL4A3", "COL4A4", "ALMS1", "ARHGDIA", 
                            "MYH9", "GLA", "ANLN", "ARHGAP24", "INF2", "PAX2", "CRB2", "MYO1E", "APOL1", "ADCK4", "ALG1", 
                            "CUBN", "PDSS2", "PMM2", "PTPRO", "SCARB2", "ZMPSTE24", "WDR73", "FN1", "NLRP3", "APOA1", "FGA", 
                            "LYZ", "B2M", "LMX1B", "PLCG2", "LAMB2")
    
    
    bardet_biedl = c("BBS1", "BBS2", "ARL6", "BBS4", "BBS5", "MKKS","BBS7", "TTC8", "BBS9", "BBS10", "TRIM32", "BBS12", "MKS1", 
                     "CEP290","WDPCP", "SDCCAG8", "LZTFL1", "BBIP1", "IFT27")
    
    setwd("~/Desktop/project_muscular/11-1-M/")
    
    S11_1_K = load_rpkm_counts("11-1-K.rpkm")
    kidneya = load_rpkm_counts("kidneya.rpkm")
    kidneyb = load_rpkm_counts("kidneyb.rpkm")
    kidneyc = load_rpkm_counts("kidneyc.rpkm")
    kidneyd = load_rpkm_counts("kidneyd.rpkm")
    
    kidneya$external_gene_name=NULL
    kidneyb$external_gene_name=NULL
    kidneyc$external_gene_name=NULL
    kidneyd$external_gene_name=NULL
    
    
    kidneys = merge_row_names(S11_1_K,kidneya)
    kidneys = merge_row_names(kidneys,kidneyb)
    kidneys = merge_row_names(kidneys,kidneyc)
    kidneys = merge_row_names(kidneys,kidneyd)
    
    write.table(kidneys,"kidneys.rpkms.txt",quote=F,sep = "\t")
    
    gene_panel = glomerular_diseases
    gene_panel = bardet_biedl
    
    panel_rpkm = kidneys[kidneys$external_gene_name %in% gene_panel,]
    row.names(panel_rpkm) = panel_rpkm$external_gene_name
    panel_rpkm$external_gene_name=NULL
    panel_rpkm = panel_rpkm[order(row.names(panel_rpkm)),]
    breaks = c(0,5,10,20,30,40,50,100,200,300,400,500,600)
    
    file = "bardet_biedl.png"
    title = "Glomerular diseases gene panel for 11-1-K and controls"
    title = "Bardet-Bield gene panel for 11-1-K and controls"
    
    
    breaks = seq(0,10)
    png(file,res=300,width=2000,height=3000)
    pheatmap(panel_rpkm,treeheight_row=0,treeheight_col=0,cellwidth = 40,
             display_number =T,cluster_rows=F, cluster_cols=T,
             main=title,
             breaks=breaks,
             colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaks)-1),
             fontsize = 8)
    dev.off()
}