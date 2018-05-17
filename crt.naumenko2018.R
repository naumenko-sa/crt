# Function for the RNA-seq article Gonorazky.Naumenko.2018
install = function()
{
    source("http://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
}

init = function()
{
    library(edgeR)
    library(RColorBrewer)
    library(org.Hs.eg.db)
    source("~/crt/crt.utils.R")
    setwd("~/Desktop/work")
    
    #gene_lengths = read.delim("~/Desktop/project_muscular/reference/gene_lengths.txt", stringsAsFactors=F, row.names=1)
}


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

dexpression = function()
{
    counts = read.feature_counts_dir(update = T)
    
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

# output expression outliers as a table eoutliers.txt
# test is within cohort, fold change reported vs cohort and gtex
# this test is too sensitive to apply for all protein coding genes
expression.outliers.table = function(for_panels = T, file.rpkms)
{
    # DEBUG
    #for_panels = T
    #setwd("~/Desktop/work")
    #file.rpkms = "rpkms.muscle.txt"
    # DEBUG
    
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
    }else{
        gene_panel_name = "protein_coding_genes"
        gene_list = protein_coding_genes$Gene_name
    }
    
    for (sample in samples)
    {
        for (gene in gene_list)
        {
            expression4gene.table(gene,sample,counts,gene_panel_name,output)
        }
    }
}

expression4gene.table = function(gene,sample,counts,gene_panel_name,output)
{
    #gene = "AQP1"
    #sample = "S08_5.1.F"
    
    gene_expression = counts[counts$external_gene_name == gene,]
    
    if (nrow(gene_expression)>1)
    {
        gene_expression = head(gene_expression,1)
    }
    
    if (nrow(gene_expression)>0)
    {
        gene_expression$external_gene_name = NULL
        v_muscle = as.numeric(gene_expression[1,])
        muscle_mean.cohort = mean(v_muscle)
    
        muscle_mean.gtex = as.numeric(gtex_rpkm[gtex_rpkm$gene_name %in% c(gene),]$GTEX)
    
        expression_sample = gene_expression[[sample]]
        #significance
        ttest = t.test(v_muscle,mu=expression_sample)
    
        if ((muscle_mean.cohort > 0) && (expression_sample > 0)){
            fold_change.cohort = (max(muscle_mean.cohort, expression_sample)/min(muscle_mean.cohort,expression_sample))
            fold_change.gtex = (max(muscle_mean.gtex, expression_sample)/min(muscle_mean.gtex,expression_sample))
        
            if (expression_sample > muscle_mean.cohort){
                regulation = "UP"
            }else{
                regulation = "DOWN"
            }
        
            #if ((fold_change.cohort > 1.5 || fold_change.gtex > 1.5) && (ttest$p.value < 0.01)){
            if (fold_change.cohort > 2 && (ttest$p.value < 0.01)){
                cat(paste(sample,gene_panel_name,gene,regulation,fold_change.cohort,fold_change.gtex,ttest$p.value,sep = ","),file = output,append=T,sep="\n")   
            }
        }
    }
}

#detect expression outliers among similar samples and plot pictures
expression_outliers.pictures()
{
    setwd("~/Desktop/work")
    #read.rpkm_counts_dir()
    counts = read.table("rpkms.muscle.txt")
    #remove second entry for CLN3
    counts = counts[setdiff(rownames(counts),"ENSG00000261832"),]
    #sample = "S12_9.1.M"
    samples=head(colnames(counts),-1)
    #not GTEX
    #samples=tail(samples,-10)
    
    for (sample in samples)
    {
        dir.create(sample)
        setwd(sample)
        for (gene_panel_name in c("congenital_myopathy","congenital_muscular_dystrophies","congenital_myastenic_syndromes",
                         "channelopathies","vacuolar_and_others","limb_girdle","distal_myopathies","muscular_dystrophies"))
        {
            gene_panel = get(gene_panel_name)
            dir.create(gene_panel_name)
            setwd(gene_panel_name)
            print(gene_panel_name)
            expression4gene_in_a_gene_panel(sample,gene_panel,counts)
            setwd("..")
        }
        setwd("..")
    }
}

# table S10 - what outliers in muscle we can detect with myotubes
expression.outliers.trio_analysis = function()
{
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


# among genes expressed in muscle, what % is expressed in myo w > 10x coverage
# trios = families 9, 12, 14-1, 14-2, 17
# Table S12. 
expression.tableS12.tissues.all_genes.average = function()
{
    trios = paste0(c("12.1","14.1","14.2","17.1","18.1","26.1","28.1","5.1","6.1","9.1"),".Myo")
    rpkms.muscle = read.csv("rpkms.muscle.txt", sep="", stringsAsFactors=F)
    rpkms.muscle$"12.1.M" = (rpkms.muscle$S20_12.1.Mpv + rpkms.muscle$S21_12.1.Mvl)/2
    rpkms.muscle$S20_12.1.Mpv = NULL
    rpkms.muscle$S21_12.1.Mvl = NULL
    
    colnames(rpkms.muscle) = gsub("S.._","",colnames(rpkms.muscle))
    colnames(rpkms.muscle) = gsub("M","Myo",colnames(rpkms.muscle))
    
    rpkms.muscle = subset(rpkms.muscle,select = c(trios,"external_gene_name"))
    rpkms.muscle = rpkms.muscle[rpkms.muscle$external_gene_name %in% protein_coding_genes$Gene_name,]
    
    rpkms.muscle$mean = rep(0,nrow(rpkms.muscle))
    
    for (sample in trios)
    {
        rpkms.muscle$mean = rpkms.muscle$mean + subset(rpkms.muscle,select=sample)
    }
    rpkms.muscle$mean = rpkms.muscle$mean / length(trios)
    
    rpkms.muscle = rpkms.muscle[rpkms.muscle$mean > 1,]
    
    print(paste0(nrow(rpkms.muscle)," genes expressed in muscle at 1RPKM"))
    
    genes_expressed_in_muscle = row.names(rpkms.muscle)
    
    #5.1 is calculating
    trios.myo = paste0(c("12.1","14.1","14.2","17.1","18.1","26.1","28.1","6.1","9.1"),".Myo")
    rpkms.myo = read.csv("rpkms.myo.txt",sep="",stringsAsFactors = F)
    colnames(rpkms.myo) = gsub("S.._","",colnames(rpkms.myo))
    rpkms.myo = subset(rpkms.myo,select = c(trios.myo,"external_gene_name"))
    rpkms.myo = rpkms.myo[rpkms.myo$external_gene_name %in% protein_coding_genes$Gene_name,]
    rpkms.myo$mean = rep(0,nrow(rpkms.myo))
    for (sample in trios.myo)
    {
        rpkms.myo$mean = rpkms.myo$mean + subset(rpkms.myo,select=sample)
    }
    rpkms.myo$mean = rpkms.myo$mean / length(trios.myo)
    rpkms.myo = rpkms.myo[rpkms.myo$mean > 1,]
    genes_expressed_in_myo = row.names(rpkms.myo)
    
    print(paste0(nrow(rpkms.myo)," genes expressed in myo at 1RPKM"))
    
    test = rpkms.myo[row.names(rpkms.myo) %in% genes_expressed_in_muscle,]
    print(paste0(nrow(test)," muscle genes expressed in myo at 1RPKM"))
    
    rpkms.fibro = read.csv("rpkms.fibro.txt", sep="", stringsAsFactors = F)
    colnames(rpkms.fibro) = gsub("S.._","",colnames(rpkms.fibro))
    colnames(rpkms.fibro) = gsub("F","Myo",colnames(rpkms.fibro))
    rpkms.fibro = subset(rpkms.fibro,select = c(trios,"external_gene_name"))
    rpkms.fibro = rpkms.fibro[rpkms.fibro$external_gene_name %in% protein_coding_genes$Gene_name,]
    rpkms.fibro$mean = rep(0,nrow(rpkms.fibro))
    for (sample in trios)
    {
        rpkms.fibro$mean = rpkms.fibro$mean + subset(rpkms.fibro,select=sample)
    }
    rpkms.fibro$mean = rpkms.fibro$mean / length(trios)
    rpkms.fibro = rpkms.fibro[rpkms.fibro$mean > 1,]
    
    print(paste0(nrow(rpkms.fibro)," genes expressed in fibro at 1RPKM"))
    
    test = rpkms.fibro[row.names(rpkms.fibro) %in% genes_expressed_in_muscle,]
    print(paste0(nrow(test)," muscle genes expressed in fibro at 1RPKM"))
    
    test = rpkms.fibro[row.names(rpkms.fibro) %in% genes_expressed_in_myo,]
    print(paste0(nrow(test)," myo genes expressed in fibro at 1RPKM"))
    
    test = test[row.names(test) %in% genes_expressed_in_muscle,]
    print(paste0(nrow(test)," muscle genes expressed in fibro and myo at 1RPKM"))
    
    # coverage
    genes_expressed_in_myo = sort(unique(myo.expressed$external_gene_name))
    
    myo_coverage = read.delim("~/Desktop/work/S13_9-1-Myo.bam.coverage", comment.char="!", stringsAsFactors=FALSE)
    S58_12.1.Myo = read.delim("~/Desktop/work/S58_12-1-Myo.bam.coverage", comment.char="!", stringsAsFactors=FALSE)
    S23_14.1.Myo = read.delim("~/Desktop/work/S23_14-1-Myo.bam.coverage", comment.char="!", stringsAsFactors=FALSE)
    S24_14.2.Myo = read.delim("~/Desktop/work/S24_14-2-Myo.bam.coverage", comment.char="!", stringsAsFactors=FALSE)
    S45_17.1.Myo = read.delim("~/Desktop/work/S45_17-1-Myo.bam.coverage", comment.char="!", stringsAsFactors=FALSE)
    
    myo_coverage = cbind(myo_coverage,S58_12.1.Myo$avg,S23_14.1.Myo$avg,S24_14.2.Myo$avg,S45_17.1.Myo$avg)
    
    myo_coverage$mean = (myo_coverage$avg + myo_coverage$`S23_14.1.Myo$avg` + myo_coverage$`S58_12.1.Myo$avg` +
                             myo_coverage$`S24_14.2.Myo$avg` + myo_coverage$`S45_17.1.Myo$avg`)/5
    
    myo_coverage = myo_coverage[myo_coverage$gene %in% genes_expressed_in_myo,]
    myo_coverage = myo_coverage[myo_coverage$mean > 10,]
    
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
gor1.expression_variability()
{
    gor1_panel = c("C2","CACNB2","CD70","CKB","CTGF","CXCL14", "EBF2", "EPB41L3", "F2R", "FABP3", 
                   "FAT3", "HES1","HSPB7", "IGFBP7", "ITGA8", "KCNE4", "MAMDC2", "NDUFA4L2","PCDHGB6", "PCSK9",
                   "PI16", "PNRC2", "PPP1R14A", "PRLR", "PRUNE2", "RHOB", "SCUBE3", "SDK2", "SPINT2", "STYK1", 
                   "TPD52L1")
    rpkm.file = "rpkms.fibro.txt"
    tissue = "fibro"
    expression_variability(gor1_panel,rpkm.file,tissue)
}

expression.variability.tableS12()
{
    panel = c("DMD","NEB","ACTA1","LMNA","RYR1","MTM1","KCNIP4","FKRP","CAV3","LAMA2")
    for (tissue in c("fibro","myo","muscle"))
    {
        rpkm.file = paste0("rpkms.",tissue,".txt")
        expression_variability(panel,rpkm.file,tissue)
    }
}

# calculate a table for expression variabiliry in a gene_panel (list of genes) for 
# file like rpkms.muscle.txt - produced by crt.utils.read.coverage2counts_dir
expression_variability = function(gene_panel,rpkm.file,tissue)
{
    #rpkm.file = "rpkms.fibro.txt"
    #tissue = "fibro"
    #gene_panel = gor1_panel
    rpkms = read.table(rpkm.file)
    
    cnames = c("Tissue","Gene","Min","Mean","Max","SD","Var")
    
    result = data.frame()
    for (gene in gene_panel)
    {
        gene_expression = rpkms[rpkms$external_gene_name==gene,]
        gene_expression$external_gene_name=NULL
        if (nrow(gene_expression) == 1)
        {
            df_temp = data.frame(tissue,gene,min(gene_expression),rowMeans(gene_expression)[[1]],
                                 max(gene_expression),sd(gene_expression),var(as.numeric(gene_expression)))
            result = rbind(result,df_temp)
        }
    }
    colnames(result) = cnames
    write.table(result,paste0("expression_variability.",tissue,".csv"),sep = ",",row.names = F)
}
