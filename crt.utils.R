library(pheatmap)
library(RColorBrewer)
library(edgeR)   
library(readr)

trios = c("12.1","14.1","14.2","17.1","18.1","26.1","28.1","5.1","6.1","9.1","40.1","40.2","4.1")

mh_panel=c("CACNA1S","RYR1","STAC3","TRDN","ASPH","JPH2","CASQ1","ATP2A1","ATP2A2","CALM1","FKBP1A")

#panels
panel_list = c("congenital_myopathy",
               "congenital_muscular_dystrophies",
               "congenital_myastenic_syndromes",
               "channelopathies", 
               "vacuolar_and_others",
               "limb_girdle",
               "distal_myopathies",
               "muscular_dystrophies")

congenital_myopathy = c("ACTA1","BIN1","CACNA1S","CCDC78","CFL2","CNTN1","DNM2","KBTBD13","KLHL40","KLHL41",
                        "LMOD3","MEGF10","MLTK","MTM1", "MTMR14","MYH2","MYH7","MYO18B","MYPN","NEB",
                        "ORAI1","PTPLA","RYR1","SEPN1","SPEG","STAC3","STIM1","TNNT1","TPM2","TPM3", 
                        "TTN")

congenital_muscular_dystrophies=c("ACTA1","ALG13","B3GALNT2","B3GNT1","CHKB","COL6A1","COL6A2","COL6A3","DNM2","DPM1",
                                  "DPM2","FHL1","FKRP","FKTN","GMPPB","ISPD","ITGA7","LAMA2","LARGE","LMNA",
                                  "POMGNT1","POMGNT2","POMK","POMT1","POMT2","SEPN1","TCAP","TMEM5","TRAPPC11")

congenital_myastenic_syndromes = c("AGRN","ALG14","ALG2","CHAT","CHRNA1","CHRNB1","CHRND","CHRNE","CHRNG","COLQ",
                                   "DOK7","DPAGT1","GFPT1","LAMB2","MUSK","PLEC","PREPL","RAPSN","SCN4A","SYT2")

channelopathies = c("ATP2A1","CACNA1A","CACNA1S","CAV3","CLCN1","CNBP","DMPK","HSPG2","KCNA1","KCNE3",
                    "KCNJ12","SCN4A")

vacuolar_and_others = c("ACVR1","CAV3","CLN3","FHL1","ISCU","LAMP2","MCOLN1","MSTN","PABPN1","PLEC",
                        "TTN","VMA21")

limb_girdle = c("ANO5","CAPN3","CAV3","DAG1","DES","DNAJB6","DPM3","DYSF","FKRP","FKTN",
                "GAA","GMPPB","HNRNPDL","ISPD","LIMS2","LMNA","MYOT","PLEC","POMGNT1","POMT1",
                "POMT2","SGCA","SGCB","SGCD","SGCG","TCAP","TNPO3","TRAPPC11","TRIM32","TTN",
                "VCP" )

distal_myopathies=c("ANO5","BAG3","CRYAB","DES","DNM2","DYSF","FLNC","FLNC","GNE","KLHL9",
                    "LDB3","LDB3","MATR3","MYH7","MYOT","NEB","SEPN1","TIA1","TRIM54","TRIM63",
                    "TTN","VCP")

muscular_dystrophies= c("DMD","DUX4","EMD","FHL1","LMNA","PTRF","SMCHD1","SYNE1","SYNE2","TMEM43",
                        "TOR1AIP1")

#gene panels for linkage region on chr19
linkage_region1 = c("GAMT", "DAZAP1", "RPS15", "APC2", "C19orf25", "PCSk4", "REEP6","ADAMTSL5","PLK5","MEX3D", 
                    "MBD3", "UQCR11", "TCF3", "ONECUT3", "ATP8B3", "REX01", "MIR1909", "KLF16", "ABHD17A", "SCAMP4", 
                    "ADAT3", "CSNK1G2", "CSNK1G2-AS1", "BTBD2", "MKNK2", "MOB3A", "IZUMO4", "AP3D1", "DOT1L", "PLEKHJ1")

linkage_region2 = c("MIR1227", "MIR6789", "SF3A2", "AMH", "MIR4321", "JSRP1", "OAZ1", "C19orf35", "LINGO3", "LSM7", 
                    "SPPL2B", "TMPRSS9", "TIMM13", "LMNB2", "MIR7108", "GADD45B", "GNG7", "MIR7850", "DIRAS1", "SLC39A3",
                    "SGTA", "THOP1", "ZNF554", "ZNF555", "ZNF556", "ZNF57", "NF77", "TLE6", "TLE2", "MIR1268A")

linkage_region3 = c("AES", "GNA11", "GNA15", "S1PR4", "NCLN", "CELF5", "NFIC", "SMIM24", "DOHH", "FZR1", 
                    "SNORD38", "FZR1", "C19orf71", "MFSD12", "HMG20B", "GIPC3", "TBXA2R", "CACTIN-AS1", "CACTIN", 
                    "PIP5K1C", "TJP3", "APBA3", "MRPL54", "RAX2", "MATK", "ZFR2", "ATCAY", "NMRK2", "DAPK3")

linkage_region4 = c("MIR637", "EEF2", "SNORD37", "PIAS4", "ZBTB7A", "MAP2K2", "CREBL3L3", "SIRT6", "ANKRD24", "EBI3",
                    "CCDC94", "SHD", "TMIGD2", "FSD1", "STAP2", "MPND", "SH3GL1", "CHAF1A","UBXN6", "MIR4746", 
                    "PLIN4", "PLIN5", "LRG1", "SEMA6B", "TNFAIP8L1", "MYDGF", "DPP9", "DPP9-AS1", "MIR7-3HG", "MIR7-3")

linkage_region5 = c("FEM1A", "TICAM1","PLIN3", "ARRDC5", "UHRF1", "MIR4747", "KDM4B", "PTPRS", "ZNRF4", "TINCR", 
                    "SAFB2", "SAFB", "C19orf70", "HSD11B1L", "RPL36", "LONP1", "CATSPERD", "PRR22", "DUS3L","NRTN",
                    "FUT6","FUT5","NDUFA11","VMAC","CAPS","RANBP3", "RFX2","ACSBG2","MLLT1","ACER1")

linkage_region6 = c("CLPP", "ALKBH7", "PSPN", "GTF2F1", "MIR6885", "MIR6790", "KHSRP", "MIR3940", "SLC25A41", "SLC25A23",
                    "CRB3", "DEND1C", "TUBB4A", "TNFSF9", "CD70", "TNFSF14", "C3", "GPR108", "MIR6791", "TRIP10", 
                    "SH2D3A", "VAV1", "ADGRE1", "MBD3L5", "MBD3L4", "MBD3L2", "MBD3L3", "ZNF557", "INSR",  "ARHGEF18")

linkage_region7 = c("PEX11G", "C19orf45", "ZNF358", "MCOLN1", "PNPLA6", "CAMPSAP3", "MIR6792", "XAB2", "PET100", "PCP2",  
                    "STXBP2", "RETN", "MCEMP1", "TRAPPC5", "FCER2", "CLEC4G", "CD209", "CLEC4M", "EVI5L", "PRR36",
                    "LRRC8E", "MAP2K7","SNAPC2","CTXN1","TIMM44","ELAVL1","CCL25","FBN3","CERS4","CD320")

linkage_region8 = c("NDUFA7", "RPS28", "KANK3", "ANGPTL4", "RAB11B-AS1", "MIR4999", "RAB11B", "MARCH2", "HNRNPM", 
                    "PRAM1", "ZNF414", "MYO1F", "ADAMTS10", "ACTL9", "OR2Z1", "ZNF558", "MBD3L1", "OR1M1", "MUC16")

protein_coding_genes <- read.table("~/cre/data/protein_coding_genes.txt", stringsAsFactors=F, header=T)

protein_coding_genes.bed = read.delim("~/cre/data/protein_coding_genes.bed", header=F, stringsAsFactors=F)
colnames(protein_coding_genes.bed) = c("chrom","start","end","gene")

protein_coding_genes.ens_ids <- read.table("~/cre/data/protein_coding_genes.ens_ids.txt", stringsAsFactors=F,header=T)

omim.file = "~/Desktop/reference_tables/omim_inheritance.csv"

if (file.exists(omim.file))
{
    omim = read.csv(omim.file, sep=";", stringsAsFactors=FALSE)
    omim = subset(omim,select=c("Gene","Omim"))
}

gtex_rpkm_file = "/home/sergey/Desktop/stories/4_RNAseq_diagnostics/rnaseq_article/figures/expression_plots/gtex.muscle_genes.rpkm.txt"
if (file.exists(gtex_rpkm_file))
{
    gtex_rpkm = read.csv(gtex_rpkm_file, sep="", stringsAsFactors = F)
}

ensembl_w_description = read.delim2("~/cre/data/ensembl_w_description.txt", row.names=1, stringsAsFactors=F)
#gene_lengths = read.delim("~/Desktop/project_RNAseq_diagnostics/reference/gene_lengths.txt", stringsAsFactors=F, row.names=1)

installation = function()
{
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    
    install.packages("pheatmap")
}


get_genes_in_panels = function()
{
    gene_list = c()
    for (gene_panel_name in panel_list)
    {
        gene_panel = get(gene_panel_name)
        gene_list = unique(c(gene_list,gene_panel))
    }
    return(gene_list)
}

# calculates RPKMs using ~/bioscripts/bam.raw_coverage.sh input - 
# usable to calculate RPKMs for exons using a bed file for page website
raw_coverage2rpkm = function(filename)
{
    #test:
    #filename="S57_32-1-M.bam.raw_coverage"
    
    counts = read_delim(filename,"\t", col_names = FALSE, 
                                 col_types = cols(X1 = col_character()))
    
    colnames(counts)=c("chrom","start","end","exon_id","coverage")
    
    Gene_lengths = counts$end - counts$start + 1
    rownames(counts)=counts$exon_id
    
    counts$chrom = NULL
    counts$start = NULL
    counts$end = NULL
    counts$exon_id = NULL
  
    counts = rpkm(counts,Gene_lengths)
  
    colnames(counts) = gsub(".bam.raw_coverage","",filename)
  
    return(counts)
}

read.raw_coverage2rpkm_dir = function()
{
    files = list.files(".","*raw_coverage")
    counts = raw_coverage2rpkm(files[1])
    for (file in tail(files,-1))
    {
        print(file)
        counts_buf = raw_coverage2rpkm(file)
        counts = merge_row_names(counts,counts_buf)
    }
    write.table(counts,"exon.rpkms.txt",quote=F)
}

# Loads coverage from bam.coverage_bamstats05.sh
# loads gene (exonic) lengths
# calculates RPKMs
# returns EXON_ID, rpkm
coverage2rpkm = function(filename)
{
    #test:
    #filename="S08_5-1-F.bam.coverage"
    library(edgeR)   
    #first line in the file is a comment
    counts = read.delim(filename, stringsAsFactors=F)
    Gene_lengths = counts$length
    counts = subset(counts, select=c("gene","mean"))
    row.names(counts)=counts$gene
    counts$gene = NULL
    
    counts = rpkm(counts,Gene_lengths)
  
    colnames(counts) = gsub(".bam.coverage","",filename)
  
    return(counts)
}

read.coverage2counts_dir = function(update=F)
{
    if(file.exists("rpkms.txt") && update == F)
    {
        counts = read.table("rpkms.txt")
    }
    else
    {
        files = list.files(".","*coverage")
        counts = coverage2rpkm(files[1])
        for (file in tail(files,-1))
        {
            counts_buf = coverage2rpkm(file)
            counts = merge_row_names(counts,counts_buf)
        }
        write.table(counts,"rpkms.txt",quote=F)
    }
    return(counts)
}

# Loads a file with counts from feature_counts
# loads gene lengths
# loads gene names
# calculates RPKMs
# returns ENS_ID, rpkm, Gene_name
feature_counts2rpkm = function(filename)
{
    #test:
    #filename="S01_1-1-B.bam.rpkm_counts.txt"
    #first line in the file is a comment
    counts = read.delim(filename, stringsAsFactors=F, row.names=1,skip=1)
    counts$Chr=NULL
    counts$Start=NULL
    counts$End=NULL
    counts$Strand=NULL
    
    Gene_lengths = counts$Length
    counts$Length=NULL
    
    counts = rpkm(counts,Gene_lengths)
    colnames(counts) = gsub(".bam","",colnames(counts))
    
    return(counts)
}

#reads all counts in the current directory
read.rpkm_counts_dir = function(update=F)
{
    if(file.exists("rpkms.txt") && update == F)
    {
        counts = read.table("rpkms.txt")
    }
    else
    {
        files = list.files(".","*rpkm_counts.txt")
        counts = feature_counts2rpkm(files[1])
        for (file in tail(files,-1))
        {
            print(file)
            counts_buf = feature_counts2rpkm(file)
            counts = merge_row_names(counts,counts_buf)
        }
        
        counts = merge(counts,ensembl_w_description,by.x="row.names",by.y="row.names")
        row.names(counts)=counts$Row.names
        counts$Row.names=NULL
        counts$Gene_description=NULL
    
        #remove a second entry of CLN3 from rpm_counts.txt
        counts = counts[!row.names(counts) %in% c("ENSG00000261832","ENSG00000267059","ENSG00000200733","ENSG00000207199","ENSG00000252408",
                                                  "ENSG00000212270","ENSG00000212377","ENSG00000167774"),]
        
        write.table(counts,"rpkms.txt",quote=F)
    }
    return(counts)
}

read.feature_counts = function(filename)
{
    #first line in the file is a comment
    counts = read.delim(filename, stringsAsFactors=F, row.names=1,skip=1)
    counts$Chr=NULL
    counts$Start=NULL
    counts$End=NULL
    counts$Strand=NULL
    counts$Length=NULL
  
    colnames(counts) = gsub(".bam","",colnames(counts))
  
    return(counts) 
}

#reads all counts in the current directory
read.feature_counts_dir = function(update=F)
{
    if(file.exists("raw_counts.txt") && update == F)
    {
        counts = read.table("raw_counts.txt")
    }
    else
    {
        files = list.files(".","*rpkm_counts.txt")
        counts = read.feature_counts(files[1])
        for (file in tail(files,-1))
        {
            print(paste0("Reading ",file))
            counts_buf = read.feature_counts(file)
            counts = merge_row_names(counts,counts_buf)
        }
        write.table(counts,"raw_counts.txt",quote=F)
    }
    return(counts)
}

# merge two dataframes by row.names and fix the row.names of the resulting df
merge_row_names = function(df1,df2)
{
    merged = merge(df1,df2,by.x='row.names',by.y='row.names')
    row.names(merged) = merged$Row.names
    merged$Row.names = NULL
    return(merged)
}

#  prepares an expression profile for GSEA
#  http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Expression_Data_Formats
#  for GSEA it is important to report all genes - genome wide
#  hopefully cpms are better than logcpms
prepare_file_4gsea = function(counts,samples,prefix)
{
  t_cpm = cpm(counts,prior.count=1,log=F)
  t_cpm = t_cpm[,samples]
  
  #remove rows with 0 expression
  t_cpm = t_cpm[rowSums(t_cpm)>0,]
  
  t_cpm = merge(protein_coding_genes,t_cpm,by.x='row.names',by.y='row.names',all = F)
  
  rownames(t_cpm)=t_cpm$Row.names
  t_cpm$Row.names = NULL
  
  result_file=paste0(prefix,".4gsea.txt")
  
  t_cpm$ensembl_gene_id = row.names(t_cpm)
  
  t_cpm = t_cpm[c("gene_name","ensembl_gene_id",paste0(samples))]
  colnames(t_cpm)[1:2]=c("NAME","DESCRIPTION")
  
  o = order(rowSums(t_cpm[,c(samples)]),decreasing = T)
  t_cpm = t_cpm[o,]
  d = duplicated(t_cpm$NAME)
  dy = t_cpm[d,]$NAME
  t_cpm = t_cpm[!d,]
  nrow(t_cpm)
  
  write.table(t_cpm,result_file,quote=F,row.names = F,sep = "\t")
}

# usually heatmap is a part of a panel - we don't need a title
# rows are not clustered to save alpabetical gene order
# cols are not clustered to save sample order
plot_heatmap = function(prefix,expression_table)
{
    rows = nrow(expression_table)
    cellheight = 10
    res = 300
    filename = paste0(prefix,".",res,"ppi.png")
  
    png(filename,res = res,height=rows * cellheight * 5+500,width=1500)
    pheatmap(expression_table,scale="row",treeheight_row=0,treeheight_col=0,
           display_numbers = T,cellheight = cellheight,cellwidth = 30,
           cluster_rows = F, cluster_cols = F)
    dev.off()
}

# plot separate heatmaps for upregulated and downregulated genes 
# parameters:
# counts - initial row count for all genes
# samples
# de_results
# prefix  
# ntop - how many top genes to plot
plot_heatmap_separate = function(counts,samples,de_results,prefix,ntop = NULL)
{
    logcpm = cpm(counts,prior.count=1,log=T)
    #cpm0  = log(cpm(y$counts+1))
    top_genes_cpm = logcpm[de_results$Ensembl_gene_id,]
    top_genes_cpm = top_genes_cpm[,samples]
    rownames(top_genes_cpm) = de_results$Gene_name
  
    #expressed higher in WNT-dependent cells.
    upregulated_genes = de_results[de_results$logFC<0,]$Gene_name
    downregulated_genes = de_results[de_results$logFC>0,]$Gene_name
  
    if (!is.null(ntop))
    {
        upregulated_genes = head(upregulated_genes,ntop)
        downregulated_genes = head(downregulated_genes,ntop)
    }
  
    #sort genes alphabetically - it is much easier to read heatmap
    #upregulated_genes = sort(upregulated_genes)
    #downregulated_genes = sort (downregulated_genes)
    # upregulated_genes has 1 gene, simple subsetting does not work
    plot_heatmap(paste0(prefix,".left"),subset(top_genes_cpm,row.names(top_genes_cpm) %in% upregulated_genes))
    plot_heatmap(paste0(prefix,".right"),subset(top_genes_cpm,row.names(top_genes_cpm) %in% downregulated_genes))
}   

# plots the heatmap of title.png for gene_panel using sample_rpkm and gtex_rpkm
plot_panel= function(gene_panel, sample_rpkm, filename,title, breaks,height=2000)
{
    #test:
    #breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,
    #  5000,6000,7000,8000,9000,10000,15000,17100)
    #gene_panel = linkage_region4
    #sample_rpkm = rpkms
    #filename = "1_congenital_myopaties.png"
    #title = "Congenital myopaties RPKM"
    
    panel_rpkm = sample_rpkm[sample_rpkm$external_gene_name %in% gene_panel,]
    gtex_panel_rpkm = gtex_rpkm[gtex_rpkm$gene_name %in% gene_panel,]
  
    all_rpkm = merge(panel_rpkm,gtex_panel_rpkm,by.x='external_gene_name',by.y='gene_name')
    row.names(all_rpkm) = all_rpkm$external_gene_name
    all_rpkm$external_gene_name=NULL
  
    png(filename,res=200,width=7000,height=height)
    pheatmap(all_rpkm,treeheight_row=0,treeheight_col=0,
           cellwidth = 40, cellheight = 10,
           display_number=T,cluster_rows=F, cluster_cols=T,
           main=title,
           breaks=breaks,
           colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaks)-1),
           fontsize = 8)
    dev.off()
}

# plot expression for 8 gene panels
plot_all_panels = function(rpkms)
{
    #rpkms = read.rpkm_counts_dir(update = F)
    all_genes = get_genes_in_panels()
    rpkms = read.table("rpkms.muscle.txt")
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,
               300,400,500,600,700,800,900,1000,2000,3000,4000,
               5000,6000,7000,8000,9000,10000,15000,17100,26000)
    plot_panel(all_genes, rpkms, 
               "muscular_genes.png",
               "Muscular genes, RPKM",
               breaks,
               height=4000)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,
               300,400,500,600,700,800,900,1000,2000,3000,4000,
               5000,6000,7000,8000,9000,10000,15000,17100,24000)
    plot_panel(congenital_myopathy, rpkms, 
               "1_congenital_myopaties.png",
               "Congenital myopaties gene panel, RPKM",
               breaks)
  
    breaks = c(0,5,10,20,30,40,50,100,220)
    plot_panel(congenital_myastenic_syndromes, rpkms, "2_congenital_myastenic_syndromes.png",
               "Congenital myastenic syndromes panel, RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,110)
    plot_panel(channelopathies, rpkms, "3_channelopathies.png","Channelopathies panel, RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000)
    plot_panel(vacuolar_and_others, rpkms, "4_vacuolar_and_others.png","Vacuolar and others panel, RPKM",breaks)
  
    breaks = c(0,5,10,50,100,500,1000,2000,5000,6000)
    plot_panel(distal_myopathies, rpkms, "5_distal_myopathies.png","Distal myopathies panel, RPKM",breaks)
  
    breaks = c(0,5,10,50,100,500,1000,2000,5000,10000,15000,17100)
    plot_panel(congenital_muscular_dystrophies, rpkms, "6_congenital_muscular_dystrophies.png","Congenital MD panel, RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,
             5000)
    plot_panel(limb_girdle, rpkms, "7_Limb_girdle_dystrophies.png","Limb_girdle_dystrophies panel, RPKM",breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,2400)
    plot_panel(muscular_dystrophies, rpkms, "8_Muscular_dystrophies.png",
               "Muscular_dystrophies panel, RPKM",breaks)
}

plot_region_panels = function(rpkms)
{
    rpkms = read.rpkm_counts_dir(update = T)
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,
             300,400)
    plot_panel(linkage_region1, rpkms, 
             "1_linkage_region.png",
             "Linkage region 1, RPKM",
             breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,
             300,400)
    plot_panel(linkage_region2, rpkms, 
             "2_linkage_region.png",
             "Linkage region 2, RPKM",
             breaks)
  
    plot_panel(linkage_region3, rpkms, 
             "3_linkage_region.png",
             "Linkage region 3, RPKM",
             breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,
               300,400,500,600,700,800)
    plot_panel(linkage_region4, rpkms, 
             "4_linkage_region.png",
             "Linkage region 4, RPKM",
             breaks)
    
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
    plot_panel(linkage_region5, rpkms, 
               "5_linkage_region.png",
               "Linkage region 5, RPKM",
               breaks)
  
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
    plot_panel(linkage_region6, rpkms, 
               "6_linkage_region.png",
               "Linkage region 6, RPKM",
               breaks)
    
    plot_panel(linkage_region7, rpkms, 
               "7_linkage_region.png",
               "Linkage region 7, RPKM",
               breaks)
    
    plot_panel(linkage_region8, rpkms, 
               "8_linkage_region.png",
               "Linkage region 8, RPKM",
               breaks)
}

#https://www.r-bloggers.com/summarising-data-using-dot-plots/
#expression_dotplot = function(rpkms,tissue_type,gene_panel_name,color)
#plot a picture like Beryl's S1
expression_dotplot = function(gene_panel_name)
{
    library(lattice)
    library(latticeExtra)
    
    #gene_panel_name = "channelopathies"
    gene_panel = get(gene_panel_name)
    
    colors = c("darkgreen","red","orange")
    
    rpkms_muscle = read.table("rpkms.muscle.txt")
    rpkms_myotubes = read.table("rpkms.myotubes.txt")
    rpkms_fibro = read.table("rpkms.fibro.txt")
    
    panel_rpkm_muscle = rpkms_muscle[rpkms_muscle$external_gene_name %in% gene_panel,]
    panel_rpkm_myotubes = rpkms_myotubes[rpkms_myotubes$external_gene_name %in% gene_panel,]
    panel_rpkm_fibro = rpkms_fibro[rpkms_fibro$external_gene_name %in% gene_panel,]
    
    row.names(panel_rpkm_muscle) = panel_rpkm_muscle$external_gene_name
    panel_rpkm_muscle$external_gene_name = NULL
    
    row.names(panel_rpkm_myotubes) = panel_rpkm_myotubes$external_gene_name
    panel_rpkm_myotubes$external_gene_name = NULL
    
    row.names(panel_rpkm_fibro) = panel_rpkm_fibro$external_gene_name
    panel_rpkm_fibro$external_gene_name = NULL
    
    df = data.frame(matrix(ncol=4,nrow=0))
    names(df) = c("gene","sample","expression","type")
    
    f_muscle_type = factor(c("muscle","myotubes","fibro"))
    
    for (m_type in f_muscle_type){
        if (m_type == "muscle"){
            panel_rpkm = panel_rpkm_muscle
        }else if (m_type == "myotubes"){
            panel_rpkm = panel_rpkm_myotubes
        }else{
            panel_rpkm = panel_rpkm_fibro   
        }
    
        f_gene = factor(sort(gene_panel,decreasing = T))  
        f_sample = factor(colnames(panel_rpkm))
    
        for (gene in f_gene){
            for (sample in f_sample){
                de = data.frame(gene,sample,panel_rpkm[gene,sample],m_type)
                names(de) = c("gene","sample","expression","type")
                df = rbind(df,de)
            }
        }
    }
 
    filename=paste0(gene_panel_name,"_dotplot.png")
    png(filename,res=300,width=2000,height=2000)

    #https://stackoverflow.com/questions/35489624/remove-grid-lines-in-dotplot-without-modifying-underlying-trellis-parameters    
    d1 <- trellis.par.get("dot.line")
    d1$lwd <- 0  ## hack -- set line width to 0
    trellis.par.set("dot.line",d1)
    
    parSettings = list(axis.line = list(lty = 3),
                       layout.widths = list(left.padding = 0,
                                            right.padding = 0,
                                            ylab.axis.padding = 0,
                                            between = 0),
                       axis.components = list (right = list (pad1 = 0, pad2 = 0),
                                               top = list(pad1 = 0, pad2 = 0))
    )
    
    #https://stat.ethz.ch/pipermail/r-help/2007-March/128502.html - removing frame borders
    plot1 = dotplot(
                gene ~ expression,
                data = df,
                col = colors,
                subset = expression <= 30,
                scales = list(x = list(relation = "sliced",
                                       tick.number = 3),
                              y = list(relation = "same")
                ),
                drop.unused.levels = F,
                group = type,
                pch = 16,
                par.settings = parSettings,
                xlab = NULL
    )
    
    #as.layer(dotplot(df$gene ~ df$expression,subset = df$expression >=100, col="blue"),x.same = F,opposite = F,)
            
    plot2 = dotplot(gene ~ expression,data=df, 
                    col=colors,
                    subset = ((expression > 30) & (expression <=100)),
                    scales = list(x = list(relation = "sliced"),
                                y = list(draw = F)),
                     drop.unused.levels = F,
                    group = type,
                    pch = 16,
                    par.settings = parSettings,
                    xlab = NULL
    )
    
    plot3 = dotplot(gene ~ expression,
                    data=df, 
                    col=colors, 
                    subset = expression >100,
                    scales = list(x=list(relation="sliced"),
                                  y=list(draw=F)),
                    drop.unused.levels = F,
                    group = type,
                    pch = 16,
                    par.settings = parSettings,
                    xlab = NULL
    )
    
    plot(plot1, split=c(1,1,3,1))
    plot(plot2, split=c(2,1,3,1),newpage = F)
    plot(plot3, split=c(3,1,3,1),newpage = F)
    #groups = expression < 100
            
    #for some reason you have to print plots in lattice
    #print(plot1)
    dev.off()
}

all_expression_dotplots = function()
{
    for (gene_panel_name in panel_list)
    {
        expression_dotplot(gene_panel_name)
    }
}

# plots coverage for every gene for all samples having sample.coverage in the current dir - output of
# bam.coverage.bamstats05.sh
coverage_plot = function ()
{
    setwd("~/Desktop/work")
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
    
    #remove ACTA1 with has coverage of 328115
    #coverage = coverage[-c(1),]
    
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
        png(paste0(bin,".png"),res=300,width=5000,height=2000)
        boxplot(t(coverage[get(bin),]),
            las=2,
            cex.axis=0.8)
        dev.off()
    }
}

decode_tissue_type = function (row)
{
    if (grepl("F",row[1])){
       tissue = "Fibroblast"
    }else if (grepl("Myo",row[1])){
       tissue = "Myotubes"
    }else{
        tissue = "Muscle" 
    }
    return(tissue)
}

splicing.read_novel_splice_events = function(file)
{
    ############################################################
    # debug
    # setwd("~/Desktop/work/splicing/muscle")
    # file = "S12_9-1-M.bam_specific_rc5_norm_rc0.05_n_gtex_184"
    ############################################################
    
    sample = strsplit(file,".",fixed=T)[[1]][1]
    print(sample)
    splice_events = read.csv(file,stringsAsFactors=FALSE)
    splice_events$sample = sample
    #splice_events = merge(splice_events,omim,by.x = "gene", by.y = "Gene",all.x = T)
    splice_events$Omim=""
    splice_events = splice_events[,c("sample","gene","pos","annotation","read_count","norm_read_count",
                                     "n_gtex_seen","total_gtex_read_count","Omim")]
    #write.csv(splice_events,file = paste0(sample,".csv"),row.names = F)
    return(splice_events)
}