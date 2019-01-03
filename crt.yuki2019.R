init = function()
{
    setwd("~/Desktop/work/wilson/")
    source("~/crt/crt.utils.R")    
}

plot_mds = function(file.png,y)
{
    png(file.png,res = 300,width=2000,height=2000)
    plotMDS(y,las=1)
    dev.off()    
}

tpm = read.tpm()
group = factor(rep(1,ncol(tpm)))

y=DGEList(counts=tpm,
          group=group,
          remove.zeros = T)

plot_mds("mds.all_samples.png",y)

# removing outliers
tpm$WL1101 = NULL
tpm$WL1723 = NULL
tpm$WL1724 = NULL

group = factor(rep(1,ncol(tpm)))

y=DGEList(counts=tpm,
          group=group,
          remove.zeros = T)

plot_mds("mds.no_outliers.png",y)

# only protein coding transcripts
protein_coding_transcripts = genes_transcripts[genes_transcripts$Ensembl_gene_id %in% protein_coding_genes.ens_ids$ENS_GENE_ID,]
tpm = tpm[row.names(tpm) %in% protein_coding_transcripts$Ensembl_transcript_id,]

group = factor(rep(1,ncol(tpm)))

y=DGEList(counts=tpm,
          group=group,
          remove.zeros = T)

plot_mds("mds.no_outliers.protein_coding_genes.png",y)

tpm$mean = rowMeans(tpm)
tpm = merge(tpm,genes_transcripts,by.x="row.names",by.y="Ensembl_transcript_id",all.x=T,all.y=F)
row.names(tpm) = tpm$Row.names
tpm$Row.names = NULL

tpm = tpm [order(tpm$mean,decreasing = T),]

tpm100 = head(tpm,100)
write.csv(tpm100,"tpm100.csv",quote = F)
