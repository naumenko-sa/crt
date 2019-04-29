# batch effect correction with SVA
# https://bioconductor.org/packages/release/bioc/html/sva.html
# https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
installation <- function() {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("sva", version = "3.8")
    BiocManager::install("bladderbatch", version = "3.8")
}

init <- function() {
    library(sva)
    #library("Biobase")
    library(tidyverse)
}

tutorial <- function(){
    #install.packages("pamr")
    library(sva)
    library(bladderbatch)
    data(bladderdata)
    library(pamr)
    library(limma)
    
    pheno <- pData(bladderEset)
    edata <- exprs(bladderEset)
    
    mod <- model.matrix(~as.factor(cancer), data = pheno)
    mod0 <- model.matrix(~1, data = pheno)
    
    n.sv <- num.sv(edata, mod, method="leek")
    svobj <- sva(edata, mod, mod0, n.sv=n.sv)
}

# input
# - rpkm counts (why rpkm not raw?)
# - phenotype.csv 
#   sample names as rows, columns: sample number, outcome, batch and disease
#   sample number is unique id of sample, can be arbitrarily assigned (1,2,3...)
#   outcome can be "normal or disease"
#   batch is a description of all the batched (ex. 'CHEO or GTEx)
#   Disease is a value of either disease or normal 
#   number of rows of phenotype csv should be the same as number of columns of raw counts
correct_known_batches <-  function(){
    counts <-  read_csv("raw_counts_joined.coding.csv")
   
    sample_names <- tibble(sample_name = colnames(counts)) %>% tail(-1) %>% pull(sample_name)
    #Log transform values
    
    counts <- counts %>% 
        mutate_at(vars(-Ensembl_gene_id), log1p) %>% 
        mutate(sd = rowSds(as.matrix(.[sample_names]))) %>% 
        filter(sd != 0)
        
    raw_counts = as.matrix(raw_counts, row.names = 1)

    pheno = read.csv("pheno_nofibro.csv", row.names = 1)
    pheno = as(pheno, "AnnotatedDataFrame") 
    pheno = pData(pheno)
  
    batch = pheno$batch
    modcombat = model.matrix(~1, data = pheno)
    
    combat_counts <- ComBat(dat = raw_counts, 
                            batch = batch, 
                            mod = modcombat, 
                            par.prior = TRUE, 
                            prior.plots = TRUE)
    
    
    #write.csv(combat_rpkm, file="combat_rpkm.csv")
  
  for (coll in colnames(combat_rpkm)) {
    combat_rpkm[,coll] = exp(combat_rpkm[,coll]) - 1 
  }
  
  combat_rpkm[combat_rpkm < 0] = 0
  combat_rpkm = round(combat_rpkm, digits=0)
  write.csv(combat_rpkm, "nofibro_rpkms_adjusted.csv")
  
}

#Limma batch effect correction
#Input is rpkm counts of all samples
#batch vector is a column vector where each entry corresponds to control or sample 

batch_effects = function() {
  raw_counts = read.csv("nofibro.original.rpkms.csv", check.names = F, row.names = 1)
  #Log transform values
  for (coll in colnames(raw_counts)) {
    raw_counts[,coll] = log(raw_counts[,coll] + 1)
  }
  # remove any rows which have 0 variance
  raw_counts = raw_counts[ - as.numeric(which(apply(raw_counts, 1, var) == 0)),]
  y <- as.matrix(raw_counts)
  cheo = rep(c("cheo"), times = 21)
  gtexblood = rep(c("gtexblood"), times=9)
  batch = c(cheo, gtexblood)
  y2 <- removeBatchEffect(y, batch)
  par(mfrow=c(1,2))
  boxplot(as.data.frame(y),main="Original", las=2)
  boxplot(as.data.frame(y2),main="Batch corrected", las=2)
  
  for (coll in colnames(y2)) {
    y2[,coll] = exp(y2[,coll]) - 1 
  }
  y2[y2 < 0] = 0
  y2 = round(y2, digits=0)
  write.csv(y2, "nofibro_rpkm_adjusted.csv")
  
}

#PCA plots of adjusted and raw data

pca = function() {
    adjusted = read.csv("nofibro_rpkms_adjusted.csv", check.names = F, row.names = 1)
    raw_counts = read.csv("nofibro.original.rpkms.csv", check.names = F, row.names = 1)

  
  raw_counts.pca <- prcomp(t(raw_counts))
  summary(raw_counts.pca)
  plot(raw_counts.pca$x[, 1], raw_counts.pca$x[, 2], xlab = "PC1", ylab = "PC2", main="RPKM counts PCA plot")
  #text(raw_counts.pca$x[, 1], raw_counts.pca$x[, 2])
  points(raw_counts.pca$x[1:21,], col="orange", pch=16, cex=1)
  points(raw_counts.pca$x[22:31,], col="blue", pch=16, cex=1)
  text(raw_counts.pca$x[, 1], raw_counts.pca$x[, 2], labels = colnames(raw_counts))
  
  adjusted.pca = prcomp(t(adjusted))
  summary(adjusted.pca)
  plot(adjusted.pca$x[, 1], adjusted.pca$x[, 2], xlab = "PC1", ylab = "PC2", main="Limma Batch effect correction PCA (rpkm)")
  points(adjusted.pca$x[1:21,], col="orange", pch=16, cex=1)
  points(adjusted.pca$x[22:30,], col="blue", pch=16, cex=1)
  text(adjusted.pca$x[, 1], adjusted.pca$x[, 2], labels = colnames(adjusted))
  
  legend("topright",
         title="Sample",
         c("GTEx blood",
           "CHEO sample"
         ),
         fill=c("blue",
                "orange"))
  
}


#sva function is to estimate number of surrogate variables
#svaseq is specifically for sequencing data not microarray
#same inputs as combat
sva = function() {
  raw_counts = read.csv("nofibro_raw_counts_original.csv", check.names = F, row.names = 1)
  #Log transform values
  for (coll in colnames(raw_counts)) {
    raw_counts[,coll] = log(raw_counts[,coll] + 1)
  }
  # remove any rows which have 0 variance
  raw_counts = raw_counts[ - as.numeric(which(apply(raw_counts, 1, var) == 0)),]
  raw_counts = as.matrix(raw_counts, row.names = 1)
  
  pheno = read.csv("pheno_nofibro.csv", row.names = 1)
  pheno = as(pheno, "AnnotatedDataFrame") 
  pheno = pData(pheno)
  mod = model.matrix(~as.factor(disease), data=pheno)
  mod0 = model.matrix(~1, data=pheno)
  n.sv = num.sv(raw_counts,mod,method="leek")
  svobj = sva(raw_counts,mod,mod0,n.sv=n.sv)
  modsv = cbind(mod, svobj$sv)
  
  svaseq = svaseq(raw_counts,mod, n.sv=n.sv)
  
}
