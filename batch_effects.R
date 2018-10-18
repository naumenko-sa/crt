normalize = function() {
  raw_counts = read.csv("test.csv", check.names = F, row.names = 1)
  
  #Log transform values
  for (coll in colnames(raw_counts)) {
    raw_counts[,coll] = log(raw_counts[,coll] + 1)
  }
  # remove any rows which have 0 variance
  raw_counts = raw_counts[ - as.numeric(which(apply(raw_counts, 1, var) == 0)),]
  raw_counts = as.matrix(raw_counts, row.names = 1)

  pheno = read.csv("pheno.csv", row.names = 1)
  pheno = as(pheno, "AnnotatedDataFrame") 
  pheno = pData(pheno)
  
  batch = pheno$batch
  modcombat = model.matrix(~1, data=pheno)
  
  combat_rpkm = ComBat(dat=raw_counts, batch=batch, mod=modcombat, par.prior = TRUE, prior.plots = TRUE)
  write.csv(combat_rpkm, file="combat_rpkm.csv")
  
  for (coll in colnames(combat_rpkm)) {
    combat_rpkm[,coll] = exp(combat_rpkm[,coll]) - 1 
  }
  
  write.csv(combat_rpkm, "log_untransformed_rpkm.csv")
  
}

