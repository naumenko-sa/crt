# modification of http://i12g-gagneurweb.informatik.tu-muenchen.de/gitlab/baderda/mitoMultiOmics_public/
# Kremer et al 2017 article on RNA-seq analysis in mitochondrial diseasesa

install = function()
{
  install.packages("beeswarm")
  install.packages("data.table")			# faster tables than data frame
  install.packages("DT")
  install.packages("gplots")
  install.packages("knitr")
  install.packages("LSD")
  install.packages("png")
  install.packages("plotrix")
  install.packages("RColorBrewer")
  install.packages("rmarkdown")
  install.packages("stringr")
  install.packages("XLConnect")
  install.packages("foreign")
  
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
}

load_packages = function(){
    library(beeswarm)
    library(data.table)			# faster tables than data frame
    library(DESeq2)
    library(DT)
    library(gplots)
    library(knitr)
    library(LSD)
    library(png)
    library(plotrix)
    library(RColorBrewer)
    library(rmarkdown)
    library(stringr)
    library(XLConnect)
}