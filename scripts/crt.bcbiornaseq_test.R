# tutorial on 
# https://github.com/hbc/bcbioRNASeq

# installation
BiocManager::install("remotes") 
# not working in every network - timeout
BiocManager::install("hbc/bcbioRNASeq")

library(bcbioRNASeq)
loadRemoteData("https://github.com/hbc/bcbioRNASeq/raw/f1000v2/data/bcb.rda")

# count data
assayNames(bcb)
raw <- counts(bcb, normalized = F)
normalized <- counts(bcb, normalized = T)
tpm <- counts(bcb, normalized = "tpm")
rlog <- counts(bcb, normalized = "rlog")
vst <- counts(bcb, normalized = "vst")
saveData(raw, normalized, tpm, rlog, vst)
writeCounts(raw, normalized, tpm, rlog, vst)

plotTotalReads(bcb)
plotMappingRate(bcb)
plotExonicMappingRate(bcb)
plotIntronicMappingRate(bcb)
plotGenesDetected(bcb)
plotGeneSaturation(bcb)


prepareRNASeqTemplate()
