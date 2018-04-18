#read the list of samples
#parse the gct file with junction counts
#output the list of non-zero junctions with sum junction counts

args = commandArgs(trailingOnly = T)

#setwd("~/Desktop/work")

#download the list of samples - works only from browser - biobank inventory:
#https://www.gtexportal.org/home/datasets

#biobank = read.delim("biobank_collection_20180417_070153.txt", stringsAsFactors=F)
biobank = read.delim(args[1],stringsAsFactors=F)
biobank = biobank[biobank$materialType=="RNA:Total RNA",]
biobank = biobank[biobank$tissueSiteDetail=="Muscle - Skeletal",]
samples = gsub("-",".",biobank$tissueSampleId)

#GTEx_junctions = read.delim("GTEx_junctions.txt", skip = 2, stringsAsFactors=FALSE)
GTEx_jucntions = read.delim(args[2], skip = 2, stringsAsFactors=F)
cnames = colnames(GTEx_junctions)
for (i in 3:length(cnames))
{
    cnames[i] = substr(cnames[i],1,15)
}

colnames(GTEx_junctions) = cnames

#not all muscular samples present in junction table
#sample present in junctions: GTEX.11ZUS.1426
#sample absent in junctions: GTEX.S33H.2226
samples = samples[samples %in% colnames(GTEx_junctions)]

junction_ids = GTEx_junctions$junction_id

GTEx_junctions = subset(GTEx_junctions,select = samples)
GTEx_counts = rowSums(GTEx_junctions)

GTEx = data.frame(row.names = junction_ids)
GTEx = cbind(GTEx,GTEx_counts)

GTEx = subset(GTEx,GTEx_counts > 4)

write.table(GTEx,"junction_counts.txt")
