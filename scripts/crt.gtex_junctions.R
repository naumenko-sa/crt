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

print(length(samples))

#GTEx_junctions = read.delim("GTEx_junctions.txt", skip = 2, stringsAsFactors=FALSE)
GTEx_junctions = read.delim(args[2], skip = 2, stringsAsFactors=F)

# there are duplicates - same junction, same counts, two different genes, 
# i.e. 1_33060825_33065687 
# ENSG00000254553.1
# ENSG00000160062.10

GTEx_junctions$Description=NULL
GTEx_junctions = unique(GTEx_junctions)

cnames = colnames(GTEx_junctions)
for (i in 2:length(cnames))
{
    cnames[i] = substr(cnames[i],1,15)
}

colnames(GTEx_junctions) = cnames

print("Columns renamed")

#not all muscular samples present in junction table
#sample present in junctions: GTEX.11ZUS.1426
#sample absent in junctions: GTEX.S33H.2226
samples = samples[samples %in% colnames(GTEx_junctions)]

print("Samples in calculated")

junction_ids = GTEx_junctions$junction_id

GTEx_junctions = subset(GTEx_junctions,select = samples)

print("Junctions subsetted")

GTEx_counts = rowSums(GTEx_junctions)

GTEx = data.frame(row.names = junction_ids)
GTEx = cbind(GTEx,GTEx_counts)

GTEx = subset(GTEx,GTEx_counts > 4)

write.table(GTEx,paste0(args[2],".junction_counts.txt"))
