# Loads a file with counts from feature_counts
# loads gene lengths
# loads gene names
# calculates RPKMs
# returns ENS_ID, rpkm, Gene_name
load_rpkm_counts = function(filename)
{
    #test:
    #filename="/home/sergey/Desktop/S01_1-1-B.bam.rpkm_counts.txt"
    library(edgeR)   
    ensembl_w_description = read.delim2("~/cre/ensembl_w_description.txt", row.names=1, stringsAsFactors=F)
    #first line in the file is a comment
    counts = read.delim(filename, stringsAsFactors=F, row.names=1,skip=1)
    counts$Chr=NULL
    counts$Start=NULL
    counts$End=NULL
    counts$Strand=NULL
    
    
    Gene_lengths = counts$Length
    
    counts$Length=NULL
    
    counts = rpkm(counts,Gene_lengths)
    
    counts = merge(counts,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(counts)=counts$Row.names
    counts$Row.names=NULL
    counts$Gene_description=NULL
    
    colnames(counts) = gsub(".bam","",colnames(counts))
    
    return(counts)
}

load_raw_counts = function(filename)
{
    #first line in the file is a comment
    counts = read.delim(filename, stringsAsFactors=F, row.names=1,skip=1)
    counts$Chr=NULL
    counts$Start=NULL
    counts$End=NULL
    counts$Strand=NULL
    counts$Length=NULL
  
    counts = merge(counts,ensembl_w_description,by.x="row.names",by.y="row.names")
    row.names(counts)=counts$Row.names
    counts$Row.names=NULL
    counts$Gene_description=NULL
  
    colnames(counts) = gsub(".bam","",colnames(counts))
  
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

merge_counts = function()
{
    samples = read.table("samples.txt", quote="\"", stringsAsFactors=F)
    samples = samples[,1]
  
    for (sample in samples)
    {
        table_name = sample
        file_name = paste0(sample,".counts")
        assign(table_name,read.delim(file_name, header=F,row.names=1,stringsAsFactors=F))
        table_itself = get(table_name)
        colnames(table_itself) = c(sample)
        assign(table_name,table_itself)
    }
  
    samples.data = get(samples[1])
  
    for (sample in tail(samples,-1))
    {
      #test
      #table_name = samples[2]
        table_name = sample
        table_itself = get(table_name)
        samples.data = merge(samples.data,get(table_name),by.x="row.names",by.y="row.names")
        row.names(samples.data)=samples.data$Row.names
        samples.data$Row.names = NULL
    }
  
    write.table(samples.data,"all_counts.txt",quote=F)
}
