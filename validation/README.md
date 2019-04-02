# Validation of small variant calling with RNA-seq data.

## Data
 RNA-Seq of GM12878, SRR307898.https://www.ncbi.nlm.nih.gov/sra/?term=SRR307898. (SRR307897 is a bad quality data).
 
## [Piskol2013](https://www.ncbi.nlm.nih.gov/pubmed/24075185) article.

## Results

|Metrics| gatk3.6, Nov 2016 |gatk 4.0.3.0|gatk3.8|
|-|-|---|---|
| snp tp-baseline      | 6655|5835|6658
| snp fp 186             | 160 |1022|89
| SNP FDR | 2%|15%|1.3%
| indels tp-baseline |81|83|89
| indels fp |75|12439|168
| Indels FDR |48%|99%|65%