# Validation of small variant calling with RNA-seq data

## Data
RNA-Seq of GM12878, SRR307898.https://www.ncbi.nlm.nih.gov/sra/?term=SRR307898. (SRR307897 is of bad quality).
This dataset was used in [Piskol2013](https://www.ncbi.nlm.nih.gov/pubmed/24075185) article. Quite old - Illumina GAII.
Better NA12878 RNA-seq?

## bcbio NA12878.yaml config
```
details:
- algorithm:
    aligner: star
    strandedness: unstranded
    variantcaller: gatk-haplotype
    tools_off:
    - gatk4
  analysis: RNA-seq
  description: NA12878_SRR307898
  files:
  - /path/NA12878/input/NA12878_SRR307898_1.fq.gz
  - /path/NA12878/input/NA12878_SRR307898_2.fq.gz
  genome_build: GRCh37
  metadata:
    batch: NA12878
fc_name: NA12878
resources:
  default:
    cores: 2
    jvm_opts:
    - -Xms750m
    - -Xmx7000m
    memory: 15G
upload:
  dir: ../final
```

## Validation

Uses two files from [cre](https://github.com/naumenko-sa/cre), required bcbio installation.

```
bcftools view -f PASS NA12878-gatk-haplotype-annotated.vcf.gz | grep -v possible_rnaedit | bgzip -c  > NA12878.pass.vcf.gz
tabix NA12878.pass.vcf.gz
cre.rtg.validate.sh NA12878.pass.vcf.gz ~/cre/data/intersect.bed
```

## Results

|date|type|bcbio|gatk|TP|FP|FDR|
|-|-|-|-|-|-|-|
|2016-12-08|SNP|0.9.9 or 1.0.0|3.6-24|6655|160|2%|
|2016-12-08|INDEL|0.9.9 or 1.0.0|3.6-24|81|75|48%|
|2018-07-14|SNP|1.0.9a0|3.8|6658|89|1%|
|2018-07-14|INDEL|1.0.9a0|3.8|89|168|65%|
|2018-07-14|SNP|1.0.9a0|4.0.3.0|5835|1022|15%|
|2018-07-14|INDEL|1.0.9a0|4.0.3.0|83|12439|99%|
|2018-04-02|SNP|1.1.3|				
|2018-04-02|INDEL|1.1.3|
