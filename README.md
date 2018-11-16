# crt

1. Run bcbio with crt.bcbio.rnaseq.yaml.
  1. create project/input
  2. name sample SX_case-N-tissue
  3. crt.prepare_bcbio_run.sh

Don't trim reads to save all data and delete fastq files.

2. Look at DNA variants
  1. cre.vcf2cre.sh
  2. ```qsub ~/cre/cre.sh -v family=project,type=rnaseq```

3. Look at gene expression
  1. crt.bam2rpkm.sh - counts for RPKM calculation in R
  2. crt.load_rpkm_counts.R - load counts into R
  3. crt.muscular.R - functions for muscular project

4. Find pathogenic splice events (Modification of [MendelianRNA-seq](https://github.com/berylc/MendelianRNA-seq))

  1. Overview

- SpliceJunctionDiscovery.py discovers splice junctions calling samtools and parsing CIGARs.
- AddJunctionsToDatabase.py load junction information to the sqlite database.
- FilterSpliceJunctions.py outputs junctions present in a sample and absent in controls.
* [The original repository of Beryl Cummings](https://github.com/berylc/MendelianRNA-seq)
* [Modification by Dennis Kao](https://github.com/dennis-kao/MendelianRNA-seq-DB)
* [Article: Cummings et al. 2017](http://stm.sciencemag.org/content/9/386/eaal5209) 
* [Manual](https://macarthurlab.org/2017/05/31/improving-genetic-diagnosis-in-mendelian-disease-with-transcriptome-sequencing-a-walk-through/)

  2. Steps
	1. Discover junctions, submit a torque job:

- `qsub [path-to-crt]/crt.splice_junction_discovery.pbs -v bam=file.bam`
- or `qsub [path-to-crt]/crt.splice_junction_discovery.pbs -v genes=my_gene_panel.bed,bam=file.bam`
- or `python3 [path-to-crt]/SpliceJunctionDiscovery.py -bam=file.bam`

Result: file.bam.junctions.txt

	2. Load GENCODE junctions to the database
```
python3 [path-to-crt]/AddJunctionsToDatabase.py \
	--addGencode \
	-transcript_model=[path-to-crt]/gencode.comprehensive.splice.junctions.txt
```
result: SpliceJunction.db

	3. Load junctions from samples to the SpliceJunctions.db database (load controls once, and copy SpliceJunctions.db for every analysis).

- `qsub [path-to-crt]/crt.load_junctions.pbs -v bam=file.bam`
- or `qsub [path2crt]/crt.load_junctions.localhd.pbs -v bam=file.bam`
- or `python3 [path-to-crt]/analysis/AddJunctionsToDatabase1.py -addBAM -transcript_file [path-to-MendelianRNA-seq-db]/all-protein-coding-genes-no-patches.txt -bamlist bamlist.list -flank 1`
- If flank was set to 1, a gencode junction was 1:100-300 and a junction in a sample was 1:99-301, the sample junction would be considered BOTH annotated. This is because both the start and stop positions fall within a +/- 1 range of that of a transcript_model's junction.

result: junctions from a bam file (files) added to SpliceJunctions.db

	4. Filter out junctions present in GTEx controls, report rare junctions in a sample

- [path2crt]/crt.filter_junctions.sh file.bam
- or 
``` 
python3 [path2crt]/FilterSpliceJunctions.py \
    --sample [sample.bam] \
    [MIN_READ_COUNT]	\
    [MIN_NORM_READ_COUNT]
```

Parameters:
- bam extension is necessary to include
- [MIN_READ_COUNT] = 5
- [MIN_NORM_READ_COUNT] = 0.05
	
With ```--sample``` option columns ```sample:read_count``` and ```sample:norm_read_count``` will not show read counts from other samples. 
For all samples use ```---all``` option.

##4.3  Output

By default the database is named SpliceJunction.db. There are 4 tables:

	1. SAMPLE_REF, a list of samples and their type (0 = GTEX or control, 1 = patient)
	2. JUNCTION_REF, a list of junctions and their frequency of appearances in samples
	3. JUNCTION_COUNTS, read counts of junctions in a sample
	4. GENE_REF, an annotation of junctions with genes, a single junction can map to multiple genes

Using one of the options of FilterSpliceJunctions.py will produce a text file containing junction information in the following format:

	gene	chromosome:start-stop	annotation	n_gtex_seen	n_patients_seen	total_patient_read_count	total_gtex_read_count	total_read_count	sample:read_count	sample:norm_read_count
	MT-ND1	MT:3540-3610	NONE	0	1	11	0	11	PATIENT.bam:11	PATIENT.bam:NULL
	AC002321.2	GL000201.1:4130-9415	NONE	1	1	32	4	36	PATIENT.bam:32	PATIENT.bam:NULL
	MT-CO1	MT:7276-13822	NONE	1	1	5	1	6	PATIENT.bam:5	PATIENT.bam:NULL
	MT-ATP6	MT:9234-9511	NONE	0	1	6	0	6	PATIENT.bam:6	PATIENT.bam:NULL
	AC002321.2	GL000201.1:9511-14322	START	1	1	70	2	72	PATIENT.bam:70	PATIENT.bam:NULL


## Dependencies

1. .bam (and .bai) files of cases and GTEx controls.
2. genes - A bed file with gene coordinates
	```
	CHROM	START	STOP	GENE
	```
	Use genes.bed - a list of protein coding genes, or your own bed file. To retrieve gene coordinates from ENSMBL biomart use 
	[genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) 
	```
3. transcript_model - text file containing a list of known canonical splice junctions. [included file](gencode.comprehensive.splice.junctions.txt). This file contains junctions from gencode v19.

4. [Python 3.5.2](https://www.python.org/downloads/) or higher

5. Python [CIGAR string library](https://pypi.python.org/pypi/cigar/0.1.3) by Brent Pedersen

6. sqlite3 Python library based off of SQLite3 version 3.11.0 or higher. You can check your library's version with:
	```
	import sqlite3
	print (sqlite3.sqlite_version_info)
	```
	


## Differences from Beryl Cumming's original MendelianRNA-seq

- junctions are stored in a database. We process GTEx controls once and reuse the database.
![alt text](./SpliceJunctionSchema.png)
- SpliceJunctionDiscovery has been rewritten in Python
- CIGAR string parsing is handled by a function called parseCIGARForIntrons() whereas before CIGAR strings were handled by piping through multiple bash tools. As a result of improper parsing using bash tools, junction start and/or stop positions were not reported properly (e.x. 1:100-200*1D30 represents an alignment that should really be 1:100-230 or 1:100-231)
- Transcript_model annotation and flanking have been implemented using database logic
- All information produced by SpliceJunctionDiscovery is stored in a database instead of text files. This allows the user to utilize previously computed results instead of having to run the entire pipeline again when a new sample needs to be analyzed.
- Transcript_model annotation now discriminates between 'START' and 'STOP' instead of 'ONE'. In addition, there is a new annotation, called 'EXON_SKIP' which denotes the event of exon skipping. This is done by checking to see if the reported 3' and 5' positions from a sample's junction belong to different transcript_model junctions.
- Normalization of annotated junctions now considers read counts from all junctions that have at least one annotated splice site as the denominator whereas before only "BOTH" annotated junctions were used

## Footnotes

### A junction in multiple gene regions
A single gene region can encompass partial or whole regions of other genes. Thus, the same junction can appear in 2 different gene text files in a sample folder generated by SpliceJunctionDiscovery. Whether a junction belongs to two or more genes is not always factually correct. However, for the sake of inclusion and for the fact that this rarely happens, this edge case has been accounted for in AddJunctionsToDatabase.py in two ways: 

	1. The mapping of a single junction to multiple genes has been done with the table GENE_REF
	2. If the script encounters the same junction in a sample more than once, it will utilize the result with the highest read count for read count and normalized read count and will discard the other.

### Splice site flanks and annotation
A +/- flanking region is considered when annotating the 5' and 3' positions of sample junctions to increase the number of annotated junctions. This value is specified by the -flank parameter (default 1). There is an option to not use flanking at all (-flank 0).

