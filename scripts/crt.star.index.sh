#PBS -l walltime=10:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

STAR --genomeDir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Mmusculus/mm10/star \
--genomeFastaFiles /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Mmusculus/mm10/seq/mm10.fa --runThreadN 10 --limitGenomeGenerateRAM 30000000000 --genomeChrBinNbits 14 --runMode genomeGenerate --genomeSAindexNbases 14