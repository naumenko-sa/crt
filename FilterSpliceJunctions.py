#!/usr/bin/python3

import os
import sys
import sqlite3
from AddJunctionsToDatabase import connectToDB, commitAndClose

def tableHeader():
	header = ['gene','pos', 'annotation', 'read_count', 'norm_read_count','n_gtex_seen','total_gtex_read_count']
	return (','.join(header) + '\n')

def countGTEX(cur):
	cur.execute('select count(*) from SAMPLE_REF where type = 0;') # 0 = GTEX, 1 = PATIENT
	return cur.fetchone()[0]

def countPatients(cur):

	"""
	Counts the number of patient samples in the database

	Args:
		cur, a cursor to a connection to a database

	Returns:
	    The number of patient files in the database

	Raises:
	    None
	"""

	cur.execute('select count(*) from SAMPLE_REF where type = 1;') # 0 = GTEX, 1 = PATIENT

	return cur.fetchone()[0]

def writeToFile(res, file):
	with open(file, "w") as out:
		out.write(tableHeader())
		for row in res:
			out.write(','.join(str(element) for element in row) + '\n')

def sampleSpecificJunctions(cur, sample, min_read, min_norm_read):

	"""
	Generates a file with junctions seen in a sample and not seen in any
	GTEx samples with a read count equal to or greater than the specified minimum read count. 

	Note that this function does not discriminate against junctions seen in other patient samples.
	Users should also note that the function does not work for GTEx samples as the database query
	relies on n_gtex_seen being 0.

	The query provides information about read counts of only the one specified sample
	in the columns 'sample:read_count' and 'sample:norm_read_count', 

		Ex. 1:100-200 PATIENT2:344

	This occurs because we are joining on only one sample in the database as opposed
	to all.

	Args:
		cur, a cursor to a connection to a database
		sample, the name of the sample file you want to investigate, must include .bam extension
		min_read, the minimum number of reads a junction must have
		min_norm_read, the minimum normalized read count a junction must have or NULL

	    not reporting junctions with NULL norm_count (mostly NONE annotated)
	"""

	count = str(countGTEX(cur))

	output = '_'.join([sample, 'specific', 'rc' + str(min_read), ('norm_rc' + str(min_norm_read)), 'n_gtex_' + count])

	cur.execute('''select gene_ref.gene,
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop),
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.total_patient_read_count,
		junction_counts.norm_read_count,
		junction_ref.n_gtex_seen,
		junction_ref.total_gtex_read_count
		from junction_counts, sample_ref, junction_ref, gene_ref
		where 
			sample_ref.sample_name = ? and
			junction_counts.read_count >= ? and
			junction_counts.norm_read_count >= ? and junction_counts.norm_read_count!='NULL' and
			junction_ref.n_gtex_seen <= 5 and
			sample_ref.rowid=junction_counts.bam_id and
			junction_counts.junction_id = junction_ref.rowid and
			junction_ref.rowid = gene_ref.junction_id;''',
		(sample, min_read, min_norm_read))

	writeToFile(cur.fetchall(), output)

def customSampleSpecificJunctions(cur, sample, min_read, min_norm_read, max_n_gtex_seen, max_total_gtex_reads):

	"""
	Generates a text file using a query in which you can discover junctions specific to a sample
	with parameters for appearances in GTEx samples

	The query provides information about read counts of only the one specified sample
	in the columns 'sample:read_count' and 'sample:norm_read_count', 

		Ex. 1:100-200 PATIENT2:344

	This occurs because we are joining on only one sample in the database as opposed
	to all.

	Args:
		cur, a cursor to a connection to a database
		sample, the name of the sample file you want to investigate, must include .bam extension
		min_read, the minimum number of reads a junction must have
		min_norm_read, the minimum normalized read count a junction must have or NULL
		max_n_gtex_seen, the maximum number of gtex samples a junction can appear in
		max_total_gtex_reads, the maximum total read count for a junction in GTEx samples

	Returns:
	    None

	Raises:
	    None
	"""

	if not max_n_gtex_seen:
		max_n_gtex_seen = 0

	if not max_total_gtex_reads:
		max_total_gtex_reads = 0

	if not min_read:
		min_read = 0

	output = '_'.join([str(sample), ('rc' + str(min_read)), ('norm_rc' + str(min_norm_read)), ('maxGTEX' + str(max_n_gtex_seen)), ('maxGTEXrc' + str(max_total_gtex_reads))])

	# does not report events with norm_count==NULL
	cur.execute('''select group_concat(gene_ref.gene),
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop),
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.n_gtex_seen, 
		junction_ref.n_patients_seen,
		junction_ref.total_patient_read_count,
		junction_ref.total_gtex_read_count,
		junction_ref.total_read_count,
		group_concat(distinct junction_counts.read_count||':'||sample_ref.sample_name),
		group_concat(distinct junction_counts.norm_read_count||':'||sample_ref.sample_name)
		from 
		sample_ref inner join junction_counts on sample_ref.ROWID = junction_counts.bam_id 
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		inner join gene_ref on junction_ref.rowid = gene_ref.junction_id 
		where
		sample_ref.sample_name = ? and
		junction_counts.read_count >= ? and
		junction_ref.n_gtex_seen <= ? and
		junction_ref.total_gtex_read_count <= ? and
		junction_counts.norm_read_count >= ?
		group by 
		junction_ref.chromosome, junction_ref.start, junction_ref.stop;''',
		(sample, min_read, max_n_gtex_seen, max_total_gtex_reads))

	writeToFile(cur.fetchall(), output)

def printSamplesInDB(cur):

	"""
	Prints to stdout all samples in the database and their experiment type

	Args:
		cur, a cursor to a connection to a database

	Returns:
	    None

	Raises:
	    None
	"""

	cur.execute('''select sample_name,
	case
		when type = 0 then 'CONTROL'
		when type = 1 then 'PATIENT'
	END
	from SAMPLE_REF;''')

	for line in cur.fetchall():
		print('\t'.join(str(i) for i in line))

def printAllJunctions(cur):

	"""
	Dumps all junction information seen in all samples to a text file.

	The query provides information about read counts of a junction across all samples
	unlike the other queries in the columns 'sample:read_count' and 'sample:norm_read_count', 

		Ex. 1:100-200 GTEx1:20,GTEx3:211,PATIENT2:344

	This occurs because we are joining and grouping all sample names in the database as opposed
	to just one name.

	Args:
		cur, a cursor to a connection to a database

	Returns:
	    None

	Raises:
	    None
	"""

	output = 'all_junctions_n_gtex_' + str(countGTEX(cur)) + '_n_paitents_' + str(countPatients(cur))

	cur.execute('''select group_concat(gene_ref.gene),
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop),
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.n_gtex_seen, 
		junction_ref.n_patients_seen,
		junction_ref.total_patient_read_count,
		junction_ref.total_gtex_read_count,
		junction_ref.total_read_count,
		group_concat(distinct junction_counts.read_count||':'||sample_ref.sample_name),
		group_concat(distinct junction_counts.norm_read_count||':'||sample_ref.sample_name)
		from 
		sample_ref inner join junction_counts on sample_ref.ROWID = junction_counts.bam_id 
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		inner join gene_ref on junction_ref.rowid = gene_ref.junction_id
		where junction_ref.total_read_count > 0
		group by 
		junction_ref.chromosome, junction_ref.start, junction_ref.stop;''')

	writeToFile(cur.fetchall(), output)

if __name__=="__main__":
	
	conn, cur = connectToDB()

	# sample = sys.argv[2]
	# min_read = int(sys.argv[3])
	# max_n_gtex_seen = int(sys.argv[4])
	# max_total_gtex_reads = int(sys.argv[5])

	if sys.argv[1] == '--printsamples':
		printSamplesInDB(cur)
	elif sys.argv[1] == '--sample':
		print(sys.argv[4])
		sampleSpecificJunctions(cur, sys.argv[2], int(sys.argv[3]), float(sys.argv[4]))
	elif sys.argv[1] == '--custom':
		customSampleSpecificJunctions(cur, sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6])
	elif sys.argv[1] == '--all':
		printAllJunctions(cur)
	else:
		print('Invalid option. Use one of the following:')
		print('--printsamples')
		print('--sample	[SAMPLE]	[MIN_READ]	[MIN_NORM_READ]')
		print('--custom [SAMPLE] [MIN_READ] [MIN_NORM_READ]	[MAX_N_GTEX_SEEN] [MAX_TOTAL_GTEX_READS]')
		print('--all')

	commitAndClose(conn)