#!/usr/bin/python3
#!/usr/bin/env bash

import os
import sys
import argparse
import multiprocessing
import subprocess
import sqlite3
import re
from datetime import datetime

def connectToDB():
	conn = sqlite3.connect('SpliceJunction.db', timeout=80)
	cur = conn.cursor()
	return conn, cur

def commitAndClose(conn):
	conn.commit()
	conn.close()

def initializeDB():

	"""
	Sets up the tables and SQLite settings required for the splice junction database

	SAMPLE_REF: Contains all BAM file names and their experiment type
		type = {0, 1} 
		GTEX, patient

	JUNCTION_REF: A collection of all junction positions seen, their transcript_model
	annotation and basic population read count statistics
		gencode_annotation = {0, 1, 2, 3, 4}
		none, only start, only stop, both, exon skipping

	JUNCTION_COUNTS: The individual read counts of junctions pertaining to a sample
	and their calculated sample specfic normalized read counts

	GENE_REF: A mapping of junction positions to certain genes. Sometimes gene positions
	can encompass multiple other smaller genes

	TRANSCRIPT_MODEL_JUNCTIONS: A storage of junctions from the user specific transcript_model
	parameter. This is table is only used as a reference for gencode annotation and normalization
	and has no actual relevance to the other tables.

	"""

	conn, cur = connectToDB()
	cur.execute('''PRAGMA foreign_keys = ON;''')
	cur.execute('''create table if not exists SAMPLE_REF (
		sample_name varchar(50) primary key, 
		type tinyint not null);''') 

		# type = {0, 1} 
		# GTEX, patient

	cur.execute('''create table if not exists JUNCTION_REF (
		chromosome tinyint not null,
		start unsigned big int not null,
		stop unsigned big int not null,
		gencode_annotation tinyint not null,
		n_patients_seen unsigned big int default 0,
		n_gtex_seen unsigned big int default 0,
		total_patient_read_count big int default 0,
		total_gtex_read_count big int default 0,
		total_read_count big int default 0,
		primary key (start, stop, chromosome));''') 
		
		# gencode_annotation = {0, 1, 2, 3, 4}
		# none, only start, only stop, both, exon skipping

	cur.execute('''create table if not exists JUNCTION_COUNTS (
		bam_id integer not null,
		junction_id integer not null,
		read_count unsigned big int not null,
		norm_read_count float,
		foreign key(bam_id) references SAMPLE_REF(ROWID),
		foreign key(junction_id) references JUNCTION_REF(ROWID),
		primary key (junction_id, bam_id));''')

	cur.execute('''create table if not exists GENE_REF (
		gene varchar(30) not null,
		junction_id integer not null,
		foreign key(junction_id) references JUNCTION_REF(ROWID),
		primary key (junction_id, gene));''')

	cur.execute('''create table if not exists TRANSCRIPT_MODEL_JUNCTIONS (
		chromosome tinyint not null,
		start unsigned big int not null,
		stop unsigned big int not null,
		primary key (chromosome, start, stop));''')

	cur.execute('''create index if not exists stopJunction
		on TRANSCRIPT_MODEL_JUNCTIONS (chromosome, stop);
		''')

	commitAndClose(conn)

def annotateJunction(cur, chrom, start, end, flank):
	"""
		Annotates a junction with gencode
	"""
	if flank > 0:
		cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
				chromosome is ? and
				start >= ? and
				start <= ? and
				stop >= ? and
				stop <= ?;''', (chrom, (start - flank), (start + flank), (end - flank), (end + flank)) )
		isBothAnnotated = cur.fetchone()

		if not isBothAnnotated:
			cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
					chromosome = ? and
					start >= ? and
					start <= ?;''', (chrom, (start - flank), (start + flank)) )
			isStartAnnotated = cur.fetchone()

			cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
					chromosome = ? and
					stop >= ? and
					stop <= ?;''', (chrom, (end - flank), (end + flank)))
			isStopAnnotated = cur.fetchone()
	else:
		cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
				chromosome = ? and
				start = ?
				stop is = ?;''', (chrom, start, end) )
		isBothAnnotated = cur.fetchone()

		if not isBothAnnotated:
			cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
					chromosome = ? and
					start = ?;''', (chrom, start))
			isStartAnnotated = cur.fetchone()

			cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
					chromosome = ? and
					stop = ?;''', (chrom, end))
			isStopAnnotated = cur.fetchone()

	if isBothAnnotated:
		annotation = 3 # both annotated
	elif isStopAnnotated and isStartAnnotated:
		annotation = 4 # exon skipping
	elif isStopAnnotated:
		annotation = 2 # only stop
	elif isStartAnnotated:
		annotation = 1 # only start
	else:
		annotation = 0 # novel junction

	return annotation

def getJunctionID(cur, chrom, start, end, flank):

	"""
	Retrieves the ROWID and annotation of a junction from the database

	If the junction does not exist in the database, then the function adds it
	and returns the appropriate values

	Args:
		cur, the cursor of a database connection
		chrom, the chromosome a junction lies on
		start, the 5' splice site of a junction
		stop, the 3' splice site of a junction
		flank, the +/- range a junction's start and stop site must fall within
		the transcript_model's start and stop site to be considered "annotated" 

	Returns:
	    ROWID (junction_id), gencode_annotation of a junction
	"""
	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only stop, both, exon skipping
	# thus, gencode junctions will always have a gencode_annotation value of 3

	# check if junction exists in the database (probably from GTEx controls or another sample)
	# if junction exists, but the annotation is different, overwrite the annotation
	
	# reannotate junction
	new_annotation = annotateJunction(cur, chrom, start, end, flank)
	
	cur.execute('''select ROWID, gencode_annotation from JUNCTION_REF where 
		chromosome is ? and
		start is ? and 
		stop is ?;''', (chrom, start, end))
	res = cur.fetchone()

	if res:
		ROWID, annotation = res
		if new_annotation != annotation:
			cur.execute('''update JUNCTION_REF 
				set gencode_annotation = ?
				where 
					chromosome is ? and
					start is ? and 
					stop is ?;''', (new_annotation,chrom, start, end))
	else:
		cur.execute('''insert into JUNCTION_REF (
				chromosome, 
				start, 
				stop, 
				gencode_annotation) 
				values (?, ?, ?, ?);''', (chrom, start, end, new_annotation))

		ROWID = cur.lastrowid
		
		
	return ROWID, new_annotation

def makeSpliceDict(junction_file,gene):
	spliceDict = {}
	with open(junction_file, "r") as gf:
		for line in gf:
			chrom, start, stop, count, j_gene = line.strip().split()
			if j_gene == gene:
			    uniqueSplice = (chrom, start, stop) 
			    spliceDict[uniqueSplice] = int(count)
	return spliceDict

def normalizeReadCount(spliceDict, junction, annotation, max_counts):
	"""
	Normalizes the read count of a splice site
	One splice site (junction start or end) is annotated, then normalize that site
	Neither site is annotated, don't perform normalization
	If both sites are annotated or there is a case of exon skipping, perform normalization
	on the site which has the largest read count.

	Args:
		gene_file, the path to a text file produced by SpliceJunctionDiscovery.py

	Returns:
	    norm_count, the ratio between an annotated splice site's read count to that of
	    the annotated splice site belonging to a junction with the largest read count 

	"""

	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only end, both, exon skipping

	chrom, start, stop = junction
	annotation = int(annotation)

	startString = makeStartString(chrom, start)
	stopString = makeStopString(chrom, stop)

	if annotation == 0: 
		return 'NULL'
	elif annotation == 1:
		key = startString
	elif annotation == 2:
		key = stopString
	elif (annotation == 3) or (annotation == 4):
		# use the bigger read count for normalization
		if max_counts[startString] > max_counts[stopString]:
			key = startString
		else:
			key = stopString
	norm_count = round((float(spliceDict[junction]) / float(max_counts[key])), 3)
	return str(norm_count)

def makeStartString(chrom, start):
	return ''.join([chrom,':','START',':',start])

def makeStopString(chrom, stop):
	return ''.join([chrom,':','STOP',':',stop])

def getMaxCounts(spliceDict):
	"""max read counts for normalization"""
	count_dict = {}
	for junction in spliceDict:
		chrom, start, stop = junction
		startString = makeStartString(chrom, start)
		stopString = makeStopString(chrom, stop)

		if startString in count_dict:
			if count_dict[startString] < spliceDict[junction]:
				count_dict[startString] = spliceDict[junction]
		else:
			count_dict[startString] = spliceDict[junction]

		if stopString in count_dict:
			if count_dict[stopString] < spliceDict[junction]:
				count_dict[stopString] = spliceDict[junction]
		else:
			count_dict[stopString] = spliceDict[junction]
	return count_dict

def addJunctionsForAGene(bam,gene,flank,bam_id,bam_type,sample):
	conn, cur = connectToDB()
	junction_file = bam+'.junctions.txt'
	spliceDict = makeSpliceDict("gene_buf.txt",gene)
	max_counts = getMaxCounts(spliceDict)
	for junction in spliceDict:
		chrom, start, end = junction
		reads = spliceDict[junction]
		junction_id, annotation = getJunctionID(cur, chrom, int(start), int(end), flank)
		try:
			norm_read_count = normalizeReadCount(spliceDict, junction, annotation, max_counts)
		except ZeroDivisionError:
			print("Zero division error when normalizing %s:%s-%s in genefile %s.txt in sample %s with annotation %d"%(chrom, start, end, gene, sample, annotation))
			norm_read_count = 'null'
		annotateJunctionWithGene(gene, junction_id, cur)
		updateJunctionInformation(junction_id, bam_id, bam_type, gene, sample, reads, norm_read_count, cur)
	del spliceDict, max_counts
	commitAndClose(conn)

def updateJunctionInformation(junction_id, bam_id, bam_type, gene, sample, new_read_count, new_norm_read_count, cur):
	"""
	Adds a junction's position and its read counts to the database. Logic for total read counts has also been implemented.

	In the case that a junction already exists for a sample, the larger read count and its corresponding normalized read count
	is used.

	Args:
		junction_id, the ROWID of a junction in JUNCTION_REF
		bam_id, the ROWID of a sample in SAMPLE_REF
		bam_type, the type of the sample in the experiment (control or disease?)

		gene, the name of the gene in which the junction position falls under. This value is the name of the text file generated
		by SpliceJunctionDiscovery.py. Because gene regions can encompass that of other genes, the program will run in to the case
		where a single junction shows up in 2 gene files. In that case, the new gene is simply added to GENE_REF.
		sample, the name of the BAM file being processed
		
		new_read_count, the reported read count from the text file
		new_norm_read_count, the calculated normalized read count from the function normalizeReadCount()
		cur, a cursor to a connection to the database

	Returns:
	    None

	Raises:
	    None
	"""
	# check if sample already has the junction in the database
	cur.execute('''select ROWID, read_count from JUNCTION_COUNTS where junction_id is ? and bam_id is ?;''', (junction_id, bam_id))
	res = cur.fetchone()

	# if it is, check if new_reads > old_reads, update JUNCTION_REF and JUNCTION_COUNTS for the appropriate sample
	if res:
		sample_junction_id, old_read_count = res

		if int(new_read_count) > int(old_read_count):

			# update entry to reflect new read count values
			cur.execute('''update JUNCTION_COUNTS set read_count = ?, norm_read_count = ? where ROWID = ?;''', (new_read_count, new_norm_read_count, sample_junction_id))

			# update total read counts
			cur.execute('''update JUNCTION_REF set total_read_count = total_read_count - ? + ? where ROWID = ?;''', (old_read_count, new_read_count, junction_id))

			# update total patient or gtex read counts
			if bam_type == 1:
				cur.execute('''update JUNCTION_REF set total_patient_read_count = total_patient_read_count - ? + ? where ROWID = ?;''', (old_read_count, new_read_count, junction_id))
			elif bam_type == 0:
				cur.execute('''update JUNCTION_REF set total_gtex_read_count = total_gtex_read_count - ? + ? where ROWID = ?;''', (old_read_count, new_read_count, junction_id))

	# if not, add it with read counts and normalized read counts, increment n_samples_seen, increment n_times_seen, increment JUNCTION_REF total times seen
	else:
		cur.execute('''insert into JUNCTION_COUNTS (bam_id, junction_id, read_count, norm_read_count) values (?, ?, ?, ?);''', (bam_id, junction_id, new_read_count, new_norm_read_count))

		# 0 = gtex, 1 = patient
		if bam_type == 1:
			cur.execute('''update JUNCTION_REF set 
				n_patients_seen = n_patients_seen + 1, 
				total_read_count = total_read_count + ?, 
				total_patient_read_count = total_patient_read_count + ? 
				where ROWID = ?;''', (new_read_count, new_read_count, junction_id))
		elif bam_type == 0:
			cur.execute('''update JUNCTION_REF set 
				n_gtex_seen = n_gtex_seen + 1, 
				total_read_count = total_read_count + ?, 
				total_gtex_read_count = total_gtex_read_count + ? 
				where ROWID = ?;''', (new_read_count, new_read_count, junction_id))

def get_bam_id_and_type(cur, bam):
	"""
	Gets the ROWID and the experiment type of a sample in the database. 
	This function works on the assumption that the sample already exists in the database.
	Before running this function, parallel_process_gene_files will have already added bam file names to the database,
	hence satisfying this assumption

	Args:
		cur, a cursor to a connection to the database
		bam, the name of a bam file in the database

	Returns:
	    bam_id, the ROWID of a sample in SAMPLE_REF
	    bam_type, 0 or 1, a number which indicates whether a sample is control or a patient

	Raises:
	    None
	"""
	cur.execute('''select ROWID, type from SAMPLE_REF where sample_name = ?;''', (bam, ))
	bam_id, bam_type = cur.fetchone()
	return bam_id, bam_type

def addSampleToDatabase(bam):
	conn, cur = connectToDB()
	if 'GTEX' in bam:
		cur.execute('''insert into SAMPLE_REF (sample_name, type) values (?, 0);''', (bam, ))
	else: # sample is a patient
		cur.execute('''insert into SAMPLE_REF (sample_name, type) values (?, 1);''', (bam, ))
	commitAndClose(conn)

def getGeneNames(bam):
	geneDict = {}
	geneList = []
	with open(bam+".junctions.txt", "r") as gf:
		for line in gf:
			gene = line.strip().split()[4]
			if not gene in geneDict:
			    geneDict[gene]=1
			    geneList.append(gene)
	return geneList

def addJunctions(bam, genes, flank):
	flank = int(flank)
	addSampleToDatabase(bam)
	geneList = getGeneNames(bam)
	conn, cur = connectToDB()
	bam_id, bam_type = get_bam_id_and_type(cur, bam)
	commitAndClose(conn)
	sample = bam[:-4]
	prev_gene=""
	gene_file = open("gene_buf.txt","w")
	with open(bam+".junctions.txt","r") as f:
		for line in f:
			gene = line.strip().split()[4]
			if (gene != prev_gene):
				if (prev_gene!= ""):
					gene_file.close()
					addJunctionsForAGene(bam,prev_gene,flank,bam_id,bam_type,sample)
					gene_file = open("gene_buf.txt","w")
			gene_file.write(line)
			prev_gene=gene
	gene_file.close()
	addJunctionsForAGene(bam,prev_gene,flank,bam_id,bam_type,sample)

def annotateJunctionWithGene(gene, junction_id, cur):
	
	"""
	Maps a junction to a gene in the database.

	Args:
		gene, the name of a gene
		junction_id, the ROWID of a junction in JUNCTION_REF
		cur, a cursor to a connection to the database

	Returns:
	    None

	Raises:
	    None
	"""
	cur.execute('''insert or ignore into GENE_REF (gene, junction_id) values (?, ?);''', (gene, junction_id))

def addTranscriptModelJunction(chrom, start, stop, cur):

	"""
	Adds a single junction from the transcript_model to the database's reference table, TRANSCRIPT_MODEL_JUNCTIONS

	Args:
		chrom, the chromosome a junction lies on
		start, the 5' splice site of a junction
		stop, the 3' splice site of a junction
		cur, a cursor to a connection to the database

	Returns:
	    None

	Raises:
	    None
	"""

	cur.execute('''insert or ignore into TRANSCRIPT_MODEL_JUNCTIONS (chromosome, start, stop) values (?, ?, ?);''', (chrom, start, stop))

def storeTranscriptModelJunctions(gencode_file):

	"""
	Adds junctions from a transcript_model to the database as a reference for annotation

	Args:
		gencode_file, a transcript_model containing known canonical junctions and their positions

	Returns:
	    None

	Raises:
	    None
	"""

	conn, cur = connectToDB()

	print ('Started adding transcript_model junctions @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	with open(gencode_file, "r") as gf:
		for line in gf:

			chrom, start, stop, gene = line.strip().split()[0:4]

			start = int(start)
			stop = int(stop)

			addTranscriptModelJunction(chrom, start, stop, cur)

	commitAndClose(conn)

	print ('Finished adding gencode annotations @ ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

def deleteSample(sample):

	"""
	Removes a sample and its read count information from the database.

	Does not remove records from JUNCTION_REF and so a juntion in JUNCTION_REF can have a read count of 0.

	Args:
		cur, a cursor to a connection to a database
		sample, the name of the sample file you want to remove, must include .bam extension
	
	Returns:
	    None

	Raises:
	    None
	"""

	conn, cur = connectToDB()

	cur.execute('select ROWID, type from sample_ref where sample_name = ?;', (sample, ))
	res = cur.fetchone()

	if not res:
		print ("Sample %s does not exist in the database!" % sample)
		exit(1)
	else:
		bam_id, bam_type = res

	cur.execute('select junction_id, read_count from JUNCTION_COUNTS where bam_id = ?;', (bam_id, ))

	for junction_id, read_count in cur.fetchall():

		if bam_type == 0:
			cur.execute('''update JUNCTION_REF set 
				n_gtex_seen = n_gtex_seen - 1,
				total_read_count = total_read_count - ?,
				total_gtex_read_count = total_gtex_read_count - ?
				where ROWID = ?;''', (read_count, read_count, junction_id))
		elif bam_type == 1:
			cur.execute('''update JUNCTION_REF set 
				n_patients_seen = n_patients_seen - 1,
				total_read_count = total_read_count - ?,
				total_patient_read_count = total_patient_read_count - ?
				where ROWID = ?;''', (read_count, read_count, junction_id))
		else:
			raise Exception ('FATAL ERROR - bam_id is not 0 or 1')

	cur.execute('''delete from JUNCTION_COUNTS where bam_id = ?;''', (bam_id, ))
	cur.execute('''delete from SAMPLE_REF where sample_name = ?;''', (sample, ))

	commitAndClose(conn)

	print ("Successfully deleted %s from database!" % sample)

if __name__=="__main__":

	print ('AddJunctionsToDatabase.py started on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	parser = argparse.ArgumentParser(description = 'Summarize the read counts of the junctions reported by SpliceJunctionDiscovery.py')
	parser.add_argument('-transcript_model',
						 help="Transcript model of canonical splicing, e.g. gencode v19. Default is set to [crt-home]/gencode.comprehensive.splice.junctions.txt",
						 action='store',
			    		default = "/home/naumenko/crt/gencode.comprehensive.splice.junctions.txt")
	parser.add_argument('-genes',
						help="The same transcript_file used in SpliceJunctionDiscovery.py",
						action='store',
			    		default = "/home/naumenko/crt/genes.bed")
	parser.add_argument('-bam',help='A bam file')
	parser.add_argument('-flank',
						help='Add a +/- flanking region for gencode annotation. Specify 0 if you don\'t want to use this feature, default=1',
						default=1)
	parser.add_argument('-sample',help='to be used with --delete, the name of the sample you want to remove from the database')
	# parser.add_argument('-db',help='The name of the database you are storing junction information in, default=SpliceJunction.db',default='SpliceJunction.db')

	mode_arguments = parser.add_mutually_exclusive_group(required=True)
	mode_arguments.add_argument('--addGencode',action='store_true',help='Populate the database with gencode junctions, this step needs to be done once before anything else')
	mode_arguments.add_argument('--addBAM',action='store_true',help='Add junction information from bamfiles found in the file bamlist.list')
	mode_arguments.add_argument('--delete',action='store_true',help='Delete a sample and its read counts from the database')
	args=parser.parse_args()

	print ('Working in directory ' + str(os.getcwd()))

	# databasePath = args.db

	initializeDB()

	if args.addGencode:
		print ('Storing junctions from the transcript model file ' + args.transcript_model)
		storeTranscriptModelJunctions(args.transcript_model)
	elif args.addBAM:
		print ('Storing junctions from bam files found in the file ' + args.bam)
		addJunctions(args.bam, args.genes, args.flank)
	elif args.delete:
		sample = args.sample

		if not sample:
			print('Please enter a sample name with its .bam extension using the parameter \'-sample SAMPLE_NAME\'')
			exit(1)

		deleteSample(sample)

	print ('AddJunctionsToDatabase.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
