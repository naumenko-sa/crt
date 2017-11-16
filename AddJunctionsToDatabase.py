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

# databasePath = ""

def connectToDB():

	"""
	Establishes a single unique connection to the database SpliceJunction.db
	and a cursor

	Args:
	    None

	Returns:
	    The connection and a cursor to that connection

	Raises:
	    None
	"""

	conn = sqlite3.connect('SpliceJunction.db', timeout=80)
	cur = conn.cursor()

	return conn, cur

def commitAndClose(conn):

	"""
	Commits any changes and closes a database connection
	
	Args:
	    conn: A single connection to a database

	Returns:
	    None

	Raises:
	    None
	"""

	conn.commit()
	conn.close()

def initializeDB():

	"""
	Sets up the tables and SQLite settings required for the splice junction database

	WAL mode - Write ahead logging, allows for reading while another worker process is
	writing to the database. Needed for concurrency reasons.

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

	Args:
		None

	Returns:
	    None

	Raises:
	    None
	"""

	conn, cur = connectToDB()

	# WAL mode is only present in SQLite versions 3.7.0 or above
	# Make sure your sqlite3 Python library is based off of the same SQLite3 or higher!
	# It is critical to having multiple writers and readers 
	cur.execute('''PRAGMA journal_mode = WAL;''') # WAL - write ahead logging - allows reading while writing. Needed for concurrency
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

def getJunctionID(cur, chrom, start, stop, flank):

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

	Raises:
	    None
	"""

	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only stop, both, exon skipping
	# thus, gencode junctions will always have a gencode_annotation value of 3

	# check if start and stop are apart of an existing gencode annotation
	cur.execute('''select ROWID, gencode_annotation from JUNCTION_REF where 
		chromosome is ? and
		start is ? and 
		stop is ?;''', (chrom, start, stop))
	res = cur.fetchone()

	if res:
		ROWID, annotation = res

	# if no such junction determine annotation of new junction: novel junction, only one annotated or a case of exon skipping?
	else:

		if flank > 0:
			cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
				chromosome is ? and
				start >= ? and
				start <= ? and
				stop >= ? and
				stop <= ?;''', (chrom, (start - flank), (start + flank), (stop - flank), (stop + flank)) )
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
					stop <= ?;''', (chrom, (stop - flank), (stop + flank)))
				isStopAnnotated = cur.fetchone()

		else:
			cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
				chromosome = ? and
				start = ?
				stop is = ?;''', (chrom, start, stop) )
			isBothAnnotated = cur.fetchone()

			if not isBothAnnotated:
				cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
					chromosome = ? and
					start = ?;''', (chrom, start))
				isStartAnnotated = cur.fetchone()

				cur.execute('''select * from TRANSCRIPT_MODEL_JUNCTIONS where 
					chromosome = ? and
					stop = ?;''', (chrom, stop))
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

		try:
			cur.execute('''insert into JUNCTION_REF (
				chromosome, 
				start, 
				stop, 
				gencode_annotation) 
				values (?, ?, ?, ?);''', (chrom, start, stop, annotation))

			ROWID = cur.lastrowid
		except sqlite3.IntegrityError: # if another worker process has inserted the same junction in between this code block's execution, then just return the junction_id from the database
			cur.execute('''select ROWID, gencode_annotation from JUNCTION_REF where 
			chromosome is ? and 
			start is ? and 
			stop is ?;''', (chrom, start, stop))

			ROWID, annotation = cur.fetchone()
		
	return ROWID, annotation

def makeSpliceDict(gene_file):

	"""
	Makes a dictionary containing junction positions and their read counts from a text file generated
	by SpliceJunctionDiscovery.py

	The format of the text file is as follows:

	chromosome	StartPos	StopPos	ReadCount

	Args:
		gene_file, the path to a text file produced by SpliceJunctionDiscovery.py

	Returns:
	    spliceDict, a dictionary of junction positions and their read counts

	Raises:
	    None
	"""

	spliceDict = {}

	with open(gene_file, "r") as gf:
		for line in gf:

			chrom, start, stop, count = line.strip().split()
			uniqueSplice = (chrom, start, stop) 

			spliceDict[uniqueSplice] = int(count)

	return spliceDict

def normalizeReadCount(spliceDict, junction, annotation, annotated_counts):

	"""
	Normalizes the read count of a splice site depending on junction's annotation

	One site is annotated, then normalize that site
	Neither site is annotated, don't perform normalization
	
	If both sites are annotated or there is a case of exon skipping, perform normalization
	on the site which has the largest read count.

	Args:
		gene_file, the path to a text file produced by SpliceJunctionDiscovery.py

	Returns:
	    norm_count, the ratio between an annotated splice site's read count to that of
	    the annotated splice site belonging to a junction with the largest read count 

	Raises:
	    None
	"""

	# gencode_annotation = {0, 1, 2, 3, 4}
	# none, only start, only stop, both, exon skipping

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
		if annotated_counts[startString] > annotated_counts[stopString]:
			key = startString
		else:
			key = stopString

	norm_count = round((float(spliceDict[junction]) / float(annotated_counts[key])), 3)

	return str(norm_count)

def makeStartString(chrom, start):

	"""
	Makes a string that distinguishes a junction position as the 5' end of a junction
	for use in a Python dictionary (hash map)

	Args:
		chrom, the chromosome a junction lies on
		start, the 5' splice site of a junction

	Returns:
	    string

	Raises:
	    None
	"""

	return ''.join([chrom,':','START',':',start])

def makeStopString(chrom, stop):

	"""
	Makes a string that distinguishes a junction position as the 3' end of a junction
	for use in a Python dictionary (hash map)

	Args:
		chrom, the chromosome a junction lies on
		stop, the 3' splice site of a junction

	Returns:
	    string

	Raises:
	    None
	"""

	return ''.join([chrom,':','STOP',':',stop])

def get_annotated_counts(spliceDict):

	"""
	Creates a dictionary containing maximum read counts for splice sites in a gene
	This function distinguishes between different ends of a junction to be consistent
	when normalizing.

	Args:
		spliceDict, a dictionary containing junctions and their read counts

	Returns:
	    count_dict, a dictionary containing maximum read counts for start and stop positions
	    of junctions in spliceDict

	Raises:
	    None
	"""

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

def summarizeGeneFile(poolArguement):

	"""
	The function each worker process must go through.

	Each process is assigned a gene, and finds the corresponding gene text file in each
	sample folder. The worker process performs transcript_model annotation, normalization, 
	and finally adds the junction information to the database.

	Args:
		poolArgument, is composed of:

			bamList, a list of bams to be used to access each sample's folder
			gene, the gene text file in which each worker process should read from each sample
			flank, the flanking region for each transcript_model junction

	Returns:
	    None

	Raises:
	    None
	"""

	bamList, gene, flank = poolArguement
	conn, cur = connectToDB()

	print ('processing ' + gene)

	for bam in bamList:

		bam_id, bam_type = get_bam_id_and_type(cur, bam)

		sample = bam[:-4]
		gene_file = ''.join([os.getcwd(), "/", sample, "/", gene, ".txt"])

		if not os.path.isfile(gene_file):
			continue

		spliceDict = makeSpliceDict(gene_file)
		annotated_counts = get_annotated_counts(spliceDict)

		for junction in spliceDict:

			chrom, start, stop = junction
			reads = spliceDict[junction]

			junction_id, annotation = getJunctionID(cur, chrom, int(start), int(stop), flank)

			try:
				norm_read_count = normalizeReadCount(spliceDict, junction, annotation, annotated_counts)
			except ZeroDivisionError:
				print("Zero division error when normalizing %s:%s-%s in genefile %s.txt in sample %s with annotation %d"%(chrom, start, stop, gene, sample, annotation))
				norm_read_count = 'null'

			annotateJunctionWithGene(gene, junction_id, cur)

			# locks are needed to maintain total read count integrity
			lock.acquire()
			updateJunctionInformation(junction_id, bam_id, bam_type, gene, sample, reads, norm_read_count, cur)
			lock.release()

		del spliceDict, annotated_counts
	
	commitAndClose(conn)

	print ('finished ' + gene)

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

def makeLockGlobal(poolLock):

	"""
	Makes a lock global for worker processes to use

	Args:
		poolLock, a multiprocessing lock to be used by worker processes

	Returns:
	    None

	Raises:
	    None
	"""

	global lock
	lock = poolLock

def addSamplesToDatabase(bam_files):

	"""
	Adds all bam files found in the text file bamList to the database.

	Returns a list of successfully added bam files. Bam files which were not successfully added
	are assumed to already exist in the database and thus their junction information should not be
	processed again. Hence they are not added to bamList.

	Args:
		bam_files, path to a file containing the names of bams to be processed in the database each on
		a seperate line

	Returns:
	    bamList, a list of successfully added bam files. This list is to be used to access folders
	    generated by SpliceJunctionDiscovery.py which are named after the bam files they contain
	    information from.

	Raises:
	    None
	"""

	conn, cur = connectToDB()
	bamList = []

	with open(bam_files, "r") as bf:
		for line in bf:

			bam = line.strip()

			try: # insert sample names into SAMPLE_REF
				if 'GTEX' in bam:
					cur.execute('''insert into SAMPLE_REF (sample_name, type) values (?, 0);''', (bam, ))
				else: # sample is a patient
					cur.execute('''insert into SAMPLE_REF (sample_name, type) values (?, 1);''', (bam, ))
			except sqlite3.IntegrityError as e:
				continue # if sample already in DB, don't process it

			bamList.append(bam) # if the script has passed over the continue statement, then the sample was succesfully added. Append to bamList.

	commitAndClose(conn)

	return bamList

def gene_file_names(transcript_file):

	"""
	Makes a list of gene names. This set is to be used by worker processes to access their corresponding
	text file.

	Args:
		transcript_file, the same transcript_file used for SpliceJunctionDiscovery.py and thus contains
		the names of each gene file produced by the script in the first column of the file

	Returns:
	    gene_set, a list of gene names to be used when accessing text files in worker processes

	Raises:
	    None
	"""

	gene_set = set() # why a set? sets can only contain unique elements

	with open(transcript_file, "r") as gf:
		for line in gf:
			gene = line.strip().split()[0]

			gene_set.add(gene)

	return gene_set

def parallel_process_gene_files(num_processes, bam_files, transcript_file, flank):

	"""
	Initializes all parameters needed for worker processes and then runs them.

	Parameters include: 
		bamList, a list of sample folder names
		gene_set, a list of gene text files
		poolLock, a global lock for worker processes

	Args:
		num_processes, the number of worker processes to run at a given time. A larger number means more ram use.
		bam_files, path to a file containing the names of all bam files to be processes
		gene_list, path to a file containing a list of gene names in its first column
		flank, the allowed +/- range for gencode annotation  

	Returns:
	    None

	Raises:
	    None
	"""

	flank = int(flank)
	poolArguements = []
	gene_set = gene_file_names(transcript_file)
	bamList = addSamplesToDatabase(bam_files)
	poolLock = multiprocessing.Lock()

	for gene in gene_set:
		poolArguements.append((bamList, gene, flank))

	print ("Creating a pool with " + str(num_processes) + " processes")
	pool = multiprocessing.Pool(initializer=makeLockGlobal, initargs=(poolLock, ), processes=int(num_processes))
	print ('pool: ' + str(pool))

	pool.map(summarizeGeneFile, poolArguements)
	pool.close()
	pool.join()

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
	parser.add_argument('-transcript_model',help="Transcript model of canonical splicing, e.g. gencode v19. Default is set to MendelianRNA-seq-DB/gencode.comprehensive.splice.junctions.txt",action='store',
			    default = "/hpf/largeprojects/ccmbio/naumenko/tools/MendelianRNA-seq-DB/gencode.comprehensive.splice.junctions.txt")
	parser.add_argument('-transcript_file',help="The same transcript_file used in SpliceJunctionDiscovery.py",action='store',
			    default = "/hpf/largeprojects/ccmbio/naumenko/tools/MendelianRNA-seq-DB/all-protein-coding-genes-no-patches.list")
	parser.add_argument('-processes',help='Number of worker processes to parse gene files, default=10.',default=10)
	parser.add_argument('-bamlist',help='A text file containing the names of bam files you want to discover splice junctions in each on a seperate line, default=bamlist.list',default='bamlist.list')
	parser.add_argument('-flank',help='Add a +/- flanking region for gencode annotation. Specify 0 if you don\'t want to use this feature, default=1',default=1)
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
		print ('Storing junctions from bam files found in the file ' + args.bamlist)
		parallel_process_gene_files(args.processes, args.bamlist, args.transcript_file, args.flank)
	elif args.delete:
		sample = args.sample

		if not sample:
			print('Please enter a sample name with its .bam extension using the parameter \'-sample SAMPLE_NAME\'')
			exit(1)

		deleteSample(sample)

	print ('AddJunctionsToDatabase.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
