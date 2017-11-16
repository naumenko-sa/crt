import sys
import os
import errno
import argparse
import multiprocessing
import subprocess
from subprocess import Popen, PIPE
from cigar import Cigar
from datetime import datetime

def run(cmd, dieOnError=True):

	"""
	Runs a single subprocess in the background. Used to execute simple bash commands like cat *.txt

	run() was taken from Andy and modified:
	http://jura.wi.mit.edu/bio/education/hot_topics/python_pipelines_2014/python_pipelines_2014.pdf
	https://stackoverflow.com/questions/13398261/python-subprocess-call-and-subprocess-popen-stdout

	Args:
		cmd, the terminal command to run
		dieOnError, used to detect errors

	Returns:
	    (exitcode, stdout, stderr), tuple containing the results of the subprocess' execution

	Raises:
	    None
	"""

	ps = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
	exitcode = ps.returncode
	stdout,stderr = ps.communicate()
	return exitcode, stdout, stderr

def printSplices(path, spliceDict):

	"""
	Prints junctions and their read counts to a specified file

	This is a striped down version of Beryl Cummings' printSplices()

	Args:
		path, the path to the output file
		spliceDict, a dictionary containing junctions and their read counts
			E.x. spliceDict[1:200-300] = 5

	Returns:
	    None

	Raises:
	    None
	"""

	for key in spliceDict:
		chrom, junctionStart, junctionEnd = key
		timesSeenInSample = str(spliceDict[key])

		with open(path, "a") as out:
			out.write("\t".join([str(chrom),str(junctionStart),str(junctionEnd),timesSeenInSample])+"\n")

def parseCIGARForIntrons(cigar):

	"""
	Parses a CIGAR string and returns values which can used to determine an intron's
	3' and 5' splice sites

	Args:
		cigar, a CIGAR string with an intron in it
			E.x. cigar='3M1D40M20N'

	Returns:
	    offset, a figure which accomodates for insertion and deletion events to adjust
	    an alignment's positions back to the reference genome

	    matchedExon, a figure to be added to the start position of an alignment which
	    forms the 5' end of a splice site. This function only considers 'M''s before 
	    an intron (N) for this figure.

		intronLength, the length of an intron as reported by the CIGAR string. This figure
		is added with matchedExon and the start of an alignment to produce the position of
		the 3' end of a splice site.

	Raises:
	    None
	"""

	if 'N' in cigar:
		cigar = cigar.split('N')[0] + 'N' #remove all information after intron
	else:
		raise Exception('No intron detected')

	offset = 0
	matchedExon = 0
	intronLength = 0
	
	for c in list(Cigar(cigar).items()): # returns list of tuples : [(20, 'N')]
		if c[1] == 'N':
			intronLength += int(c[0])
		elif c[1] == 'D':
			offset += int(c[0])
		elif c[1] == 'I':
			offset -= int(c[0])
		elif c[1] == 'M':
			matchedExon += int(c[0])
		## soft clipping is ignored
		## hard clipping is ignored too

	return offset, matchedExon, intronLength

def intronDiscovery(poolArguement):

	"""
	The function a worker process goes through. Produces a file containing junction
	positions and their read counts for each gene:

		E.x. NPHS1.txt
			1	200	300	5
			1	344	355	2

	Each worker is assigned a single gene region to work on. The worker loops through each BAM file,
	counts the number of alignments and produces a single gene text file in the corresponding BAM folder. 

	Args:
		poolArguement, the single argument for each worker process, which can be broken down
		in to these components:

			bamFiles, a list of each bam file, this is used to access the corresponding sample folder
			gene, the gene the worker must process
			chrom, the chromosome the gene lies on
			start and stop, the 3' and 5' locations on a chromosome in which samtools should begin looking for alignments in
			cwd, path to the current working directory. This is used to create the path of a sample folder and a gene text file

	Returns:
	    None

	Raises:
	    None
	"""

	bamFiles, gene, chrom, start, stop, cwd = poolArguement

	print ('processing ' + gene)

	pos = ''.join([chrom, ':', start, '-', stop])

	for bam in bamFiles:

		spliceDict = {}
		geneFilePath = (cwd + "/" + bam[:-4] + "/" + gene + ".txt")

		try:
			exitcode, stdout, stderr = run(' '.join(['samtools view', bam, pos]))
		except Exception as e:
			print ('Exception message: ' + str(e))
			print ("Exception occured while running \"samtools view\" on " + bam + " for position " + pos + " Skipping.")
			continue

		if not stdout:
			#print ('No introns found for ' + gene + ' at ' + pos + ' in ' + bam)
			continue

		for line in stdout.splitlines():

			elems = line.decode().split()

			alignmentStart = int(elems[3])
			cigar = str(elems[5])
			alignmentScore = int(elems[1])
 
			if 'N' not in cigar:  	#only get introns
				continue

			if (alignmentScore >= 256):  	#only primary alignments
				continue

			if not ((alignmentStart > int(start)) and (alignmentStart < int(stop))):  	#check if alignment start is after known junction start but before known junction end 
				continue

			try:
				offset, matchedExon, intronLength = parseCIGARForIntrons(cigar)
			except Exception as e:
				print ('Error message: ' + str(e))
				print ('Error trying to parse CIGAR string: ' + cigar +  ' with the bam file ' + bam +  ' and the position: ' + pos + ' Skipping.')
				continue

			junctionStart = alignmentStart + matchedExon + offset
			junctionEnd = junctionStart + intronLength

			# Beryl Cummings' Code, taken from makeUniqSpliceDict()
			# uniqueSplice = ':'.join([chrom, str(junctionStart), str(junctionEnd)])
			uniqueSplice = (chrom, str(junctionStart), str(junctionEnd))
			
			if uniqueSplice not in spliceDict:
				spliceDict[uniqueSplice] = 1
			else:
				spliceDict[uniqueSplice] += 1

		del stdout # saves ram in between samtool calls

		if spliceDict:
			printSplices(geneFilePath, spliceDict)
			del spliceDict

	print ('finished ' + gene)

def makeBamListAndDirectories(bamList):

	"""
	Makes a list of bam files for each worker process to use and
	the directories in which each gene file will reside in

	Args:
		bamList, text file containing the names of all .bam files
		each on a seperate line

	Returns:
	    bamFiles, a list of bamFiles to be processed

	Raises:
	    None
	"""

	bamFiles = []

	with open(bamList) as bl:
		for i in bl:

			i = i.strip()

			bamLocation = os.getcwd() + '/' + i
			if not os.path.isfile(bamLocation):
				print ('bam file: ' + i + ' does not exist in CWD! Skipping.')
				continue

			outputDirectory = bamLocation[:-4]

			os.system("mkdir " + outputDirectory)
			bamFiles.append(i)

	return bamFiles

def processGenesInParallel(transcriptFile, bamList, numProcesses):

	"""
	Sets up the parameters for each worker process and then runs them.

	Args:
		transcriptFile, path to a file which contains a list of genes and locations of investigation
		bamList, a list of bam files you want to discover splice sites in
		numProcesses, the number of worker processes to run at a given time

	Returns:
	    None

	Raises:
	    None
	"""

	cwd = os.getcwd()
	bamFiles = makeBamListAndDirectories(bamList)
	poolArguements = []

	with open(transcriptFile) as tf:
		for line in tf:

			elems = line.strip().split()
			try:
				gene, gene2, plus, chrom, start, stop, gene_type = elems #edit the transcript file so that you only deal with junction coordinates
			except Exception as e:
				print ('Error while parsing transcript file named: ' + str(transcriptFile) + "\n" + 'Error message: ' + str(e) + "\nExiting.")
				exit (3)

			poolArguements.append((bamFiles, gene, chrom, start, stop, cwd))

	print ("Creating a pool with " + str(numProcesses) + " processes")
	pool = multiprocessing.Pool(int(numProcesses))
	print ('pool: ' + str(pool))

	pool.map(intronDiscovery, poolArguements) # run the worker processes
	pool.close()
	pool.join()
	
if __name__=="__main__":

	print ('SpliceJunctionDiscover.py started on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))

	parser = argparse.ArgumentParser(description = 'Discover splice junctions from a list of bam files')
	parser.add_argument('-transcript_file',help="A list of positions that you want to discover junctions in",action='store',
			    default = "/hpf/largeprojects/ccmbio/naumenko/tools/MendelianRNA-seq-DB/all-protein-coding-genes-no-patches.list")
	parser.add_argument('-bam_list',help='A text file containing the names of bam files you want to discover splice junctions in each on a seperate line',default='bamlist.txt')
	parser.add_argument('-processes',help='number of processes to run multiple instances of: "samtools view", default=10',default=10)
	args=parser.parse_args()

	print ('Working in directory' + str(os.getcwd()))
	print ('Transcript file is ' + str(args.transcript_file))
	print ('Identifying splice junction is ' + str(args.bam_list))

	processGenesInParallel(args.transcript_file, args.bam_list, args.processes)
	
	# transcriptFile = str(args.transcriptFile).rsplit('/')[-1] #remove paths

	print ('SpliceJunctionDiscover.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
