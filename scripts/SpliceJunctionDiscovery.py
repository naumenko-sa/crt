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

def printJunctions(path, spliceDict,gene):
	"""
		spliceDict, a dictionary containing junctions and their read counts
			E.x. spliceDict[1:200-300] = 5
	"""
	for key in spliceDict:
		chrom, junctionStart, junctionEnd = key
		timesSeenInSample = str(spliceDict[key])

		with open(path, "a") as out:
			out.write("\t".join([str(chrom),str(junctionStart),str(junctionEnd),timesSeenInSample,gene])+"\n")

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

def getJunctionsForAGene(bam,gene,chrom,start,end):
	"""
	Produces a file containing junction
	positions and their read counts for each gene:
		E.x. NPHS1.txt
			1	200	300	5
			1	344	355	2
	"""
	pos = ''.join([chrom, ':', start, '-', end])

	spliceDict = {}
	geneFilePath = (bam + ".junctions.txt")

	try:
		exitcode, stdout, stderr = run(' '.join(['samtools view', bam, pos]))
	except Exception as e:
		print ('Exception message: ' + str(e))
		print ("Exception occured while running \"samtools view\" on " + bam + " for position " + pos + " Skipping.")

	for line in stdout.splitlines():
		elems = line.decode().split()

		alignmentStart = int(elems[3])
		cigar = str(elems[5])
		alignmentScore = int(elems[1])
 
		if 'N' not in cigar:  	#only get introns
			continue

		if (alignmentScore >= 256):  	#only primary alignments
			continue

		if not ((alignmentStart > int(start)) and (alignmentStart < int(end))):  	#check if alignment start is after known junction start but before known junction end 
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
		printJunctions(geneFilePath,spliceDict,gene)
		del spliceDict

def getJunctions(genes, bam):
	with open(genes) as tf:
		for line in tf:
			elems = line.strip().split()
			chrom, start, end, gene = elems #edit the transcript file so that you only deal with junction coordinates
			getJunctionsForAGene(bam, gene, chrom, start, end)

if __name__=="__main__":
	print ('SpliceJunctionDiscover.py started on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
	parser = argparse.ArgumentParser(description = 'Discover splice junctions in a bam file')
	parser.add_argument('-genes',help="A bed file of gene coordinates to discover junctions in",
						action='store',
			    		default = "/home/naumenko/crt/genes.bed")
	parser.add_argument('-bam',help='A bam file')
	args=parser.parse_args()

	getJunctions(args.genes, args.bam)
	
	# transcriptFile = str(args.transcriptFile).rsplit('/')[-1] #remove paths

	print ('SpliceJunctionDiscover.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
