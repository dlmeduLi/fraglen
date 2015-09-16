#!/usr/bin/env python

from __future__ import print_function

import re
import os
import os.path
import optparse
import sys
import pysam
import subprocess

# global regular expressions

readKeyRe = re.compile('(.*)[\s\/][1-2](.*)')
readNumRe = re.compile('.*[\s\/]([1-2]).*')

# sort bam file with pysam.sort function
# the bam file should be sorted by qname

def SortSam(inBam, outBam):
	pysam.sort("-n", inBam, outBam)

def TrimReadSeq(seq, cigar):
	trimedSeq = ''
	pos = 0
	for m in cigar :
		if(m[0] == 0) : # M
			trimedSeq += seq[pos : pos + m[1]]
			pos += m[1]
		elif(m[0] == 1): # I
			pos += m[1]
		elif(m[0] == 2): # D
			trimedSeq += ('N' * m[1])
		elif(m[0] == 3): # N
			trimedSeq += ('N' * m[1])
		elif(m[0] == 4): # S
			pos += m[1]
	return trimedSeq

def FormatPosList(posList):
	if(len(posList) == 0):
		return 'NA'
	return ','.join(map(str, posList))


##### write fragement length #####
def WriteFragLen(fragType,chrname, pos, fragLen, outFile):
#	outFile.write('%s\t%s\t%s\t%15ld\t%6d\n' % (frg, fragType, chrname, pos, fragLen))
	outFile.write('%s\t%s\t%15ld\t%6d\n' % (fragType, chrname, pos, fragLen))

##### get fragement length from SE resd #####
def FragOfSingleReads(dictSingle, bamFile, outFile):

	for read in dictSingle.itervalues():
		readSeq = TrimReadSeq(read.seq, read.cigar)
		fragLen = len(readSeq)
		readPos = read.pos
		tag = read.qname
		chrname = bamFile.getrname(read.rname)
		WriteFragLen('S', chrname, readPos, fragLen, outFile)

##### get fragement length from PE resd #####
def FragOfPairedReads(dictPaired, bamFile, outFile):
	for pair in dictPaired.itervalues():
		if(not (pair[0] and pair[1])):
			return		

		chrname = bamFile.getrname(pair[0].rname)
		readSeq1 = TrimReadSeq(pair[0].seq, pair[0].cigar)
		readSeq2 = TrimReadSeq(pair[1].seq, pair[1].cigar)
		readSeqLen1 = len(readSeq1)
		readSeqLen2 = len(readSeq2)
		fragPos = pair[0].pos
		fragLen = pair[1].pos + readSeqLen2 - fragPos
				
		WriteFragLen('P', chrname, fragPos, fragLen, outFile)

def main():

	# parse the command line options

	usage = 'usage: %prog [options] input.bam -o output.csv'
	parser = optparse.OptionParser(usage=usage, version='%prog 0.1.0')
	parser.add_option('-o', '--output-file', dest='outputfile',
						help='write the result to output file')
	parser.add_option('-s', '--sort', 
						action="store_true", dest="sort", default=False,
						help='sort the input BAM file before data processing')

	(options, args) = parser.parse_args()
	if(len(args) != 1):
		parser.print_help()
		sys.exit(0)
	
	inputBamFileName = args[0]
	bamFileName = inputBamFileName
	baseFileName = os.path.splitext(os.path.basename(inputBamFileName))[0]
	outputFileName =  baseFileName + '.fraglen.csv'
	logFileName = baseFileName + '.fraglen.log'

	# sort the input bam file if -s option is set

	rmTemp = False
	if(options.sort):
		print('[*] Sorting by QNAME...')
		bamFileName = 'sorted.' + baseFileName
		SortSam(inputBamFileName, bamFileName)
		bamFileName += '.bam'
		rmTemp = True

	# load input files

	print('[*] Initializing...')

	if(not os.path.exists(bamFileName)):
		print('error: Failed to open file "', bamFileName, '"')
		sys.exit(-1)
	bamFile = pysam.AlignmentFile(bamFileName, "rb")

	# prepare output files
	if(options.outputfile):
		outputFileName = options.outputfile
	try:
		outFile = open(outputFileName, 'w')
	except IOError:
		print('error: write to output file failed!')
		sys.exit(-1)
	outFile.write('type\tchr\tpos\tlen\n')

	# analyse algnments

	print('[*] Analyzing...')
	
	dictSingle = {}
	dictPaired = {}
	currentGroupKey = ''
	readCount = 0
	for read in bamFile.fetch(until_eof = True):
		try:
			groupKey = ''.join(readKeyRe.findall(read.qname)[0])
		except:
			groupKey = read.qname

		if(groupKey != currentGroupKey):
			currentGroupKey = groupKey

			# handle and write results
			FragOfSingleReads(dictSingle, bamFile, outFile)
			FragOfPairedReads(dictPaired, bamFile, outFile)
			dictPaired.clear()
			dictSingle.clear()

		# check if it is properly paired

		if(read.is_proper_pair):

			# paired reads

			chrpos = bamFile.getrname(read.rname).strip() + str(read.pos)
			chrposNext = bamFile.getrname(read.rnext).strip() + str(read.pnext)
			if(chrpos < chrposNext):
				readKey = groupKey + ':' + chrpos + ':' + chrposNext
			else:
				readKey = groupKey + ':' + chrposNext + ':' + chrpos

			if(read.is_reverse):
				readType = 1
			else:
				readType = 0
			
			if(readType == 0 or readType == 1):
				if(readKey in dictPaired):
					dictPaired[readKey][readType] = read
				else:
					dictPaired[readKey] = [None, None]
					dictPaired[readKey][readType] = read
		else:

			# single read

			try:
				readKey = groupKey + ':' + readNumRe.findall(read.qname)
			except:
				readKey = groupKey
			dictSingle[readKey] = read

		# progress

		readCount += 1
		sys.stdout.write('\r    read: #%ld' % (readCount))
		sys.stdout.flush()

	FragOfSingleReads(dictSingle, bamFile, outFile)
	FragOfPairedReads(dictPaired, bamFile, outFile)

	sys.stdout.write('\n')
	
	# release resources

	bamFile.close()
	outFile.close()
	
	if(rmTemp):
		os.unlink(bamFileName)

	# ggplot2 
	subprocess.call (["Rscript","plot_fraglen.R",outputFileName,baseFileName])

	print('[*] Complete')

if __name__ == '__main__':
	main()
