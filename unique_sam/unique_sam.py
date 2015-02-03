#!/usr/bin/python

from __future__ import print_function

import re
import sys
import os
import os.path
import subprocess
from optparse import OptionParser
import shelve
from libsam import samparser

# Sort SAM file with system sort program
# sed -i '/^@/d;/!^@/q' file

def SortSamFile(inputFile, outputFile):
	if(not os.path.exists(inputFile)):
		return False

	# Make a copy of the original input file

	tempFile = 'tmp.' + os.path.basename(inputFile) + '.~'
	os.system('cp ' + inputFile + ' ' + tempFile)

	# Extract the header of the input SAM file

	os.system('sed -n \'/^@/p;/!^@/q\' ' + inputFile + ' > ' + outputFile)

	# Remove the header in the temp file

	os.system('sed -i \'/^@/d;/!^@/q\' ' + inputFile)

	# Sort temp file and redirect the output 

	os.system('sort ' + tempFile + ' >> ' + outputFile)

	# Remove temp files

	os.unlink(tempFile)

	return True

# Count file lines

def opcount(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

# Calculate the key an alignment
# The algrithm should make sure that one key identify one unique read

def AlignmentKey(alignment):
	# Get the tag of QNAME without read number
	key = alignment.qname
	if(re.match('.*\/[1-2]$', alignment.qname)):
		tags = re.split('/', alignment.qname)
		key = tags[0]

	if(alignment.flag & 0x40):
		key += (':' + str(alignment.pos) + ':' + str(alignment.pnext))
	elif(alignment.flag & 0x80):
		key += (':' + str(alignment.pnext) + ':' + str(alignment.pos))

	return key

def AlignmentGroupKey(alignment):
	key = alignment.qname
	if(re.match('.*\/[1-2]$', alignment.qname)):
		tags = re.split('/', alignment.qname)
		key = tags[0]

	return key

#
# EvaluateAlignmentCigar:
#		Calculate the goodness of an alignment by CIGAR and return a score for it,
#		this value will be used to sort alignments of the read.
#		
#		The retur value of this function should be in range [0, 2^-16 -1]

def EvaluateAlignmentCigar(cigar):	
	if(cigar == '*'):

		# cigar unavilable

		return 255
	else:

		# calculate a score from cigar

		nums = re.findall('[0-9]+', cigar)
		tags = re.findall('[MIDNSHPX=]', cigar)

		if(len(nums) != len(tags)):
			return 255

		i = 0
		ntotal = 0
		nmatch = 0
		while (i < len(nums)):
			if(tags[i] == 'M'):
				nmatch = nmatch + int(nums[i])
			ntotal = ntotal + int(nums[i])
			i = i + 1

		if(ntotal != 0):
			return int(round(100.0 * nmatch / ntotal))
		else:
			return 255

#
# EvaluateAlignmentMD:
#		Calculate the goodness of an alignment by MD tag and return a score for it,
#		this value will be used to sort alignments of the read.
#		
#		The retur value of this function should be in range [0, 2^-16 -1]

def EvaluateAlignmentMD(mdstr):	
	if(not re.compile('[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*').match(mdstr)):

		# unavaliable MD value

		return 255
	else:

		# get match numbers

		nums = re.findall('[0-9]+', mdstr)
		nmis = len(re.findall('[A-Z]', mdstr))
		
		i = 0
		ntotal = 0
		nmatch = 0
		while (i < len(nums)):
			nmatch += int(nums[i])
			i = i + 1
		ntotal = nmatch + nmis

		if(ntotal != 0):
			return int(round(100.0 * nmatch / ntotal))
		else:
			return 255

#
# EvaluateAlignment: -- By MD
#		Calculate the goodness of an alignment by MD and return a score for it,
#		this value will be used to sort alignments of the read.
#		
#		The retur value of this function should be in range [0, 2^-16 -1]

def EvaluateAlignment(alignment):
	if(alignment.mapq == 255):
		# tha mapq value is unavilable

		if('MD' in alignment.tags):
			mdstr = alignment.tags['MD'].value
			return EvaluateAlignmentMD(mdstr)
		else:
			return EvaluateAlignmentCigar(alignment.cigar)
	else:
		# use existing MAPQ value

		return int(round(alignment.mapq))

##
## ReadPair Object
##

class ReadPair(object):
	def __init__(self):
		self.read1 = None
		self.read2 = None
		self.score = 255

	def updateScore(self):
		self.score = 0
		if(self.read1):
			self.score  += EvaluateAlignment(self.read1)
		
		if(self.read2):
			self.score += EvaluateAlignment(self.read2)
	
	def add(self, alignment):
		if(alignment.flag & 0x40):
			self.read1 = alignment
			self.updateScore()
		elif(alignment.flag & 0x80):
			self.read2 = alignment
			self.updateScore()

	def str(self):
		alignmentStr = ''
		if(self.read1):
			alignmentStr += self.read1.str()

		if(self.read1 and self.read2):
			alignmentStr += '\n'

		if(self.read2):
			alignmentStr += self.read2.str()

		return alignmentStr

#
# Unique Strategy:
# Following strategies are used to find the unique & the best alignment
#
# 1. Keep the alignment pair that has the highest score. If more than one pairs 
# are found to have the same "Highest Score", these pairs will be removed. 
#

def UniquePairs(pairs, outfile):
	bestPair = None
	bestScore = -1
	scoreCount = 0

	for key in pairs :
		pairRead = pairs[key]
		if(pairRead.score > bestScore):
			bestScore = pairRead.score
			scoreCount = 1
			bestPair = pairRead
		elif(pairRead.score == bestScore):
			scoreCount += 1

	# Rule No. 1: keep the best pair

	if(scoreCount <= 0 or scoreCount >= 2):
		return 0

	if(not bestPair):
		return 0

	# Rule No. 2: 

	# Output this pair

	outfile.write(bestPair.str() + '\n')

	return 2


def main():
	# parse the command line options
	usage = 'usage: %prog [options] input.sam -o output.sam'
	parser = OptionParser(usage=usage, version='%prog 1.0.0')
	parser.add_option('-o', '--output-file', dest='outputfile',
						help='write the result to output file')
	parser.add_option('-s', '--sort', 
						action="store_true", dest="sort", default=False,
						help='sort the input SAM file before further processing')

	(options, args) = parser.parse_args()
	if(len(args) != 1):
		parser.print_help()
		sys.exit(0)
	inputFileName = args[0]
	samFileName = inputFileName
	isSamFileTemp = False

	if(options.sort):
		print('* Sorting...')
		samFileName = 'sorted.' + os.path.basename(inputFileName)
		if(not SortSamFile(inputFileName, samFileName)):
			print('error: Failed to sort file "', inputFileName, '"')
			sys.exit(-1)
		isSamFileTemp = True

    # Load the sam file

	if(not os.path.exists(samFileName)):
		print('error: Failed to open file "', samFileName, '"')
		sys.exit(-1)

	# prepare the output file

	outputFileName = 'unique.' + os.path.basename(inputFileName)
	if(options.outputfile):
		outputFileName = options.outputfile

	try:
		outfile = open(outputFileName, 'w')
	except IOError:
		print('error: write to output file failed!')
		sys.exit(-1)

	# processing file

	print('* Analyzing...')
	totalLineCount = opcount(samFileName)
	print('  %ld lines found.' % totalLineCount)
	lineCount = 0
	writtenLineCount = 0
	print('* Processing...')
	pairs = {}
	currentGroupKey = ''
	with open(samFileName) as samFile:
		for line in samFile:
			# build alignment dictionary 

			if(re.match('^\s*$', line)):
				continue
			elif(line[0] == '@'):
				outfile.write(line)
			else:
				alignment = samparser.SamAlignment()
				if(alignment.parse(line.strip())):
					
					# Write result

					groupKey  = AlignmentGroupKey(alignment)
					if(groupKey != currentGroupKey):
						currentGroupKey = groupKey
						writtenLineCount += UniquePairs(pairs, outfile)
						pairs.clear()

					# Pair up

					key = AlignmentKey(alignment)
					if(key in pairs):
						readPair = pairs[key]
					else:
						readPair = ReadPair()
						pairs[key] = readPair

					readPair.add(alignment)

				else:
					print('error: Encountered unknown alignment line: "', line, '"')
					sys.exit(-1)
			
			# progress information 

			lineCount = lineCount + 1
			if(totalLineCount == 0):
				percentage = 0
			else:
				percentage = lineCount * 1.0 / totalLineCount
			sys.stdout.write('\r  line %ld (%.2f%%)' % (lineCount + 1, percentage * 100))
			sys.stdout.flush()

	print('\n  %ld lines written' % (writtenLineCount))
	outfile.close()
	
	# Clear resources

	if(isSamFileTemp):
		os.unlink(samFileName)

	print('* Complete')


if __name__ == '__main__':
	main()