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

	# Extract the header of the input SAM file

	os.system('sed -n \'/^@/p;/^[^@]/q\' ' + inputFile + ' > ' + outputFile)

	# Remove the header in the temp file

	os.system('sed -n \'/^@/d;/^[^@]/p\' ' + inputFile + ' | sort >> ' + outputFile)
	
	return True

# Count file lines

def opcount(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

# Get the qname tag

def QNameKey(qname, keyRe):
	key = qname
	if keyRe :
		if(keyRe.match(qname)):
			tags = keyRe.findall(key)
			if(len(tags) > 0):
				key = ''
				if(type(tags[0]) is str):
					key = tags[0]
				elif(type(tags[0]) is tuple):
					for s in tags[0] :
						key += s
	return key

# Calculate the key an alignment
# The algrithm should make sure that one key identify one unique read

def AlignmentKey(alignment, keyRe):
	# Get the tag of QNAME without read number
	key = QNameKey(alignment.qname, keyRe)

	if(alignment.flag & 0x40):
		key += (':' + str(alignment.pos) + ':' + str(alignment.pnext))
	elif(alignment.flag & 0x80):
		key += (':' + str(alignment.pnext) + ':' + str(alignment.pos))

	return key

def AlignmentGroupKey(alignment, keyRe):
	return QNameKey(alignment.qname, keyRe)

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

		nums = samparser.tagNumRe.findall(cigar)
		tags = samparser.cigarTagRe.findall(cigar)

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
	if(not samparser.mdRe.match(mdstr)):

		# unavaliable MD value

		return 255
	else:

		# get match numbers

		tags = samparser.mdTagRe.findall(mdstr)
		

		ntotal = 0
		nmatch = 0
		for tag in tags:
			if(tag.isdigit()):
				tagLen = int(tag)
				nmatch += tagLen
				ntotal += tagLen
			elif(tag.isalpha()):
				ntotal += len(tag)

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
		else:
			return False
		return True

	def str(self):
		alignmentStr = ''
		if(self.read1):
			alignmentStr += self.read1.str()

		if(self.read1 and self.read2):
			alignmentStr += '\n'

		if(self.read2):
			alignmentStr += self.read2.str()

		return alignmentStr

	def segmentLen(self):
		readLen = 0
		if(self.read1):
			readLen = len(self.read1.seq)
		if(self.read2):
			readLen = len(self.read2.seq)
		if(readLen == 0):
			return 0

		# If read1 or read2 is not paired properly, return 
		# the read length of either of them

		if(not(self.read1 and self.read2)):
			return readLen

		# If the reads are paired, we define the RC read as - strand
		# and non RC as the + strand. We will find the start positions
		# of both + read and - read.

		if((self.read1.flag & 0x10) == (self.read2.flag & 0x10)):
			# read1 and read2 on the same strand
			return 0

		posL = 0
		posR = 0
		if(self.read1.flag & 0x10):
			posR = self.read1.pos + readLen
			posL = self.read2.pos
		else:
			posR = self.read2.pos + readLen
			posL = self.read1.pos

		return (posR - posL)

#
# [*] Unique Strategy:
# Following strategies are used to find the unique & the best alignment
#
# 1. Keep the alignment pair that has the highest score. If more than one pairs 
#    are found to have the same "Highest Score", these pairs will be removed. 
#
# 2. Read1 and Read2 should be mapped on different strands.
#
# 3. The segment length decided by the read pairs should be longer than 
#    0.7 * read length
#
# [*] Log file Specification
# 
# Symbol	Description
# -------------------------------------------------
# !		Error lines
# <		Low score alignments
# =		Pairs with more than one best score
# ~		Read pair mapped on the same strand
# ?		Segment length too short
# 

def UniquePairs(pairs, outfile, logfile):
	bestPair = None
	bestScore = -1
	scoreCount = 0

	writeCount = 0

	for key in pairs :
		pairRead = pairs[key]
		if(pairRead.score > bestScore):
			bestScore = pairRead.score
			scoreCount = 1
			bestPair = pairRead
		elif(pairRead.score == bestScore):
			scoreCount += 1

	# no best pair found

	if(not bestPair or scoreCount <= 0):
		return 0

	# Rule No. 1: keep the best pair

	for key in pairs:
		pairRead = pairs[key]
		if((pairRead.score == bestScore) and (scoreCount >= 2)):
			if(pairRead.read1):
				logfile.write('= ' + pairRead.read1.str() + '\n')
			if(pairRead.read2):
				logfile.write('= ' + pairRead.read2.str() + '\n')
		elif(pairRead.score < bestScore):
			if(pairRead.read1):
				logfile.write('< ' + pairRead.read1.str() + '\n')
			if(pairRead.read2):
				logfile.write('< ' + pairRead.read2.str() + '\n')

	if(scoreCount >= 2):
		return 0

	if(bestPair.read1 and bestPair.read2):

		# Rule No. 2: r1 and r2 should be on different strands
		
		if((bestPair.read1.flag & 0x10) == (bestPair.read2.flag & 0x10)):
			logfile.write('~ ' + pairRead.read1.str() + '\n')
			logfile.write('~ ' + pairRead.read2.str() + '\n')
			return 0

		# Rule No. 3: segment length >= 0.7 * read length

		if(bestPair.segmentLen() <= int(round(0.7 * len(bestPair.read1.seq)))):
			logfile.write('? ' + pairRead.read1.str() + '\n')
			logfile.write('? ' + pairRead.read2.str() + '\n')
			return 0

	# Output this pair

	if(bestPair.read1):
		outfile.write(bestPair.read1.str() + '\n')
		writeCount += 1
	if(bestPair.read2):
		outfile.write(bestPair.read2.str() + '\n')
		writeCount += 1

	return writeCount


def main():
	# parse the command line options
	usage = 'usage: %prog [options] input.sam -o output.sam'
	parser = OptionParser(usage=usage, version='%prog version 0.1.6')
	parser.add_option('-o', '--output-file', dest='outputfile',
						help='write the result to output file')
	parser.add_option('-s', '--sort', 
						action="store_true", dest="sort", default=False,
						help='sort the input SAM file before further processing')
	parser.add_option('-k', '--key-reg', dest="keyreg",
						help='qname regular expression to extract the alignment key')

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

	# Prepare the qname key

	alignmentKeyReg = None
	if(options.keyreg):
		alignmentKeyReg = re.compile(options.keyreg)

	# prepare the output file

	outputFileName = 'unique.' + os.path.basename(inputFileName)
	if(options.outputfile):
		outputFileName = options.outputfile

	try:
		outfile = open(outputFileName, 'w')
	except IOError:
		print('error: write to output file failed!')
		sys.exit(-1)

	logFileName = os.path.basename(inputFileName) + '.log'
	try:
		logfile = open(logFileName, 'w')
	except IOError:
		print('error: create log file failed!')
		sys.exit(-1)

	# processing file

	print('* Analyzing...')
	totalLineCount = opcount(samFileName)
	print('  %ld lines found' % totalLineCount)
	lineCount = 0
	writtenLineCount = 0
	print('* Processing...')
	pairs = {}
	currentGroupKey = ''
	with open(samFileName) as samFile:
		for line in samFile:
			# build alignment dictionary 

			if(samparser.blankLineRe.match(line)):
				continue
			elif(line[0] == '@'):
				outfile.write(line)
			else:
				alignment = samparser.SamAlignment()
				if(alignment.parse(line.strip())):
					
					# Write result

					groupKey  = AlignmentGroupKey(alignment, alignmentKeyReg)
					if(groupKey != currentGroupKey):
						currentGroupKey = groupKey
						writtenLineCount += UniquePairs(pairs, outfile, logfile)
						pairs.clear()

					# Pair up

					key = AlignmentKey(alignment, alignmentKeyReg)
					if(key in pairs):
						readPair = pairs[key]
					else:
						readPair = ReadPair()
						pairs[key] = readPair

					# alignment without map information (read1 or read2)

					if(not readPair.add(alignment)):
						logfile.write('- ' + line)

				else:
					logfile.write('! ' + line)
					print('error: Encountered unknown alignment line: "', line, '"')
					sys.exit(-1)

			# progress information 

			lineCount = lineCount + 1
			if(totalLineCount == 0):
				percentage = 0
			else:
				percentage = lineCount * 1.0 / totalLineCount
			sys.stdout.write('\r  line %ld (%.2f%%)' % (lineCount, percentage * 100))
			sys.stdout.flush()

		# write the alignments in the list

		writtenLineCount += UniquePairs(pairs, outfile, logfile)

	print('\n  %ld alignments written' % (writtenLineCount))
	samFile.close()
	outfile.close()
	logfile.close()
	
	# Clear resources

	if(isSamFileTemp):
		os.unlink(samFileName)

	print('* Complete')


if __name__ == '__main__':
	main()