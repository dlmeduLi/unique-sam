#!/usr/bin/python

from __future__ import print_function

import re
import sys
import os
import os.path
from optparse import OptionParser
import shelve
from libsam import samparser

# Count file lines

def opcount(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

# Calculate the key an alignment
# The algrithm should make sure that one key identify one unique read

def AlignmentKey(alignment):
	rCode = ''
	if(alignment.flag & 0x40):
		rCode = '1'
	elif(alignment.flag & 0x80):
		rCode = '2'
	if(rCode == ''):
		return alignment.qname
	else:
		return alignment.qname + '/' + rCode

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

def main():
	# parse the command line options
	usage = 'usage: %prog [options] samfile'
	parser = OptionParser(usage=usage, version='%prog 1.0.0')
	parser.add_option('-o', '--output-file', dest='outputfile',
						help='write the output to outputfile insted of stdio')

	(options, args) = parser.parse_args()
	if(len(args) != 1):
		parser.error('incorrect number of arguments')
	samFileName = args[0]

    # Load the sam file

	if(not os.path.exists(samFileName)):
		print('error: Failed to open file "', samFileName, '"')
		sys.exit(-1)

	sam = samparser.Sam()
	alignments = {}
	dbFileName = '.tmp.' + os.path.basename(samFileName) + '.db.~'
	db = shelve.open(dbFileName)
	print('* Analyzing...')
	totalLineCount = opcount(samFileName)
	print('  %ld lines found.' % totalLineCount)
	linecount = 0
	print('* Processing...')
	with open(samFileName) as samFile:
		for line in samFile:
			# build alignment dictionary 

			if(re.match('^\s*$', line)):
				continue
			elif(line[0] == '@'):

				# parse and  store header line
				
				sam.parseLine(line)

			else:
				alignment = samparser.SamAlignment()
				if(alignment.parse(line.strip())):
					
					# build dictionary

					key = AlignmentKey(alignment)
					score = EvaluateAlignment(alignment)
					alignment.mapq = score

					if(key in alignments):
						# compare the mapq value

						# oldAlignment = alignments[key]
						# currentScore = (oldAlignment.mapq & 0xFF)
						value = alignments[key]
						currentScore = (value & 0xFF)
						if(score > currentScore):
							
							# replace current record with the higher score record
							
							#alignments[key] = alignment
							alignments[key] = score
							db[key] = line
						elif(score == currentScore):
							
							# update the time repeat field of mapq value

							count = ((value & 0xFFFF0000) >> 16)
							count = count + 1
							alignments[key] = ((count << 16) ^ currentScore)
					else:
						alignments[key] = score
						db[key] = line
				else:
					print('error: Encountered unknown alignment line: "', line, '"')
					sys.exit(-1)
			linecount = linecount + 1
			if(totalLineCount == 0):
				percentage = 0
			else:
				percentage = linecount * 1.0 / totalLineCount
			sys.stdout.write('\r  line %ld (%.2f%%)' % (linecount + 1, percentage * 100))
			sys.stdout.flush()

	# Output the result
	print('\n* Writing results...')
	outputFileName = 'unique.' + os.path.basename(samFileName)
	if(options.outputfile):
		outputFileName = options.outputfile
		
	# write to the output file
	try:
		outfile = open(outputFileName, 'w')
	except IOError:
		print('error: write to output file failed!')
		sys.exit(-1)

	# Output header

	for header in sam.header :
		outfile.write(header.str().strip())

	# Output alignments

	linecount = 0
	for key in alignments :

		# Apply the unique strategy

		score = alignments[key]
		if(((score & 0xFFFF0000) >> 16) == 0):
			outfile.write(db[key])
			linecount += 1
	print('  %ld alignments written.' % (linecount))
	# Clear resources

	db.close()
	os.unlink(dbFileName)

	print('* Finished')


if __name__ == '__main__':
	main()