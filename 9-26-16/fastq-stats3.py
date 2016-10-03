#!/usr/bin/python

#example usage ./fastq-stats.py FASTQ
#example output:
#Total Reads: 120000
#Total Base Pairs: 12234984bp
#Mean Read (Min, Max): 107.33 (33, 210)
#Read Lengths:
#    33:44
#    210: 20
#Mean Quality (Min, Max) Per Read: 36.67 (13, 38)
#Mean Quality Per Base:
#    1: 36
#    2: 37
#    210: 24

#Update with new functions
	#Convert FASTQ to FASTA as optional argument
	#removes reads with ambiguous nucleotides(N) by default
		#Add option to keep these reads
	#Removes low quality reads based on mean quality of the read as optional argument
	#Hard trim reads as an optional argument(left trim, right trim, both)
	#Trim reads based on per base quality
	#Output FASTQ/FASTA to file required positional argument
	#Output FASTQ stats before and after filtering to screen
		#optional argument to output to a file x

import argparse
from Bio import SeqIO
import numpy as np
import sys

parser = argparse.ArgumentParser("commandline arguments")
parser.add_argument("fastq", help="name of FASTQ file")
parser.add_argument("output", help="name of output file for fasta/q")
parser.add_argument("--statsout", help="name of output file for stats")
parser.add_argument("--convert", help="convert fastq to fasta", action="store_true")
parser.add_argument("--ltrim", help="number of characters to trim at beginning", type=int)
parser.add_argument("--rtrim", help="number of characters to trim at end", type=int)
parser.add_argument("--lrtrim", help="number of characters to trim from both ends", type=int)
parser.add_argument("--keepN", help="keep reads with ambiguious bases", action="store_true")
parser.add_argument("--meantrim", help="removes reads based on mean quality of read", type=int)
parser.add_argument("--basetrim", help="trim reads based on per base quality", type=int)

args = parser.parse_args()

if args.ltrim and args.lrtrim:
	sys.exit("Error: lrtrim cannot be called with ltrim or rtrim")
elif args.rtrim and args.lrtrim:
	sys.exit("Error: lrtrim cannot be called with ltrim or rtrim")

#Fastq before any processing
tot_reads = 0 
read_len = []
for rec in SeqIO.parse(args.fastq, "fastq"):
	tot_reads += 1
	read_len.append(len(rec))

read_distrib, read_freq = np.unique(read_len, return_counts=True)
tot_bp = sum(read_len)
min_read = min(read_len)
max_read = max(read_len)
mean_read = np.mean(read_len)

ncol, nrow = max_read, 2
qual_matrix = [[0 for x in range(ncol)] for y in range(nrow)]
avg_qual = []
for rec in SeqIO.parse(args.fastq, "fastq"):
	qual = rec.letter_annotations['phred_quality']
	avg_qual.append(np.mean(qual))
	for i, base in enumerate(qual):
		qual_matrix[0][i] += base
		qual_matrix[1][i] += 1

mean_qual = np.mean(avg_qual)
min_qual = min(avg_qual)
max_qual = max(avg_qual)

pbase_qual = []
for r in range(max_read):
	pbase_qual.append(qual_matrix[0][r]/float(qual_matrix[1][r]))

#Processing 
if args.basetrim and any(pbase_qual) < args.basetrim:
	pbase_array = np.array(pbase_qual)
	pbase_index = np.where(pbase_array <args.basetrim)[0]
	if len(pbase_index) == 0:
		args.basetrim = False

read_len2 = []
tot_qual = []
with open(args.output, "w") as wh:
	for rec in SeqIO.parse(args.fastq, "fastq"):
		seq = rec.seq
		qual = rec.letter_annotations['phred_quality']
		if args.keepN:
			pass
		elif 'N' in seq:
			seq = []
		if not len(seq) == 0:
			if args.ltrim:
				seq = seq[args.ltrim:len(seq)]
				qual = qual[args.ltrim:len(qual)]
			if args.rtrim:
				seq = seq[0:len(seq)-args.rtrim]
				qual = qual[0:len(qual)-args.rtrim]
			if args.lrtrim:
				seq = seq[args.lrtrim:len(seq)-args.lrtrim]
				qual = qual[args.lrtrim:len(qual)-args.lrtrim]
			if args.basetrim and ltrim and rtrim:
				trim_index = [i for i in pbase_index if i >= lrtrim][0]-ltrim
				if trim_index > len(seq):
					pass
				else:
					seq = seq[0:trim_index]
					qual = qual[0:trim_index]
			elif args.basetrim and lrtrim:
				trim_index = [i for i in pbase_index if i >= lrtrim][0]
				seq = seq[0:trim_index]
				qual = qual[0:trim_index]
			elif args.basetrim and ltrim:
				trim_index = [i for i in pbase_index if i >= ltrim][0]
				seq = seq[0:trim_index]
				qual = qual[0:trim_index]
			elif args.basetrim and rtrim:
				if trim_index > len(seq):
					pass
				else:
					seq = seq[0:trim_index]
					qual = qual[0:trim_index]
			if np.mean(qual) < args.meantrim:
				seq = []
		if not len(seq) == 0:
			if not args.convert:
				wh.write("@%s\n" % (str(rec.id)))
				wh.write("%s\n+\n" % (seq))
				qual2 = []
				for i in qual:
					qual2.append(chr(i+33))
				qual2 = ''.join(qual2)
				wh.write("%s\n" % (qual2))
			else: 
				wh.write(">%s\n" % (str(rec.id)))
				wh.write("%s\n" % (seq))
			read_len2.append(len(seq))
			tot_qual.append(qual)

ncol2, nrow2 = max(read_len2), 2
qual_matrix2 = [[0 for x in range(ncol2)] for y in range(nrow2)]
avg_qual2 = []
for qualz in tot_qual:
	avg_qual2.append(np.mean(qualz))
	for i, base in enumerate(qualz):
		qual_matrix2[0][i] += base
		qual_matrix2[1][i] += 1

read_distrib2, read_freq2 = np.unique(read_len2, return_counts=True)
	#First, check if there are Ns in sequence
	#Next, check if hard trims were set(left, right, both)
	#Next check if per base trim was set
	#Next, check mean quality and throw away those reads
	#Last, check if need to convert to fasta

#if output file name given
if args.statsout:
	with open(args.statsout, "w") as w:
		w.write("Stats pre-processing\n")
		w.write("Total Reads: %s\n" % (tot_reads))
		w.write("Total Base Pairs: %s\n" % (tot_bp))
		w.write("Mean Read (Min, Max): %s (%s, %s)\n" % (mean_read, min_read, max_read))
		w.write("Read Lengths:\n")
		for i, read in enumerate(read_distrib):
			w.write("\t%s: %s\n" % (read, read_freq[i]))
		w.write("Mean Quality (Min, Max) Per Read: %s (%s, %s)\n" % (mean_qual, min_qual, max_qual))
		w.write("Mean Quality Per Base:\n")
		for r in range(1, max_read+1):
			qual_pbase = qual_matrix[0][r-1]/float(qual_matrix[1][r-1])
			w.write("\t%s: %s\n" % (r, qual_pbase))
		w.write("Stats post-processing\n")
		w.write("Total Reads: %s\n" % (len(read_len2)))
		w.write("Total Base Pairs: %s\n" % (sum(read_len2)))
		w.write("Mean Read (Min, Max): %s (%s, %s)\n" % (np.mean(read_len2), min(read_len2), max(read_len2)))
		w.write("Read Lengths:\n")
		for i, read in enumerate(read_distrib2):
			w.write("\t%s: %s\n" % (read, read_freq2[i]))		
		w.write("Mean Quality (Min, Max) Per Read: %s (%s, %s)\n" % (np.mean(avg_qual2), min(avg_qual2), max(avg_qual2)))
		w.write("Mean Quality Per Base:\n")
		for r in range(1, ncol2+1):
			qual_pbase2 = qual_matrix2[0][r-1]/float(qual_matrix2[1][r-1])
			w.write("\t%s: %s\n" % (r, qual_pbase2))			
else: 
	print "Total Reads: %s" % (tot_reads)
	print "Total Base Pairs: %s" % (tot_bp)
	print "Mean Read (Min, Max): %s (%s, %s)" % (mean_read, min_read, max_read)
	print "Read Lengths:"
	for i, read in enumerate(read_distrib):
		print "\t%s: %s" % (read, read_freq[i])
	print "Mean Quality (Min, Max) Per Read: %s (%s, %s)" % (mean_qual, min_qual, max_qual)
	print "Mean Quality Per Base:"
	for r in range(1, max_read+1):
		qual_pbase = qual_matrix[0][r-1]/float(qual_matrix[1][r-1])
		print "\t%s: %s" % (r, qual_pbase)
	print "Stats post-processing\n"
	print "Total Reads: %s\n" % (len(read_len2))
	print "Total Base Pairs: %s\n" % (sum(read_len2))
	print "Mean Read (Min, Max): %s (%s, %s)\n" % (np.mean(read_len2), min(read_len2), max(read_len2))
	print "Read Lengths:\n"
	for i, read in enumerate(read_distrib2):
		print "\t%s: %s\n" % (read, read_freq2[i])	
	print "Mean Quality (Min, Max) Per Read: %s (%s, %s)\n" % (np.mean(avg_qual2), min(avg_qual2), max(avg_qual2))
	print "Mean Quality Per Base:\n"
	for r in range(1, ncol2+1):
		qual_pbase2 = qual_matrix2[0][r-1]/float(qual_matrix2[1][r-1])
		print "\t%s: %s\n" % (r, qual_pbase2)

