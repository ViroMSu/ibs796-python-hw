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

import argparse
from Bio import SeqIO
import numpy as np

parser = argparse.ArgumentParser("commandline arguments of fastq and output file name")
parser.add_argument("fastq", help="name of FASTQ file")
parser.add_argument("--output", help="name of output file")

args = parser.parse_args()

#Get total reads and read lengths and quality scores
tot_reads = 0 
read_len = []
for rec in SeqIO.parse(args.fastq, "fastq"):
	tot_reads += 1
	read_len.append(len(rec))

#Statistics: total bp, mean read, min read, max read
read_distrib, read_freq = np.unique(read_len, return_counts=True)
tot_bp = sum(read_len)
min_read = min(read_len)
max_read = max(read_len)
mean_read = np.mean(read_len)

#Matrix for base qual stats
ncol, nrow = max_read, 2
qual_matrix = [[0 for x in range(ncol)] for y in range(nrow)]
avg_qual = []
for rec in SeqIO.parse(args.fastq, "fastq"):
	qual = rec.letter_annotations['phred_quality']
	avg_qual.append(np.mean(qual))
	for i, base in enumerate(qual):
		qual_matrix[0][i] += base
		qual_matrix[1][i] += 1

#Base qual stats
mean_qual = np.mean(avg_qual)
min_qual = min(avg_qual)
max_qual = max(avg_qual)

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

#if output file name given
if args.output:
	with open(args.output, "w") as w:
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


