#!/usr/bin/python

#gunzip fastq.gz
#example usage ./fastq-stats.py FASTQ

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

import sys
from collections import Counter

fastq = sys.argv[1]
with open(fastq, "r") as fq:
	lines = fq.readlines()

tot_reads = 0
tot_bp = 0
read_len = []
header = []	
qual = {}
mean_qual = []
for i in range(len(lines)):
	lines[i] = lines[i].rstrip('\n')

for i in range(len(lines)):
	#Name of read
	if lines[i].startswith('@') and i < len(lines)-2 and lines[i+2] == '+':
		tot_reads += 1
		header.append(i)
	#Sequence
	elif lines[i].isalpha() and i-1 in header:
		read_len.append(len(lines[i]))
		tot_bp += len(lines[i])
	#Quality Phred +33 offset
	elif i-3 in header:
		qual[i] = []
		#Dictionary read: quality scores
		for x in lines[i]:
			qual[i].append(ord(x)-33)
		mean_qual.append(sum(qual[i])/float(len(qual[i])))	

min_read = min(read_len)
max_read = max(read_len)
avg_read = sum(read_len)/float(len(read_len))

min_qual = min(mean_qual)
max_qual = max(mean_qual)
avg_qual = sum(mean_qual)/float(len(mean_qual))

#For distribution of read lengths
read_tbl = sorted(Counter(read_len).items())	

base_qual = {}
#For mean quality per base, dictionary for each base
for t in range(1, max_read+1):
	base_qual[t] = []

for v in qual.values():
	for q in range(1, len(v)+1):
		base_qual[q].append(v[q-1])

qual_pbase = []
for z in base_qual.values():
	qual_pbase.append(sum(z)/float(len(z)))

#Output
print "Total Reads: %s" % (tot_reads)
print "Total Base Pairs: %s" % (tot_bp)
print "Mean Read (Min, Max): %s (%s, %s)" % (avg_read, min_read, max_read)
print "Read Lengths:"
for w in read_tbl:
	print "\t%s: %s" %(w[0], w[1])
print "Mean Quality (Min, Max) Per Read: %s (%s, %s)" % (avg_qual, min_qual, max_qual)
print "Mean Quality Per Base:"
for r in range(1, len(qual_pbase)+1):
	print "\t%s: %s" % (r, qual_pbase[r-1])

with open("fastq-stats.output.txt", "w") as w:
	w.write("Total Reads: %s\n" % (tot_reads))
	w.write("Total Base Pairs: %s\n" % (tot_bp))
	w.write("Mean Read (Min, Max): %s (%s, %s)\n" % (avg_read, min_read, max_read))
	w.write("Read Lengths:\n")
	for t in read_tbl:
		w.write("\t%s: %s\n" %(t[0], t[1]))
	w.write("Mean Quality (Min, Max) Per Read: %s (%s, %s)\n" % (avg_qual, min_qual, max_qual))
	w.write("Mean Quality Per Base:\n")
	for r in range(1, len(qual_pbase)+1):
		w.write("\t%s: %s\n" % (r, qual_pbase[r-1]))
