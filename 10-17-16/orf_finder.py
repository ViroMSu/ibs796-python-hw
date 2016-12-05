#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import re
import timeit

fasta = sys.argv[1]

def my_function(fasta):
	record = SeqIO.read(fasta, "fasta")

	seq = record.seq
	cseq = seq.reverse_complement()

	pattern = re.compile(r'(?=(ATG(?:...)*?)(?=TAG|TGA|TAA))')
	orfs = pattern.findall(str(seq))
	orf_ind = pattern.finditer(str(seq))

	corfs = pattern.findall(str(cseq))
	corf_ind = pattern.finditer(str(cseq))

	for i, match in enumerate(orf_ind):
		seq_len = len(orfs[i]) + 3
		end = match.start() + seq_len
		stop = seq[end-3:end]
	#	print ">ORF %s, strand +, start:%d, end:%d" % (i, match.start(), end)
	#	print "%s%s" % (orfs[i], stop)

	for i, match in enumerate(corf_ind):
		seq_len = len(corfs[i]) + 3
		end = match.start() + seq_len
		stop = seq[end-3:end]
	#	print ">ORF %s, strand -, start:%d, end:%d" % (i, match.start(), end)
	#	print "%s%s" % (corfs[i], stop)


times = timeit.Timer(lambda: my_function(sys.argv[1])).repeat(repeat=100, number=1)
print(sum(times) / 100.0)
