#!/usr/bin/python

#Example Usage:
#        ./random-seq.py TOTAL GC READ_LENGTH --seed 12345
#Example output:
#        >sequence_001 gc=0.30 length=20
#        ATGCCATTATATGCCATTAT
#		 >sequence_002 gc=0.25 length=20
#        ATGCCATTATGAAAATTTCA


import argparse 
import random
from Bio.Seq import Seq
from Bio.SeqUtils import GC

parser = argparse.ArgumentParser(description="command line arguments of sequence info and seed")
parser.add_argument("total", help="total number of reads to generate", type=int)
parser.add_argument("gc", help="gc content to aim for in sequence", type=float)
parser.add_argument("length", help="length of sequences generated", type=int)
parser.add_argument("--seed", help = "random seed for reproduction", type=int)

args = parser.parse_args()

#Set random seed
if args.seed:
	random.seed(args.seed)

#Given gc content and length, generate random sequence
def seq_gen(gc, length):
	rseq = Seq('')
	for i in range(length):
		if random.random() < gc:
			rseq += Seq(random.choice('GC'))
		else:
			rseq += Seq(random.choice('AT'))
	return rseq

for i in range(1, args.total+1):
	rseq = seq_gen(args.gc, args.length)
	gc_actual = GC(rseq)/100
	header = ">sequence_%s gc=%s length=%s" %("%03d" % i, gc_actual, args.length)
	print header
	print str(rseq)

#Open file to write to
with open("random-seq.py.out", "w") as w:
	#For number of reads requested, call seq_gen function
	for i in range(1, args.total+1):
		rseq = seq_gen(args.gc, args.length)
		gc_actual = GC(rseq)/100
		header = ">sequence_%s gc=%s length=%s\n" %("%03d" % i, gc_actual, args.length)
		w.write(header)
		w.write(str(rseq))
		w.write('\n')


