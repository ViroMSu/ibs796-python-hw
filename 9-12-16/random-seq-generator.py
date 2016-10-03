#!/usr/bin/python
#chmod +x script.py to make executable
import sys
import random

#sys.argv is a list of commandline arguments passed to script
#index 0 is name of script
#example usage ./random-seq-generator.py TOTAL GC READ_LENGTH 
#Example output:
#        >sequence_001 gc=0.30 length=20
#       ATGCCATTATATGCCATTAT
#GC is given as decimal
num_reads = int(sys.argv[1])
gc = float(sys.argv[2])
read_len = int(sys.argv[3])

#for number of sequences, initiate sequence and header
for i in range(1, num_reads+1):
	seq = []
	line1 = []
	for c in '>sequence_':
		line1.append(c)

#for length of sequence, if random number is less than gc content, add either G or C
	for b in range(read_len):
		if random.random() < gc:
			seq.append(random.choice('GC'))
		else:
			seq.append(random.choice('AT'))

#Compute actual gc content plus lots of string pasting
	gc_frac = (seq.count("G") + seq.count("C"))/float(read_len)
	line1.extend(["%03d" % i, ' gc=', str(gc_frac), ' length=', str(read_len)])

#print to terminal
	print ''.join(line1)
	print ''.join(seq)

	if i == 1:	
		with open("random-seq-generator.output.txt", "w") as w:
			w.write(''.join(line1))
			w.write('\n')
			w.write(''.join(seq))
			w.write('\n')
	else:
		with open("random-seq-generator.output.txt", "a") as w:
			w.write(''.join(line1))
			w.write('\n')
			w.write(''.join(seq))
			w.write('\n')
