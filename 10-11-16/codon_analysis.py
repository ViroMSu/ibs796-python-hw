#!/usr/bin/python

import sys

fasta = sys.argv[1]

codontable = {
    'ATA':0, 'ATC':0, 'ATT':0, 'ATG':0,
    'ACA':0, 'ACC':0, 'ACG':0, 'ACT':0,
    'AAC':0, 'AAT':0, 'AAA':0, 'AAG':0,
    'AGC':0, 'AGT':0, 'AGA':0, 'AGG':0,
    'CTA':0, 'CTC':0, 'CTG':0, 'CTT':0,
    'CCA':0, 'CCC':0, 'CCG':0, 'CCT':0,
    'CAC':0, 'CAT':0, 'CAA':0, 'CAG':0,
    'CGA':0, 'CGC':0, 'CGG':0, 'CGT':0,
    'GTA':0, 'GTC':0, 'GTG':0, 'GTT':0,
    'GCA':0, 'GCC':0, 'GCG':0, 'GCT':0,
    'GAC':0, 'GAT':0, 'GAA':0, 'GAG':0,
    'GGA':0, 'GGC':0, 'GGG':0, 'GGT':0,
    'TCA':0, 'TCC':0, 'TCG':0, 'TCT':0,
    'TTC':0, 'TTT':0, 'TTA':0, 'TTG':0,
    'TAC':0, 'TAT':0, 'TAA':0, 'TAG':0,
    'TGC':0, 'TGT':0, 'TGA':0, 'TGG':0,
    }
tot_codon = 0 

with open(fasta, "r") as ff:
	for line in ff:
		if not line.startswith('>'):
			line = line.strip()
			tot_codon += len(line)/3
			for i in range(0,len(line),3):
				codon = line[i:i+3]
				if len(codon) == 3:
					codontable[codon] += 1

for q in codontable:
	print "%s: %.2f" % (q, codontable[q]/float(tot_codon))
