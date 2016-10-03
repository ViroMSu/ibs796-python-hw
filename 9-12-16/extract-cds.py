#!/usr/bin/python

#Example usage ./extract-cds.py GENBANK_FILE

#Example output(multi-fasta)
#>SA_RS00145 product=chromosomal replication initiator protein DnaA
#MGDAVLDQYVRTYIVLKLKSKPNKLHQMSKKYVSAKSQAQTLEYLMEQEWFTDEEM


import sys
import collections
import re
import string

gbank = sys.argv[1]
with open(gbank, "r") as gb:
	lines = gb.readlines()

genes = collections.OrderedDict()
seq = []
chars = ['(', ')', '>', '<']
#Strip new line character and get lines with query, locus tags, products 
for i in range(len(lines)):
	line = lines[i].strip()
	if line.startswith('CDS') and not ':' in line and not 'pseudo' in lines[i-1].strip():
		line = line.translate(None, ''.join(chars))
		line = line.replace('..', ' ')
		elem = line.split()
		genes[i] = [elem[1], elem[2]]
		#get locus tag
		if 'locus' in lines[i+1]:		
			locus_line = lines[i+1].strip()
		elif 'locus' in lines[i+2]:
			locus_line = lines[i+2].strip()
		locus = re.findall(r'"([^"]*)"', locus_line)
		locus = ''.join(locus)
		genes[locus] = genes.pop(i)
	#get product tag
	elif line.startswith('/product'):
		if line.count('"') == 2:
			product = re.findall(r'"([^"]*)"', line)
		else:
			x = 1
			product = re.findall(r'"([^"]*)', line)
			z=i+1
			while x <2:
				line3 = lines[z].strip()
				if line3.count('"') == 1:
					product.append(re.findall(r'([^"]*)"', line3)[0])
					x +=1
				else:
					product.append(line3)
					z +=1
		if len(genes[locus]) == 2:
			genes[locus].append(product)
	elif line.startswith('ORIGIN'):
		seq_indx = i
		break

#Get entire sequence based on origin index
for p in lines[seq_indx+1:]:
	line2 = p.strip()
	for x in line2.split():
		if x.isalpha():
			for z in x:
				seq.append(z)

gen_code = {'TTT': 'F', 'TTC': 'F', 'TTA': 'F', 'TTG': 'F',
'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
'ATT':'I', 'ATC': 'I', 'ATA': 'I',
'ATG': 'M',
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
'TAT': 'Y', 'TAC': 'Y',
'TAA': '', 'TAG': '', 'TGA': '',
'CAT': 'H', 'CAC': 'H',
'CAA': 'Q', 'CAG': 'Q',
'AAT': 'N', 'AAC': 'N',
'AAA': 'K', 'AAG': 'K',
'GAT': 'D', 'GAC': 'D',
'GAA': 'E', 'GAG': 'E',
'TGT': 'C', 'TGC': 'C',
'TGG': 'W',
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'} 


#Extract sequence for each gene
for d in genes:
	if 'complement' in genes[d][0]:
		start = int(genes[d][0].replace('complement', '')) -1
		stop = int(genes[d][1])
		comp_seq = ''.join(seq[start:stop]).upper()
		gene_seq = comp_seq.translate(string.maketrans("ATCG", "TAGC"))[::-1]
		genes[d].append(gene_seq)
	else:
		start = int(genes[d][0]) -1
		stop = int(genes[d][1]) 
		gene_seq = ''.join(seq[start:stop]).upper()
		genes[d].append(gene_seq)

stops = ['TAA', 'TAG', 'TGA']

#Translate protein
for k in genes:
	protein = []
	dna = genes[k][3]
	for g in range(0,len(dna),3):
		if g == 0:
			codon = 'ATG'
		else:		
			codon = dna[g:g+3]
		protein.append(gen_code[codon])
		if codon in stops:
			break
	genes[k].append(protein)

for e in genes:
	print '>%s product=%s' % (e, genes[e][2])
	print ''.join(genes[e][4])

for index, e in enumerate(genes):
	if index == 0:
		with open("extract-cds.output.txt", "w") as w:
			w.write('>%s product=%s\n' % (e, genes[e][2]))
			w.write(''.join(genes[e][4]))
			w.write('\n')
	else:
		with open("extract-cds.output.txt", "a") as w:
			w.write('>%s product=%s\n' % (e, genes[e][2]))
			w.write(''.join(genes[e][4]))
			w.write('\n')

