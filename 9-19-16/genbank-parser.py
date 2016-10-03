#!/usr/bin/python

#Example usage ./extract-cds.py GENBANK_FILE DNA_OUTPUT AA_OUTPUT --ordered

import argparse
from Bio import SeqIO
from collections import OrderedDict

parser = argparse.ArgumentParser("commandline of genbank file, dna output, aa output, optional ordered out")
parser.add_argument("gb", help="genbank file")
parser.add_argument("dna_out", help="name of dna fasta output")
parser.add_argument("aa_out", help = "name of aa fasta output")
parser.add_argument("--ordered", help = "option to see if output should be ordered as in genbank", action="store_true")

args = parser.parse_args()

#Read in genbank file
record = SeqIO.read(args.gb, "genbank")

features = {}
features2 = OrderedDict()

#Check if ordered was set
if args.ordered:
	gene_dict = features2
else:
	gene_dict = features

#For each feature(CDS(+/- pseudo), tRNA, rRNA)
for i in record.features:
	if i.type == 'CDS' or i.type == 'tRNA' or i.type == 'rRNA':
		feat_type = i.type
		locus = i.qualifiers['locus_tag'][0]
		product = i.qualifiers['product'][0]
		dna_info = i.location
		if dna_info.strand == 1:
			seq = record.seq[dna_info.start:dna_info.end]
		else:
			seq = record.seq[dna_info.start:dna_info.end].reverse_complement()
		gene_dict[locus] = [product, seq]
		if feat_type == 'CDS' and not 'pseudo' in i.qualifiers:
			gen_code = i.qualifiers['transl_table'][0]
			protein = seq.translate(table=gen_code, cds=True)
			gene_dict[locus].append(protein)

with open(args.dna_out, "w") as wd:
	for z in gene_dict:
		wd.write("%s product=%s\n" %(z, gene_dict[z][0]))
		wd.write("%s\n" %(gene_dict[z][1]))

cds = [x for x in gene_dict if len(gene_dict[x]) ==3]
with open(args.aa_out, "w") as wa:
	for z in cds:
		wa.write("%s product=%s\n" %(z, gene_dict[z][0]))
		wa.write("%s\n" %(gene_dict[z][2]))

