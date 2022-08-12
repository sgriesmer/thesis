#! /usr/bin/python

###################################################################################################################################################
#
#   pair_prob.py - computes Shannon entropy and base-pair distance for a set of alignment files
#   Input:      seq_file - the set of sequence alignments
#   Output:     Shannon entropy and base-pair distance for each alignment on standard output
#
###################################################################################################################################################

import re
import math
import os

seq_file = "file2-abbrev.txt"
seq_dir = "RF00004_60_blastclust_sub1"
root = "c:/sjg/2010/thesis/"
sample_type = "positive"

# patterns

first_line_pattern = re.compile(r'sequence; length of alignment')
base_data_pattern = re.compile(r'^\s*[0-9]*\s*[0-9]')

print("sample type,sample name,shannon entropy,base-pair distance")

# read aln file

for file in open(root + seq_dir + "/" + seq_file):

# initialize sums
	sum_shannon = 0
	sum_dbp = 0
	combo = file.replace('-fin.txt','').strip()

# run through rnaalifold

	rnafold_cmd = "rnaalifold -p < " + root + seq_dir + "/" + file
	rnafoldout = os.system(rnafold_cmd)

# get alignment length and shannon and base pair distance

	for line in open(root + "alifold.out"):
     
# get alignment length

		if first_line_pattern.search(line.strip()) != None:
			alignment = line.split(';')
			length = int(alignment[1].replace('length of alignment ',''))
		elif base_data_pattern.search(line.strip()) != None:
			base_data = line.split()
			prob = float(base_data[3].replace('%', ''))
			prob = prob/100.0
   
# compute shannon probability and base pair distance

			if prob != 0:
				shannon = (-1)*prob*math.log(prob)/0.6931
			sum_shannon += shannon
			dbp = prob - pow(prob,2)
			sum_dbp += dbp
   
	norm_shannon = sum_shannon/length
	norm_dbp = sum_dbp/length
 
 # print output lines separated by commas
  
	out = [sample_type, combo, str(norm_shannon), str(norm_dbp)]
	print(','.join(out))
