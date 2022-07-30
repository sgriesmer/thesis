#! /usr/bin/python

import os
import re

# parameters for run


slice_dir = "srprna_ali-fasta_60-blastclust-sub1-nh/negative"
seq_file = "rnazout-neg-file2.txt"
root = "c:/sjg/2009/personal/thesis"
seq_no = 0
ref_seq = 0
brucei = 1

# sequence pattern match

sequence_pattern = re.compile(r'^[GUAC_]*$')

print("seq_no,number_sequences,col_no,reading_direction,aspi,mfe_avg,mfe_consensus,energy_contrib,cov contrib,combo,z_score,sci,svm_decision_value,probability,prediction,sequence")

cell = 0
value_list = []
value_array = []
for line in open(root + "/" + slice_dir + "/" + seq_file):
	if line.__contains__(':') and not line.__contains__('WARNING'):
		cleanline = line.strip()
		entry = cleanline.split(':')
		value_list.append(entry[1].strip())
		cell += 1
	elif sequence_pattern.search(line.strip()) != None and line.strip() != "" and cell == 14:
		value_list.append(line.strip())
		value_array.append(value_list)
		cell = 0
		value_list = []
	
for (i,record) in enumerate(value_array,start = 1):
	print(i, ','.join(record), sep=',')
		
		