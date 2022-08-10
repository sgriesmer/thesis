#! /usr/bin/python

###################################################################################################################################################
#
#
#   cvt_aln_fasta.py - convert CLUSTAL alignment to FASTA file format to input to rnafold and calculate mean free energy ot he consensus structure
#
#   Input:      seq_file - list of alignment files to process in CLUSTAL format
#   Output:     mean free energy of the consensus structure output to standard output
#
#
###################################################################################################################################################


import re
from statistics import mean
import os

# patterns

blankline_pattern = re.compile(r'^\w*$')
similarity_line_pattern = re.compile(r'^\s\*')
sequence_line_pattern = re.compile(r'^[a-z]')
mfe_line_pattern = re.compile(r'^[\.\(\)]* ')

# open directory

seq_dir = "RF00004_60_blastclust_sub1"
neg_dir = "z_score"
err_file = "error"
seq_file = "file2a.txt"
root = "/home/azureuser/"
neg_path = root + seq_dir + "/" + neg_dir

# os.mkdir(neg_path)

for file in open(root + seq_dir + "/" + seq_file):
    
	# initialize arrays

	name = []
	d = {}
	sequence = ""

	# strip file entry of whitespace and nl

	file = file.strip()

	# Read Clustal file

	i = 0
	for line in open(root + seq_dir + "/" + file):
		line = line.strip()
		if blankline_pattern.search(line) != None:
			continue
		elif line.__contains__("CLUSTAL"):
			continue
		elif similarity_line_pattern.search(line) != None:
			continue
		elif sequence_line_pattern.search(line) != None:
			nm, sequence = line.split()
			nm = nm.strip()
			if nm not in name:
				name.append(nm)
				d[nm] = ""
			d[nm] += sequence

	fasta_file = file.replace('-fin','-fasta')
	mfe_file = file.replace('-fin', '-mfe')

	fasta_file_handle = open(root + seq_dir + "/" + fasta_file, 'w')

	for key,value in d.items():
		fasta_file_handle.write('>' + key +'\n')
		fasta_file_handle.write(value + '\n')

	cmd = 'RNAfold ' + '< ' + root + seq_dir + "/" + file + ' >' + root + seq_dir + "/" + mfe_file
	os.system(cmd)

	mfe_list = []

	for line in open(root + seq_dir + "/" + mfe_file):
		if mfe_line_pattern.search(line) != None:
			mfe = line[-9:]
			mfe = mfe.replace('(', '')
			mfe = mfe.replace(')', '')
			mfe = float(mfe.strip())
			mfe_list.append(mfe)
	mfe_mean = mean(mfe_list)
	print(",".join([mfe_file, str(mfe_mean)]))

	for f1 in [root + seq_dir + "/" + fasta_file, root + seq_dir + "/" + mfe_file]:
		os.remove(f1)