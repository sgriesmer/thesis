#! /usr/bin/python

###################################################################################################################################################
#
#
#   rnaz_score.py - creates RNAz output file from a list of sequence alignments
#   Input:      seq_file - the set of sequence alignments
#   Output:     output of RNAz command in a single file
#
#
###################################################################################################################################################

import os

# parameters for run

root = "c:/sjg/2010/thesis/"
seq_dir = "RF00004_60_blastclust_sub1"
seq_file = "file4.txt"

# create output filename

rnaz_file = "rnazout2-" + seq_file

# create output of rnaz scores and sequences

for line in open(root + seq_dir + "/" + seq_file):
	score = "rnaz" + " " + root + seq_dir + "/" + line.rstrip() + " >> " + root + seq_dir + "/" + rnaz_file
	print(score)
	os.system(score)
