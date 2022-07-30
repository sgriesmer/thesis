import os
import re

seq_file = "file4.txt"
seq_dir = "RF00004_60_blastclust_sub1"
root = "c:/sjg/2010/thesis/"

# secondary structure pattern

ss_pattern = re.compile(r'^[\.\(\)]*$')


# print header

print("sample type,sample name,rnafold structure,rnafold base pairs,rnafold number of loops,rnafold largest loop,rnafold largest bulge,rnafold hairpin length,rnafold max consecutive base pairs,rnafold length,centroidfold structure,centroidfold base pairs,centroidfold number of loops,centroidfold largest loop,centroidfold largest bulge,centroidfold hairpin length,centroidfold max consecutive base pairs,centroidfold length\n")

# read aln file

for file in open(root + seq_dir + "/" + seq_file):
	file.strip()
	rnafold_file = file.replace("-fin","-rnafold")
	centroidfold_file = file.replace("-fin", "-centroidfold")
	combo = file.replace("-fin.txt","")

# run through rnaalifold

	rnafold_cmd = "rnaalifold <" + root + seq_dir + "/" + file
	rnafoldout = os.system(rnafold_cmd)
	break

for line in open("ttt.txt"):
	if ss_pattern.search(line) != None:
		print(line)


