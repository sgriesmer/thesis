###################################################################################################################################################
#
#
#   rna_structure.py - parses RNA alignment files for key secondary structural components (e.g., hairpin length, bulge number) using RNAalifold
#
#   Input:      seq_file - list of alignment files to process
#   Output:     key secondary structural components for each alignment in stdout
#
#
###################################################################################################################################################

import os
import re

# RNA alignment files and directories

seq_file = "file2-abbrev.txt"
seq_dir = "RF00004_60_blastclust_sub1"
root = "c:/sjg/2010/thesis/"
sample_type = "positive"

# secondary structure patterns

ss_pattern = re.compile(r'^[\.\(\)]*$')
sequence_pattern = re.compile(r'^[GUAC_]*$')

# print header

print("sample type,sample name,rnafold structure,rnafold base pairs,rnafold number of loops,rnafold largest loop,rnafold largest bulge,rnafold hairpin length,rnafold max consecutive base pairs,rnafold length\n")

# read list of alignment files

for file in open(root + seq_dir + "/" + seq_file):
    # file = file.strip()
    rnafold_file = file.replace("-fin","-rnafold").strip()
    sample_name = file.replace("-fin.txt","").strip()

# run alignment files through rnaalifold

    rnafold_cmd = "rnaalifold < " + root + seq_dir + "/" + file + " >" + root + seq_dir + "/" + rnafold_file + " 2>" + root + seq_dir + "/" + "err_file"
    rnafoldout = os.system(rnafold_cmd)

# extract secondary structure

    for record in open(root + seq_dir + "/" + rnafold_file):
        if sequence_pattern.search(record) == None:
            ss_struct = re.sub(r' .*$', '', record)

# parse secondary structure for key characteristics (e.g., number of hairpins, bulges)

    i = 0
    s0 = ""
    s1 = ""
    s2 = ""
    s3 = ""
    s4 = ""
    ss_length = 0
    largest_loop_length = 0
    largest_bulge_length = 0
    loop_number = 0
    bulge_number = 0
    blech_number = 0
    consecutive_base_pairs = 0
    max_consecutive_base_pairs = 0
    rnafold_hairpin = ""
    rnafold_hairpin_length = 0
    rnafold_length = 0

    # find structure length

    rnafold_length = len(ss_struct)

    # find hairpin length

    rnafold_hairpin = ss_struct
    rnafold_hairpin = re.sub(r'^\.*','',rnafold_hairpin)
    rnafold_hairpin = re.sub(r'\.*$','',rnafold_hairpin)
    rnafold_hairpin_length = len(rnafold_hairpin)

    # convert ss_struct to an array

    rnafold_array = list(ss_struct)

    current = rnafold_array[0]
    state = 's0'
    previous_state = 'start'

    while i < rnafold_length-1:
        if state == 's0' and current == '.':
            s0 += current
        elif state == 's0' and current == "(":
            s1 += current
            state = 's1'
            previous_state = 's0'
            consecutive_base_pairs += 1
        elif state == 's1' and current == "(":
            s1 += current
            consecutive_base_pairs += 1
        elif state == 's1' and current == '.':
            s2 += current
            state = 's2'
            previous_state = 's1'
            ss_length += 1
            if consecutive_base_pairs >= max_consecutive_base_pairs:
                max_consecutive_base_pairs = consecutive_base_pairs
            consecutive_base_pairs = 0
        elif state == 's2' and current == '.':
            s2 += current
            ss_length += 1
        elif state == 's2' and current == '(':
            s2 += current
            state = 's1'
            previous_state = 's2'
            if ss_length >= largest_bulge_length:
                largest_bulge_length = ss_length
            bulge_number += 1
            ss_length = 0
            consecutive_base_pairs += 1
        elif state == 's2' and current == ')':
            s3 += current
            state = 's3'
            previous_state = 's2'
            if ss_length >= largest_loop_length:
                largest_loop_length = ss_length
            loop_number += 1
            ss_length = 0
        elif state == 's3' and current == '.':
            s4 += current
            state = 's4'
            previous_state = 's3'
            ss_length += 1
        elif state == 's3' and current == '(':
            s1 += current
            state = 's1'
            previous_state = 's3'
            ss_length = 0
        elif state == 's3' and current == ')':
            s3 += current
        elif state == 's4' and current == '.':
            s4 += current
            ss_length += 1
        elif state == s4 and current == '(':
            s1 += current
            state = 's1'
            previous_state = 's4'
            blech_number += 1
            ss_length = 0
        elif state == 's4' and current == ')':
            s3 += current
            if ss_length != 0:
                if ss_length >= largest_bulge_length:
                    largest_bulge_length = ss_length
                bulge_number += 1
            ss_length = 0
        i += 1
        current = rnafold_array[i]

    bp = len(s1)

# print key characteristics from parsing
        
    outstring = [sample_type, sample_name, ss_struct.strip(), str(bp), str(loop_number), str(largest_loop_length), str(largest_bulge_length), str(rnafold_hairpin_length), str(max_consecutive_base_pairs), str(rnafold_length)]
    print(",".join(outstring))


