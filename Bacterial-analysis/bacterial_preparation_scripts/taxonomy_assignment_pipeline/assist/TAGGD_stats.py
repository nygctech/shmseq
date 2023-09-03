#!/usr/bin/env python

'''
Run as this: python TAGGD_stats.py \
    *_results.tsv \
   /output_file_path/output_file_NAME \

'''
import sys
import os

# The input file with TAGGD resukts - tsv output from TAGGD 
input_file = sys.argv[1]

# Output file to write to, will be a fastq file
output = sys.argv[2]

# Get the directory to put the output file in and the file name
output_path = os.path.dirname(output)
output_name = os.path.basename(output)

output_file = os.path.join(output_path, output_name+'_TAGGD-stats.log') 
#################################################################
# try to parse taggd output results file

with open(input_file, "r" ) as mR1:
    # First top line
    header=mR1.readline()
    
    perfect_count = 0
    ambig_count = 0
    unmatched_count = 0
    
    # For every line in file except the first top line
    for line in mR1:
        
        if line.split('\t')[1].startswith('MATCHED_UNAMBIGUOUSLY'):
            ambig_count += 1
        
        elif line.split('\t')[1].startswith('MATCHED_PERFECTLY'):
            perfect_count += 1
            
        elif line.split('\t')[1].startswith('UNMATCHED'):
            unmatched_count += 1
    
    tot = ambig_count + perfect_count + unmatched_count
    
    perfect_count_per = perfect_count / tot
    ambig_count_per = ambig_count / tot
    unmatched_count_per = unmatched_count / tot

    print(tot)
    print(perfect_count_per)
    print(ambig_count_per)
    print(unmatched_count_per)
    
mR1.close()

# Open output file to write to
with open(output_file, 'w+') as f:
    print('MATCHED_PERFECTLY %:\t' + str(perfect_count_per), file= f)
    print('MATCHED_UNAMBIGUOUSLY %:\t' + str(ambig_count_per), file= f)
    print('UNMATCHED %:\t' + str(unmatched_count_per), file= f)
    
f.close()

