#!/usr/bin/env python

"""
Script to make a simulated R2 by taking fastq headers, strand and quality score froma a real R2 but replacing the sequence with part of a 16S gene.
This simulated R2 is then used as a simulated fastq file with 16S reads in it. 
Outputs a simulated R2 fastq file and a txt file with the names (genus+species) of those whose 16S gene was collected (to have the truth of what species are actually in the simulated fastq file). 

From which species the 16S gene is taken from is chosen randomly based on a list of species for which there are 16S gene reference fasta available. The following has to be performed before running this script:
1. Select 16S gene fasta sequences (here: 16SMicrobial.fasta from NCBI) based on a list with interesting species. 
	--> use script prepare_custom_16S_ref_from_species_list.py
2. Run blastn on command line to align the 16S surface probe to these fasta references --> get a list of the E-values of the best alignments. 
3. Run this script to get a simluated R2 with 'fake' read 2 based on where the 16S surface probe aligned in the gene.  
	The newly created 'fake' R2 will contain the same number of fastq header and quality score as the original R2 which it is made from. 

Run script as:
	python prepare_sim_16S_read2.py E-value.txt output_folder name_on_output_files

"""
import os
import pandas as pd
import sys
import random
import numpy as np
import gzip
import shutil

fragment_length = 400 # Estimated fragment size of the 16S gene captured on the slide surface
max_R2_length = 300 # Max R2 length

######################################################################################################
######################################################################################################
# Read the fasta file with all the species.

fasta_file = '16SMicrobial.fasta'

fasta_dict = {}

with open(fasta_file) as fp:
    fasta = fp.readlines()
fp.close()

species = 'tmp'

for line in fasta:
    # Get speceis name
    if line.startswith('>'):
        if len(species.split(' ')) == 2:
            fasta_dict[species] = [long_species, ''.join(new_seq)]
            long_species = line.rstrip()
            species = ' '.join([line.split(' ')[1], line.split(' ')[2]])
            new_seq = []
        else:
            long_species = line.rstrip()
            species = ' '.join([line.split(' ')[1], line.split(' ')[2]])
            new_seq = []            
    # Get sequence
    else:
        new_seq.append(line.rstrip())
        
######################################################################################################
######################################################################################################
# Snippet to find the best query (16S_surface_probe) alignment in fasta reference and collect the fasta seq upstream of that alignment. 
# The e-value alignment file was created using blastn command line. 
# (q) = query
# (s) = subject, ie fasta reference

# read e-value alignment file from blastn command line
eval_file = sys.argv[1]

evals = pd.read_csv(eval_file, sep = '\t', comment='#', header=None, names = ['qacc', 'sacc', 'evalue', 'qstart', 'qend', 'sstart', 'send'])

# Get the sacc in a list 
sacc_list = list(set(evals['sacc'].tolist()))

# Get the start position in the fasta ref of the query alignment for each sacc
sacc_startpos_dict = {}

for sacc in sacc_list:
    
    # Sort to get the lowest evalue per sacc
    for index, row in evals.sort_values(['evalue']).groupby('sacc').head().iterrows():
        if sacc == row['sacc']:
            sacc_startpos_dict[sacc] = int(row['sstart']) # To include the 16S probe sequence which is on the surface 

# Create a dictionary with species ID as key and position in the fasta reference where the surface probes aligns as value
sp_read_length_dict = {}

for ki,vi in sacc_startpos_dict.items():
    
    for kj, vj in fasta_dict.items():
        
        if ki == vj[0].split('|')[3].split('.')[0]:
            
            sp_read_length_dict[kj] = ''.join(list(vj[1])[vi-1-fragment_length:vi-1])

######################################################################################################
######################################################################################################
# Replace the sequence in a R2 fastq file from a real run with these sequences to create a fake fastq file

real_sample_R2_file = 'S1_C2_S1_R2.fastq.gz'

# Gunzip input files
decompressed_real_sample_R2_file = real_sample_R2_file[:-3]

# Create decompressed file
with gzip.open(real_sample_R2_file, 'rb') as f_in:
    with open(decompressed_real_sample_R2_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# Store R2 file 
with open(decompressed_real_sample_R2_file) as r2:
    r2_lines = r2.readlines()

# Fastq header
r2_header_lines = r2_lines[0::4]

# Only look at the actual read
r2_read_lines = r2_lines[1::4]

# But also the quality line
r2_qual_lines = r2_lines[3::4]

os.remove(decompressed_real_sample_R2_file)

######################################################################################################
######################################################################################################
# Building a simulated R2
## For each fastq header I want to replace the sequence with randomly picked 16S gene
## R2 on this sequence also starts randomly somewhere on the sequence (this is dependent on the fragment length)

random_species_list = [] # To store what genus+species got picked

output_folder = sys.argv[2]
name = sys.argv[3]

r2_output_file_name = os.path.join(output_folder, name + '_simulated_R2.fastq')

with open(r2_output_file_name, 'a+') as f:

    for i, header in enumerate(r2_header_lines):
        
        # Randomly pick sequence from part of the 16S gene
        random_key = random.choice(list(sp_read_length_dict))
        random_value = sp_read_length_dict[random_key]
        
        # Store the genus+species which was randomlly selected
        random_species_list.append(random_key)
        
        # randomly select where the R2 read is gonna start from a normal distibuiton with mean = fragment_length
        # and standard deviation = stdev 
        stdev = 44
        # start position is end position (total length) - randomly picked fragment length using the normal distribution
        start = len(random_value) - int(np.random.normal(fragment_length, stdev))
        
        if start+max_R2_length > max_R2_length: # Fragment length was longer than the read length
            random_value_chopped = random_value[start:start+max_R2_length]
        elif start+max_R2_length < max_R2_length: # Fragment length was shorter than the read length read will go in 16S probe and further down towards the surface (read 2 includes the surface probe). However, in this simulation I will just shorten the read since I don't have anything beyond the end of the 16S probe.
            random_value_chopped = random_value[abs(start):abs(start)+max_R2_length] 
        elif start+max_R2_length == max_R2_length: # Fragment length equals read length
            random_value_chopped = random_value[start:start+max_R2_length]

        # Get the sequence length of the randomly picked key
        random_value_length = len(random_value_chopped)
        
        # Create a quality sequence as long as the randomly picked value, based on the existing quailty sequence
        if random_value_length == len(r2_qual_lines[i].rstrip()): # If same length as ramdonly picked sequence, do nothing
            qual_seq = r2_qual_lines[i].rstrip()
        elif random_value_length < len(r2_qual_lines[i].rstrip()): # If longer than ramdonly picked sequence, chop quality seq length
            qual_seq = r2_qual_lines[i].rstrip()[0:random_value_length]
        elif random_value_length > len(r2_qual_lines[i].rstrip()): # If shorter than ramdonly picked sequence, extend quality seq length
            extension = ''.join(random.choice(r2_qual_lines[i].rstrip()) for _ in range(random_value_length-len(r2_qual_lines[i].rstrip())))
            qual_seq = r2_qual_lines[i].rstrip() + extension
        
        # print to output file
        print(header.rstrip(), file= f) # fastq header
        print(random_value_chopped, file= f) # randomly picked sequence from part of the 16S gene
        print('+', file= f) # strand
        print(qual_seq, file= f) # quality sequence 

######################################################################################################
######################################################################################################
# Print the genus+species list to output file
gsp_output_file_name = os.path.join(output_folder, name + '_genus_species.txt')

with open(gsp_output_file_name, 'a+') as f:
    for item in random_species_list:
        print(item, file= f)
        



