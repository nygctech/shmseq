#!/usr/bin/env python

"""
Script to collect the full 16S sequence from a fasta reference file (used here: 16SMicrobial.fasta downloaded from NCBI https://ftp.ncbi.nlm.nih.gov/blast/db/ on 2021-07-24) for
species which are named 'genus specie' in a list (65sp.txt).

Run script as:
	python prepare_custom_16S_ref_from_species_list.py input_fasta.fasta sp_list.txt output_file.fasta


"""
import os
import pandas as pd
import sys

######################################################################################################
# Read the fasta file with all the species.

fasta_file = sys.argv[1]

fasta_dict = {}

with open(fasta_file) as fp:
    fasta = fp.readlines()
fp.close()

species = 'tmp'


### When using 16SMicrobial.fasta as reference
for line in fasta:
    # Get speceis name
    if line.startswith('>'):
        if len(species.split(' ')) == 2:
            fasta_dict[species] = [long_species, ''.join(new_seq)]
            long_species = line.rstrip()
            try:
                species = ' '.join([line.split('| ')[1].split(' ')[0], line.split('| ')[1].split(' ')[1]])
            except IndexError:
                print(line)
            new_seq = []
        else:            
            long_species = line.rstrip()
            species = ' '.join([line.split('| ')[1].split(' ')[0], line.split('| ')[1].split(' ')[1]])
            new_seq = []
    # Get sequence
    else:
        new_seq.append(line.rstrip())

######################################################################################################
# read in species list

species_file = sys.argv[2]

with open(species_file) as fp:
    sptmp = fp.readlines()
fp.close()

# Remove newline at end of each line
sps = [line.rstrip('\n') for line in sptmp]

######################################################################################################
output_file = sys.argv[3]
missing = []
with open(output_file, "a") as f:
    
    for sp in sps:
        if sp in fasta_dict:
            print(fasta_dict[sp][0], file=f)
            print(fasta_dict[sp][1], file=f)
        else:
            missing.append(sp)

print('Missing 16S sequences')
print(missing)
