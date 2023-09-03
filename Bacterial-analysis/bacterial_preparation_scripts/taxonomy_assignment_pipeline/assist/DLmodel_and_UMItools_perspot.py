#!/usr/bin/env python

"""
Script to create tsv file with taxa counts per spot
Used for when doing UMI collapsing PER spot AND taxa

Can only be used for ST slides

Run like this:
	python DLmodel_and_UMItools_perspot.py kraken2_report.txt *headers_wBarcode_TaxaClass.txt *R2_unaligned_matchBarcode.fastq DL-model DL-encoder

Both inputs are output files from taxonomy assignmnet pipeline with kraken2
   - kraken2_report.txt      
   - *headers_wBarcode_TaxaClass.txt
   - *R2_unaligned_matchBarcode.fastq

"""
import sys
import os
import numpy as np
import pandas as pd
import pickle
from functions_bac_pipeline import taxa_assignment, tax_order_sp, parser, dna_encode_embedding_table, stack_padding, predict_taxa_in_model, extract_taxa_info
from umi_tools import UMIClusterer
from collections import Counter
import tensorflow.keras as keras
from tensorflow.keras import Model
    
#########################################################################
#########################################################################
def umi_grouping(row, umi_groups):
    """
    Function to define which UMI has been collapsed to which UMI.
    
    Returns the collapsed UMI
    """
    d = {}
    
    umi = bytes(row['UMI'][5:], 'utf-8')

    for umi_group in umi_groups:
        top_umi = umi_group[0]
    
        for raw_umi in umi_group:
            d[raw_umi] = top_umi
        
    return d[umi]

#########################################################################
#########################################################################
# Load kraken2 report file

report_file = sys.argv[1]

report = pd.read_csv(report_file, sep='\t', header=None, names=['fragments covered by clade (%)', 'fragments covered by clade (num)', 'fragments assigned (num)', 'torder', 'ncbi taxID', 'sci-name'])

# Create a dictionary with the taxa order with taxID as key and taxa order as value + sp_dict = only taxID as value and specie name as key (for species only) 
tax_order_dict, sp_dict = taxa_assignment(report)

#########################################################################
#########################################################################
# Load fastq files

# read input fasta file
fastq = parser(sys.argv[3])

#########################################################################
#########################################################################
# Load fastq header - taxa text file
ft_file = sys.argv[2]
path = os.path.dirname(os.path.abspath(ft_file))
base = os.path.basename(os.path.abspath(ft_file)).split('_headers')[0]

ft = pd.read_csv(ft_file, sep='|', header = None, usecols=[0, 3, 4, 5, 6], names = ['fastq', 'x', 'y', 'UMI', 'taxa'])

ft.set_index('fastq', inplace=True)

# Remove fastq headers with short UMIs
ft = ft[ft['UMI'].str.len().gt(11)]

# Place whole taxa order per taxID
ft['taxa_orderTMP'] = ft['taxa'].str.split(' ').str[-1].str[:-1].astype(int)
ft['sp_short'] = ft['taxa'].str.split(' ').str[0:2]

ft['taxa_order'] = ft.apply(lambda row: tax_order_sp(row, tax_order_dict, sp_dict), axis = 1)

# To keep assignments on GENUS level
ft['to_TMP'] = ft['taxa_order'].str.split(',').str[5] + ',' + ft['taxa_order'].str.split(',').str[4] + ',' + ft['taxa_order'].str.split(',').str[3] + ',' + ft['taxa_order'].str.split(',').str[2] + ',' + ft['taxa_order'].str.split(',').str[1] + ',' + ft['taxa_order'].str.split(',').str[0]

# Drop columns
ft.drop(columns= ['taxa_orderTMP', 'sp_short', 'taxa_order'], inplace=True)
ft.rename(columns={'to_TMP':'taxa_order'}, inplace=True)

# Pair with fastq headers
ft['read'] = ft.index.map(fastq)

#########################################################################
#########################################################################
# Load DL model and encoder

model = keras.models.load_model(sys.argv[4])
encoder = pickle.load(open(sys.argv[5], 'rb'))

#########################################################################
#########################################################################
########### DL reassignment #############

unassigned_string = ','.join(['Unassigned'] * 6)

dl_l = []

info_coordXY = ft.loc[:,['x','y']]
info_coordXY.drop_duplicates(inplace=True)
info_coordXY['tuple'] = list(zip(info_coordXY['x'], info_coordXY['y']))

for tup in info_coordXY['tuple'].tolist():
    
    # Select spot      
    info_xy = ft[(ft['x'] == tup[0]) & (ft['y'] == tup[1])]
     
    if info_xy.shape[0] >0:
        Xpad = stack_padding(info_xy) # Stacking and padding/masking of reads
        y_taxaorder, fastqH, umi = extract_taxa_info(info_xy, 'taxa_order', 'ST') # Collect taxa info
        predictions = predict_taxa_in_model(Xpad, model) # Predict assignments using model

        rv_predictions = encoder.inverse_transform(predictions.argmax(axis=1)) # Predict taxa using encoder
        
        # Reassign taxa based on prediciton
        new_taxaorder = []
        
        for i, taxon in enumerate(y_taxaorder):
            if taxon.startswith('Unassigned'):
                new_taxa = rv_predictions[i]
                
                # If unassigned by DL model
                if new_taxa == '':
                    new_taxa = unassigned_string
                
                new_taxaorder.append(new_taxa)
            else:
                new_taxaorder.append(taxon)

        # Store in df
        p = pd.DataFrame({'fastq':fastqH, 'taxa_order':new_taxaorder, 'UMI':umi})
        # Add spot coord
        p['x'] = tup[0].split(':')[-1]
        p['y'] = tup[1].split(':')[-1] 

        dl_l.append(p)

df = pd.concat(dl_l)

# Some stats, before dropping unassigned
stat_file = os.path.join(path, base + "UMIfilt_stats.tsv")
with open(stat_file, 'a') as f:
    print('Number of reads BEFORE dropping unassigned: ' + str(df.shape[0]), file=f)

# Drop 'Unassigned' *6
df = df[df['taxa_order'] != unassigned_string]

with open(stat_file, 'a') as f:
    print('Number of reads AFTER dropping unassigned: ' + str(df.shape[0]), file=f)

#########################################################################
#########################################################################
##### Do UMI collapsing with 1 allowed mismatch PER spot AND taxa  ######

# Some stats - before UMI collapsing
with open(stat_file, 'a') as f:
    print('Number of reads BEFORE UMI collapsing: ' + str(df.shape[0]), file=f)

spot_dfs = []

for spot, group in df.groupby(['x', 'y']):
    
    # All UMIs per spot. Remove the "B3..."
    umi_list = [bytes(item[5:], 'utf-8') for item in group['UMI'].tolist()] 
 
    # First, count number of occurances for each UMI
    umi_raw_counts = dict(Counter(umi_list))

    # Run UMI-tools
    clusterer = UMIClusterer(cluster_method="directional")
    umi_groups = clusterer(umi_raw_counts, threshold=1)
    
    # Identify the collapsed UMI per UMI
    group['UMIcollapsed'] = group.apply(lambda row: umi_grouping(row, umi_groups), axis=1)
    
    # Remove rows where taxa and UMI are identical (per spot)
    group = group.drop_duplicates(subset=['taxa_order', 'UMIcollapsed'])

    # Add spots (x, y) to group
    group['x'] = spot[0]
    group['y'] = spot[1]

    spot_dfs.append(group)
    
dff = pd.concat(spot_dfs)
dff.set_index('fastq', inplace=True)
# Some stats after UMI collapsing
with open(stat_file, 'a') as f:
    print('Number of reads AFTER UMI collapsing: ' + str(dff.shape[0]), file=f)
    print('\n', file=f)
    print('Reads lost in UMI collapsing: ' + str((df.shape[0] - dff.shape[0])/df.shape[0]) + ' %', file=f)
    
#########################################################################
#########################################################################
# Create an output tsv file

# Spot coordinate xXy
dff['spot_coord'] = dff['x'].str.split(':').str[-1] + 'x' + dff['y'].str.split(':').str[-1]

# Count taxa per spot coordinate

dff.drop(columns=['x', 'y'], inplace=True)
dff.reset_index(inplace=True)
print(dff.head())

tsv = dff.groupby(['spot_coord','taxa_order']).count()[['fastq']].reset_index().pivot(index='spot_coord', columns='taxa_order', values='fastq')
tsv.fillna(0, inplace = True)

#print to output file
tsv_file = os.path.join(path, base + "_Bacstdata.tsv")

tsv.to_csv(tsv_file, sep = '\t')
