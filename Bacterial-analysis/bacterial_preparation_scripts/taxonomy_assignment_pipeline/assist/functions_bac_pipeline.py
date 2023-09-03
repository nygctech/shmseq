#!/usr/bin/env python

# Script which stores functions used in other python scripts.

from itertools import islice,tee
import io
from collections import defaultdict
import re
import numpy as np



def taxa_assignment(report):
    """
    Create a dictionary with the taxa order with taxID as key and taxa order as value
    Takes the kraken2 report as input and Returns the dictionary.
    """
    tax_order_dict = {}
    sp_tax_dict = {} # This is an extra dictionary with only species and theri corresponding NCBI taxID

    for index, row in report.iterrows():
        
        # Domain
        if row['torder'] == 'D':
            domain = row['sci-name'].lstrip()
            dvalue = ','.join([domain] + ['Unassigned'] * 6)
            tax_order_dict[row['ncbi taxID']] = dvalue
        # Subdomain
        elif row['torder'] == 'D1':
            tax_order_dict[row['ncbi taxID']] = dvalue

        # Phylum
        elif row['torder'] == 'P':
            phylum = row['sci-name'].lstrip()
            pvalue = ','.join([domain] + [phylum] + ['Unassigned'] * 5)
            tax_order_dict[row['ncbi taxID']] = pvalue

        elif ((row['torder'].startswith('P')) & (len(row['torder']) > 1)): # if x1 or x2 or xn
            tax_order_dict[row['ncbi taxID']] = pvalue

        # Class
        elif row['torder'] == 'C':
            Class = row['sci-name'].lstrip()
            cvalue = ','.join([domain] + [phylum] + [Class] + ['Unassigned'] * 4)
            tax_order_dict[row['ncbi taxID']] = cvalue

        # Order
        elif ((row['torder'] == 'O') | (row['torder'] == 'C1')): 
            order = row['sci-name'].lstrip()
            ovalue = ','.join([domain] + [phylum] + [Class] + [order] + ['Unassigned'] * 3)
            tax_order_dict[row['ncbi taxID']] = ovalue            

        # Family
        elif ((row['torder'] == 'F') | (row['torder'] == 'O1') | (row['torder'] == 'C2')):
            family = row['sci-name'].lstrip()
            fvalue = ','.join([domain] + [phylum] + [Class] + [order] + [family] + ['Unassigned'] * 2)
            tax_order_dict[row['ncbi taxID']] = fvalue

        # Genus
        elif ((row['torder'] == 'G') | (row['torder'] == 'F1') | (row['torder'] == 'O2') | (row['torder'] == 'C3')):
            genus = row['sci-name'].lstrip()
            gvalue = ','.join([domain] + [phylum] + [Class] + [order] + [family] + [genus] + ['Unassigned'] * 1)
            tax_order_dict[row['ncbi taxID']] = gvalue

        # species
        elif ((row['torder'] == 'S') | (row['torder'] == 'G1') | (row['torder'] == 'F2') | (row['torder'] == 'O3') | (row['torder'] == 'C4')):
            specie = row['sci-name'].lstrip()
            svalue = ','.join([domain] + [phylum] + [Class] + [order] + [family] + [genus] + [row['sci-name'].lstrip()])
            tax_order_dict[row['ncbi taxID']] = svalue
            sp_tax_dict[specie] = row['ncbi taxID']
        
        # Subspecies    
        elif row['torder'] == 'S1':
            tax_order_dict[row['ncbi taxID']] = svalue

    return tax_order_dict, sp_tax_dict
   
def tax_order_sp(row, d, d_sp):
    """
    Based on taxID, return the taxa order found in d.
    If not taxID found in d (it is a subspecies), look for the genus+species in the d_sp.
    Returns taxa order.
    """
    if row['taxa_orderTMP'] == 10090:
        return 'Animalia,Chordata,Mammalia,Rodentia,Muridae,Mus,Mus musculus'
    
    elif row['taxa_orderTMP'] == 131567: # Cellular organism
        return ','.join(['Cellular organism'] + ['Unassigned'] * 6)
    
    else:
        try:
            taxorder = d[row['taxa_orderTMP']]
            return taxorder
        except KeyError:
            sp = ' '.join([row['sp_short'][0], row['sp_short'][1]])
            
            try:        
                new_taxid = d_sp[sp]
                return d[new_taxid]
            except KeyError:
                return ','.join(['Unassigned'] * 7)
           
def parser(filename):
    """
    Generate a list of tuples (header, read)
    """
    fastq_parsed = {}
    try:

        with open(filename) as fq:
            header = next(fq)
            read = next(fq)
            fastq_parsed[header[1:-1].split(' ')[0]] = read[:-1] # we don't want new-line chars or @ sign in reads
            while True:
                next(fq) # skip + line
                next(fq) # skip qscore
                header = next(fq) # store next read header
                read = next(fq) # store next read seq
                fastq_parsed[header[1:-1].split(' ')[0]] = read[:-1]
    except:
        StopIteration # when we run out of reads, stop
    return fastq_parsed
   
def dna_encode_embedding_table(dna_input, name="dna_encode"):
    """
    DNA embedding.
    """

    embedding_values = np.zeros([len(dna_input), 5], np.float32)
    values = ("A", "C", "G", "T", "N")
    for j, b in enumerate(dna_input):
        if b in values:
            embedding_values[j, values.index(b)] = 1
    return embedding_values

def stack_padding(info_xy):
    # Stack reads into one tensor
    info_xy['one hot tensor'] = info_xy.apply(lambda row: dna_encode_embedding_table(row['read']), axis=1)
    X = np.array(info_xy['one hot tensor'].tolist())

    # Padding to the same sequence length
    masking_value = -1
    max_seq_len = max(len(x) for x in info_xy['one hot tensor'].tolist())
    N = X.shape[0]
    dimension = 5

    Xpad = np.full((N, max_seq_len, dimension), fill_value=masking_value)
    for s, x in enumerate(X):
        seq_len = x.shape[0]
        Xpad[s, 0:seq_len, :] = x
        
    return Xpad

def extract_taxa_info(info_xy, column_header, slide_type):
    y_taxaorder = info_xy[column_header].tolist()
    y_fastqH = info_xy.index.tolist()

    if slide_type == 'ST':
        y_umi = info_xy['UMI'].tolist()
        return y_taxaorder, y_fastqH, y_umi
    elif slide_type == 'QC':
        return y_taxaorder, y_fastqH

def predict_taxa_in_model(Xpad, model):
    predictions = model.predict(Xpad)
    
    return predictions


