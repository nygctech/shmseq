#!/usr/bin/env python coding: utf-8

"""
This script is to train the DL model

Simulated reads here 150 bp since it reflects the actual read length in a real sample

Includes shorten of read length and random base error
"""

### set parameters
epochs = 10
batches = 10
tot_reads = 150000

# Sequence error rate 
errorR = 0.001

# To save model in folder
output_path = ''

# Simulated R2
input_fq = 'files/simulated_R2_150.fastq'

# Read in header ref file, ie fastq header is paired with taxa the sequence was made from
ref_input = 'files/header_taxa.txt'

# To save output figures
output_fig = ''

# Naming of output files
naming = str(tot_reads)+'_totreads_'+str(epochs)+'epochs_'+str(batches)+'batches_'+str(errorR)+'errorRate'

#####################################################################################################################################################################################
#####################################################################################################################################################################################

import pickle
import os
import numpy as np
import sys
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
import tensorflow.keras as keras
import tensorflow as tf
from tensorflow.keras.utils import to_categorical
from sklearn.metrics import roc_auc_score
from tensorflow.keras.callbacks import Callback
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report
from collections import Counter
import random
import sklearn.metrics
from numpy import argmax
import json
from functions_bac_pipeline import dna_encode_embedding_table, parser, tax_order_sp, taxa_assignment

#####################################################################################################################################################################################
#####################################################################################################################################################################################

class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

gpus = tf.config.list_physical_devices('cuda:0')
if gpus:
  # Restrict TensorFlow to only use the first GPU
  try:
    tf.config.experimental.set_visible_devices(gpus[0], 'GPU')
    logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPU")
  except RuntimeError as e:
    # Visible devices must be set before GPUs have been initialized
    print(e)

tf.test.is_gpu_available()

#####################################################################################################################################################################################
#####################################################################################################################################################################################

def shorten_reads(row):
    """
    Based on avergae read length in real ST experiment, shorten these simulated reads.
    """
    # Avergae read length + std of ST reads (sequenced 150bp)
    mu = 143
    sigma = 13
    random_length = int(np.random.normal(mu, sigma, 1)[0])

    # Read length cannot be more than 150 bp
    if random_length > 150:
        random_length = 150
    
    # New read length cannot be longer than what it already is
    if len(row['read']) <= random_length:
        new_length = len(row['read'])
    else:
        new_length = random_length
    
    new_read = row['read'][0:new_length]
    
    return new_read

def introducing_seq_error(row, num_errors):
    """
    Introducing a single base error (replacement). 
    """
    inds = [i for i,_ in enumerate(row) if not row.isspace()]

    sam = random.sample(inds, num_errors)

    lst = list(row)
    for ind in sam:
        letts =  iter(random.sample(['A', 'C', 'T', 'G'], 1))
        lst[ind] = next(letts)
    return "".join(lst)

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# # Load and prepare simulated data

# read input fasta file
fasta_d_tmp = parser(input_fq)

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# Read in header ref file, ie fastq header is paired with taxa the sequence was made from

ref_d = {}

with open(ref_input) as f:
    lines = f.readlines()
    
    for line in lines:
        header = line.split(' ')[0][1:]
        sp = line.split('|')[1].rstrip()
        
        ref_d[header] = sp

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# taxa names
ge_sp_dict = {'Blautia argi': ['Blautia argi', 'Blautia', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Anaerocolumna cellulosilytica': ['Anaerocolumna cellulosilytica', 'Anaerocolumna', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Akkermansia muciniphila': ['Akkermansia muciniphila', 'Akkermansia', 'Akkermansiaceae', 'Verrucomicrobiales', 'Verrucomicrobiae', 'Verrucomicrobia', 'Bacteria'], 
              '[Clostridium] scindens': ['[Clostridium] scindens', 'Lachnoclostridium', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Selenomonas ruminantium': ['Selenomonas ruminantium', 'Selenomonas', 'Selenomonadaceae', 'Selenomonadales', 'Negativicutes', 'Firmicutes', 'Bacteria'], 
              'Clostridium beijerinckii': ['Clostridium beijerinckii', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Parabacteroides distasonis': ['Parabacteroides distasonis', 'Parabacteroides', 'Tannerellaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              '[Clostridium] sphenoides': ['[Clostridium] sphenoides', 'Lacrimispora', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Flavonifractor plautii': ['Flavonifractor plautii', 'Flavonifractor', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Acetobacterium woodii': ['Acetobacterium woodii', 'Acetobacterium', 'Eubacteriaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Phocaeicola vulgatus': ['Phocaeicola vulgatus', 'Phocaeicola', 'Bacteroidaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Anaerobutyricum hallii': ['Anaerobutyricum hallii', 'Anaerobutyricum', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Enterococcus faecalis': ['Enterococcus faecalis', 'Enterococcus', 'Enterococcaceae', 'Lactobacillales', 'Bacilli', 'Firmicutes', 'Bacteria'], 
              'Hungatella hathewayi': ['Hungatella hathewayi', 'Hungatella', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium perfringens': ['Clostridium perfringens', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Muribaculum gordoncarteri': ['Muribaculum gordoncarteri', 'Muribaculum', 'Muribaculaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Butyrivibrio fibrisolvens': ['Butyrivibrio fibrisolvens', 'Butyrivibrio', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Roseburia hominis': ['Roseburia hominis', 'Roseburia', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Oscillibacter valericigenes': ['Oscillibacter valericigenes', 'Oscillibacter', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Roseburia intestinalis': ['Roseburia intestinalis', 'Roseburia', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Enterococcus faecium': ['Enterococcus faecium', 'Enterococcus', 'Enterococcaceae', 'Lactobacillales', 'Bacilli', 'Firmicutes', 'Bacteria'], 
              'Bacteroides thetaiotaomicron': ['Bacteroides thetaiotaomicron', 'Bacteroides', 'Bacteroidaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Blautia pseudococcoides': ['Blautia pseudococcoides', 'Blautia', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Massilistercora timonensis': ['Massilistercora timonensis', 'Massilistercora', 'Eubacteriales incertae sedis', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium chauvoei': ['Clostridium chauvoei', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Anaerostipes caccae': ['Anaerostipes caccae', 'Anaerostipes', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Duncaniella dubosii': ['Duncaniella dubosii', 'Duncaniella', 'Muribaculaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Blautia producta': ['Blautia producta', 'Blautia', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Ruminococcus bicirculans': ['Ruminococcus bicirculans', 'Ruminococcus', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Anaerostipes rhamnosivorans': ['Anaerostipes rhamnosivorans', 'Anaerostipes', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Pseudobutyrivibrio xylanivorans': ['Pseudobutyrivibrio xylanivorans', 'Pseudobutyrivibrio', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium baratii': ['Clostridium baratii', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Lachnospira eligens': ['Lachnospira eligens', 'Lachnospira', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium saccharoperbutylacetonicum': ['Clostridium saccharoperbutylacetonicum', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium septicum': ['Clostridium septicum', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium isatidis': ['Clostridium isatidis', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium gasigenes': ['Clostridium gasigenes', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Bacteroides ovatus': ['Bacteroides ovatus', 'Bacteroides', 'Bacteroidaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Bacteroides caecimuris': ['Bacteroides caecimuris', 'Bacteroides', 'Bacteroidaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Bacteroides xylanisolvens': ['Bacteroides xylanisolvens', 'Bacteroides', 'Bacteroidaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Sodaliphilus pleomorphus': ['Sodaliphilus pleomorphus', 'Sodaliphilus', 'Muribaculaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Lactobacillus johnsonii': ['Lactobacillus johnsonii', 'Lactobacillus', 'Lactobacillaceae', 'Lactobacillales', 'Bacilli', 'Firmicutes', 'Bacteria'], 
              '[Clostridium] hylemonae': ['[Clostridium] hylemonae', 'Lachnoclostridium', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Lacrimispora saccharolytica': ['Lacrimispora saccharolytica', 'Lacrimispora', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Lachnoclostridium phocaeense': ['Lachnoclostridium phocaeense', 'Lachnoclostridium', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Coprococcus comes': ['Coprococcus comes', 'Coprococcus', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium taeniosporum': ['Clostridium taeniosporum', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Blautia obeum': ['Blautia obeum', 'Blautia', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium butyricum': ['Clostridium butyricum', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Longicatena caecimuris': ['Longicatena caecimuris', 'Longicatena', 'Erysipelotrichaceae', 'Erysipelotrichales', 'Erysipelotrichia', 'Firmicutes', 'Bacteria'], 
              'Eubacterium limosum': ['Eubacterium limosum', 'Eubacterium', 'Eubacteriaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Muribaculum intestinale': ['Muribaculum intestinale', 'Muribaculum', 'Muribaculaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Barnesiella viscericola': ['Barnesiella viscericola', 'Barnesiella', 'Barnesiellaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Enterocloster bolteae': ['Enterocloster bolteae', 'Enterocloster', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Coprococcus catus': ['Coprococcus catus', 'Coprococcus', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridioides difficile': ['Clostridioides difficile', 'Clostridioides', 'Peptostreptococcaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Anaerocolumna sedimenticola': ['Anaerocolumna sedimenticola', 'Anaerocolumna', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Turicibacter sanguinis': ['Turicibacter sanguinis', 'Turicibacter', 'Turicibacteraceae', 'Erysipelotrichales', 'Erysipelotrichia', 'Firmicutes', 'Bacteria'], 
              'Bacteroides fragilis': ['Bacteroides fragilis', 'Bacteroides', 'Bacteroidaceae', 'Bacteroidales', 'Bacteroidia', 'Bacteroidetes', 'Bacteria'], 
              'Lachnoclostridium phytofermentans': ['Lachnoclostridium phytofermentans', 'Lachnoclostridium', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Butyrivibrio hungatei': ['Butyrivibrio hungatei', 'Butyrivibrio', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Anaerostipes hadrus': ['Anaerostipes hadrus', 'Anaerostipes', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Eubacterium maltosivorans': ['Eubacterium maltosivorans', 'Eubacterium', 'Eubacteriaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Intestinimonas butyriciproducens': ['Intestinimonas butyriciproducens', 'Intestinimonas', 'Eubacteriales incertae sedis', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium bornimense': ['Clostridium bornimense', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Enterocloster clostridioformis': ['Enterocloster clostridioformis', 'Enterocloster', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Clostridium botulinum': ['Clostridium botulinum', 'Clostridium', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Ruminococcus albus': ['Ruminococcus albus', 'Ruminococcus', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              '[Ruminococcus] torques': ['[Ruminococcus] torques', 'Mediterraneibacter', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Dysosmobacter welbionis': ['Dysosmobacter welbionis', 'Dysosmobacter', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Acutalibacter muris': ['Acutalibacter muris', 'Acutalibacter', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Faecalibacterium prausnitzii': ['Faecalibacterium prausnitzii', 'Faecalibacterium', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              '[Ruminococcus] gnavus': ['[Ruminococcus] gnavus', 'Mediterraneibacter', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Ruthenibacterium lactatiformans': ['Ruthenibacterium lactatiformans', 'Ruthenibacterium', 'Oscillospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'], 
              'Herbinix luporum': ['Herbinix luporum', 'Herbinix', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'],
              '[Clostridium] hathewayi': ['Hungatella hathewayi', 'Hungatella', 'Clostridiaceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'],
              '[Clostridium] clostridioforme': ['Enterocloster clostridioformis', 'Enterocloster', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria'],
              '[Eubacterium] hallii': ['Anaerobutyricum hallii', 'Anaerobutyricum', 'Lachnospiraceae', 'Eubacteriales', 'Clostridia', 'Firmicutes', 'Bacteria']}

# Add taxa order for the species names in dict

ref_df = pd.DataFrame.from_dict(ref_d, orient='index', columns=['taxa'])

ref_df['taxa_order'] = ref_df['taxa'].map(ge_sp_dict)
ref_df['taxa_order'] = ref_df['taxa_order'].str.join(',')

ref_d = ref_df.to_dict('index')

#####################################################################################################################################################################################
#####################################################################################################################################################################################
#### PREPARE SIMULATED READS ###

df = pd.DataFrame.from_dict(ref_d, orient='index', columns=['taxa', 'taxa_order'])
df.reset_index(inplace=True)
df = df.sample(frac=1).reset_index(drop=True) # Random shuffle

# Keep genera with at least x examples
MIN_COUNT = 100

df_keep = df.groupby('taxa').head(MIN_COUNT)
df_keep.reset_index(inplace=True)

# Randomly select the missing values
num_to_add = tot_reads - df_keep.shape[0]
print('Num to add:', num_to_add)

df_add = df[~df['index'].isin(df_keep['index'].tolist())]
df_add = df_add.sample(n = num_to_add)

# Merge df_keep and df_add
df_sim = pd.concat([df_keep, df_add])

# Match taxa with fastq header
df_sim['read'] = df_sim['index'].map(fasta_d_tmp)
df_sim.fillna('Unknown', inplace=True)

### Shorten read length based on real ST read lengths
df_sim['shortened read'] = df_sim.apply(lambda row: shorten_reads(row), axis=1) 

### INTRODUCE RANDOM ERRORS  ###

# Sum read length for all reads
tot_read_length = df_sim['shortened read'].str.len().sum()

tot_num_errors = int(tot_read_length * errorR)

if tot_num_errors > df_sim.shape[0]:
    # Randomly select reads which shall include a simulated seq error
    df_subset = df_sim.copy()
    num_errors = int(tot_num_errors/df_sim.shape[0])
else:
    # Randomly select reads which shall include a simulated seq error, to simplify: one error per read
    df_subset = df_sim.sample(n=tot_num_errors, replace=False)
    num_errors = 1

df_subset['shortened read'] = df_subset.apply(lambda row: introducing_seq_error(row['shortened read'], num_errors), axis=1)

# Delete subset from main df
new_sim = df_sim[~df_sim['index'].isin(df_subset['index'].tolist())]

df_sim = pd.concat([new_sim, df_subset])

# # SELECTINGA ##
merging_cols = ['index', 'taxa_order', 'shortened read'] 

df_merge = df_sim.loc[:,merging_cols]
df_merge['taxa_order_GE'] = df_merge['taxa_order'].str.split(',').str[1:]
df_merge['taxa_order_GE'] = df_merge['taxa_order_GE'].str.join(',')

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# # MODEL

# To save model/embedder
output_file = 'model_'+naming+'.h5'
encoder_file = 'encoder_'+naming+'.h5'

output_model = os.path.join(output_path, output_file)
output_encoder = os.path.join(output_path, encoder_file)

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# Stack into one tensor
df_merge['one hot tensor'] = df_merge.apply(lambda row: dna_encode_embedding_table(row['shortened read']), axis=1)
X = np.array(df_merge['one hot tensor'].tolist())

# Padding to the same sequence length
masking_value = -1
max_seq_len = max(len(x) for x in df_merge['one hot tensor'].tolist())
N = X.shape[0]
dimension = 5

Xpad = np.full((N, max_seq_len, dimension), fill_value=masking_value)
for s, x in enumerate(X):
    seq_len = x.shape[0]
    Xpad[s, 0:seq_len, :] = x

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# ONE-HOT ENCOIDNG OF THE TAXA ##

y = df_merge['taxa_order_GE']
encoder = LabelEncoder()
encoder.fit(y)
encoded_y = encoder.transform(y)

# convert integers to dummy variables (i.e. one hot encoded)
dummy_y = to_categorical(encoded_y)

# Save Encoder
pickle.dump(encoder, open(output_encoder, 'wb'))

assert set(dummy_y.sum(axis=1)) == {1}

#####################################################################################################################################################################################
#####################################################################################################################################################################################

from tensorflow.keras import Sequential, Model, callbacks
from tensorflow.keras.layers import Dense, LSTM, Dropout, Conv1D, SimpleRNN, Input, Concatenate, Bidirectional, Masking
from tensorflow.keras.callbacks import ModelCheckpoint, Callback
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from tensorflow.keras.optimizers import SGD

class FragmentClassifier(object):

    def __init__(self, input_dna, input_taxa, test_fraction,
                 epochs, batch_size, n_classes, Max_seq_len, Masking_value):
        self.input_dna = input_dna
        self.input_taxa = input_taxa
        self.test_fraction = test_fraction
        self.epochs = epochs
        self.batch_size = batch_size
        self.n_classes = n_classes
        self.max_seq_len = Max_seq_len
        self.masking_value = Masking_value
        self.x_train, self.x_test, self.y_train, self.y_test = self.load_and_split()
        
        # Initialize architecture
        self.classifier = None
        self.embedder = None
        self.architecture()
        
    def load_and_split(self):
        dna_reads = self.input_dna
        print('total number of reads in input:', len(dna_reads))
        print('number of classes:', self.input_taxa.shape[1])

        x = np.array(dna_reads)
        y = np.array(self.input_taxa)

        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=self.test_fraction)

        print(x_train.shape, x_test.shape, y_train.shape, y_test.shape)
        print("Training classes:", y_train.sum(axis=0).astype(int))
        print("Test classes:", y_test.sum(axis=0).astype(int))
               
        return x_train, x_test, y_train, y_test

    def architecture(self):
        
        inp = Input(shape=(None , 5))        
        mask = Masking(mask_value=self.masking_value, input_shape=(self.max_seq_len, 5))(inp)
        
        conv5 = Conv1D(filters=64, kernel_size=15, activation='relu', padding="same")(mask) 
        conv7 = Conv1D(filters=64, kernel_size=17, activation='relu', padding="same")(mask) 
        conv11 = Conv1D(filters=64, kernel_size=19, activation='relu', padding="same")(mask)
        conv13 = Conv1D(filters=64, kernel_size=23, activation='relu', padding="same")(mask)
        mrg = Concatenate()([conv5, conv7, conv11, conv13])
        mrg = Dropout(rate=0.5)(mrg)

        # Add 2 bidirectional LSTMs
        x = Bidirectional(LSTM(64, return_sequences=True))(mrg)
        x = Bidirectional(LSTM(64))(x)
        dp1 = Dropout(rate=0.2)(x)
        
        emb1 = Dense(64, activation='relu')(dp1) #16
        dp2 = Dropout(rate=0.1)(emb1)
        emb2 = Dense(32, activation='relu')(dp2) #8
        decision = Dense(self.input_taxa.shape[1], activation='softmax')(emb2)
        
        model = Model(inputs=inp, outputs=decision)
        model.compile(optimizer='Adam', 
                      loss='categorical_crossentropy', 
                      metrics=['categorical_accuracy']) 
                
        ### Point class to models
        self.embedder = embedder
        self.classifier = model

    def train(self):
        classifier = self.classifier
        embedder = self.embedder
        x_train, y_train = self.x_train, self.y_train
        x_test, y_test = self.x_test, self.y_test
        
        train_loss = []
        val_loss = []
        train_acc = []
        val_acc = []

        # Optimize number of epochs
        earlystopping = callbacks.EarlyStopping(monitor ="loss", 
                                        mode ="min", patience = 5, 
                                        restore_best_weights = True, verbose=1)
        
        # Define pairs in the training set
        training_set  = list(zip(x_train, y_train))
        
        num_training_loops = self.epochs
        for e in range(num_training_loops):
            
            ### it is a good idea to shuffle the training set before each training loop
            np.random.shuffle(training_set)
            
            # Compute derivatives for training set once
            for x, y in training_set:
                xt = x.reshape((1, x.shape[0], x.shape[1]))   
                yt = y.reshape((1, y.shape[0]))

                ### Do one update step with one sequence
                ca = classifier.fit(xt, yt, batch_size=self.batch_size, 
                                    epochs=1, verbose=0, validation_data=None, callbacks =[earlystopping]) 
                    
            ### evaluate once per epoch
            ### Computes metrics in the order they are entered in Model() above
            eval_tr = classifier.evaluate(x_train, y_train, verbose=0)
            eval_te = classifier.evaluate(x_test, y_test, verbose=0)

            
            ### Compute once per epoch ####
            ### Store only the last point ####
            train_loss.append(eval_tr[0])
            train_acc.append(eval_tr[1])
            val_loss.append(eval_te[0]) 
            val_acc.append(eval_te[1]) 
            
            ### Report progress
            print("Computed epoch %d/%d, training loss: %.3f, validation loss: %.3f, training acc: %.3f, validation acc: %.3f" % 
                  (e, num_training_loops, train_loss[-1], val_loss[-1], train_acc[-1], val_acc[-1]))
        
        ########### LOSS vs EPOCHS ################
        fig, ax = plt.subplots()
        plt.plot(train_loss)
        plt.plot(val_loss)
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'validation'])
        plt.show()
        fig.savefig(os.path.join(output_fig, 'loss_vs_epochs_'+naming+'.pdf'))
        
        ########### Acc vs EPOCHS ################
        fig, ax = plt.subplots()
        plt.plot(train_acc)
        plt.plot(val_acc)
        plt.title('model accuracy')
        plt.ylabel('acc')
        plt.xlabel('epoch')
        plt.legend(['train', 'validation'])
        plt.show()   
        fig.savefig(os.path.join(output_fig, 'acc_vs_epochs_'+naming+'.pdf'))
        
        ########### Save model and architecture to single file  ###########
        classifier.save(output_model)

#####################################################################################################################################################################################
#####################################################################################################################################################################################

# Separate plots from model.

def evaluation(model):
    classifier = model.classifier
    embedder = model.embedder
    x_test, y_test = model.x_test, model.y_test
    n_classes = model.n_classes

    y_pred_keras = classifier.predict(x_test)

    ########### ROC ################

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_pred_keras[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_pred_keras.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # Plot ROC 
    fig, ax = plt.subplots()
    ax.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]), color='deeppink', linestyle=':', linewidth=4)

    ax.plot([0, 1], [0, 1], 'k--', lw=2)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC-AUC')
    #ax.legend(loc="lower right")
    plt.show()
    fig.savefig(os.path.join(output_fig, 'ROC_AUC_'+naming+'.pdf'))

    ########### STATS ################

    test_F1 = sklearn.metrics.f1_score(np.argmax(y_test, axis=1), np.argmax(y_pred, axis=1), average='micro') 

    acc = np.mean(np.argmax(y_test, axis=1) == np.argmax(y_pred, axis=1))

    with open(os.path.join(output_fig, "output_metrics_"+naming+'.txt'), "a") as f:
        print("F1: " + str(test_F1), file=f)
        print("Accuracy: " + str(acc), file=f)

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# Initialize model and printout summary.

# Number of taxa to classify
n_classes = dummy_y.shape[1]

model = FragmentClassifier(Xpad, dummy_y, test_fraction=0.2, 
                           epochs=epochs, batch_size=batches, n_classes=n_classes, Max_seq_len=max_seq_len, Masking_value=masking_value)
model.classifier.summary()


#####################################################################################################################################################################################
#####################################################################################################################################################################################
# Train model

model.train()

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# Plot evaluation metrics.

evaluation(model)


