#!/usr/bin/env python
"""
Run script:
   python decompression.py FW.fq RV.fq tempfolder
   
Takes compressed R1 and R2 fastq files as input and puts decompressed files in tempfolder. Doesn't affect the original Fw and RV.

"""

import sys
import os
import subprocess

FW_input = sys.argv[1]
RV_input = sys.argv[2]
tmp_folder = sys.argv[3]
##############################################################
# Deal with compressed input files
# Check if input fastq files are compressed
temp_r1_name = os.path.join(tmp_folder, "R1_TMP.fq")
temp_r2_name = os.path.join(tmp_folder, "R2_TMP.fq")
      
if FW_input.endswith(".gz"):
    r1_decompression_command = "gzip --decompress --stdout {} > {}".format(FW_input,temp_r1_name)
else:
    r1_decompression_command = None
    
if RV_input.endswith(".gz"):
    r2_decompression_command = "gzip --decompress --stdout {} > {}".format(RV_input,temp_r2_name)
else:
    r2_decompression_command = None

if r1_decompression_command:
    r1r = subprocess.Popen(r1_decompression_command, shell=True, preexec_fn=os.setsid)
if r2_decompression_command:
    r2r = subprocess.Popen(r2_decompression_command, shell=True, preexec_fn=os.setsid)
r1r.wait()
r2r.wait()
