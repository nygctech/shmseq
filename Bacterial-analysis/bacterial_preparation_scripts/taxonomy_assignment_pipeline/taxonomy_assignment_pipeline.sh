#!/bin/bash

### Taxonomy assignment pipeline ###
## Dependent on:
# TAGGD (https://github.com/JoelSjostrand/taggd.git))
# SAMTOOLS (http://www.htslib.org/)
# BWA (http://bio-bwa.sourceforge.net/bwa.shtml)
# Fastq-pair (https://github.com/linsalrob/fastq-pair.git)
# Kraken2 (https://ccb.jhu.edu/software/kraken2/)
# UMI tools (https://github.com/CGATOxford/UMI-tools.git)

################################################################################################################################################
################################################################################################################################################
# Assign shortest allowed read length after quality trimming
R_LENGTH=100

# Set FW and RV reads
# Example name of FW: 10015CN00_E2_R1_trimmed.fastq.gz
# Example name of RV: 10015CN00_E2_R2_trimmed.fastq.gz
FW=
RV= 

# Set experiment name
EXP=

# Output folder. Subfolders with the sample names will be created under this folder
OUTPUTFOLDER=

# Barcode identifier file
IDS=SHM-seq/id_file/10015_barcodes.txt

# Keep all temporary files (usually no): yes/no
KEEP_ALL_FILES=no

# Genome reference mouse BWA index
MOUSEREF=GCF_000001635.27_GRCm39_genomic.fna

################################################################################################################################################
################################################################################################################################################
source activate /ahg/regevdata/projects/ST_SpatialTranscriptomics/britta_bac/conda_env/umi-tool

# Start time of script
START=$(date +%s)

# Set output folder
OUTPUT=${OUTPUTFOLDER}/${EXP}

# Create sub output folder
mkdir -p $OUTPUT

# Create temp folder
mkdir -p $OUTPUT/tmp

echo "######### INPUT FILES #########" >> $OUTPUT/${EXP}_BACpipeline.log
echo $FW >> $OUTPUT/${EXP}_BACpipeline.log
echo $RV >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n" >> $OUTPUT/${EXP}_BACpipeline.log
echo "######### MOUSE REFERENCE #########" >> $OUTPUT/${EXP}_BACpipeline.log
echo $MOUSEREF >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n" >> $OUTPUT/${EXP}_BACpipeline.log

##############################################################
# Deal with compressed files
assist/decompression.py $FW $RV $OUTPUT/tmp

##############################################################
# Quality Filter --> use of ST-pipeline's quality filtering
echo "######### Quality filtering R2... #########" >> $OUTPUT/${EXP}_BACpipeline.log

# R_LENGTH = shortest read to save after trimming
python assist/filterInputReadsNEW.py $OUTPUT/tmp/R1_TMP.fq $OUTPUT/tmp/R2_TMP.fq $OUTPUT/${EXP} $R_LENGTH

# Write the quality-trim-log to pipeline log
sed -n '1,10 p' $OUTPUT/${EXP}_R2_quality_trimmed.log >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n" >> $OUTPUT/${EXP}_BACpipeline.log

# Make bam to sam
samtools view -h $OUTPUT/${EXP}_R2_quality_trimmed.bam > $OUTPUT/${EXP}_R2_quality_trimmed.sam

# From sam to fastq + adding the UMI to the fastq header after 'B3:Z:'
cat $OUTPUT/${EXP}_R2_quality_trimmed.sam | grep -v ^@ | awk '{$13=substr($13,6); print "@"$1"_B3:Z:"$13"\n"$10"\n+\n"$11}' > $OUTPUT/${EXP}_R2_quality_trimmed.fastq    

# To save space, rm
rm $OUTPUT/tmp/R2_TMP.fq
rm $OUTPUT/${EXP}_R2_quality_trimmed.sam
rm $OUTPUT/${EXP}_R2_quality_trimmed.bam

##############################################################
##############################################################
# Run R1 in taggd 
#--> *_matched.fastq 
#--> *_results.tsv with fq header and ambi/matched results
echo "######### TAGGD R1 #########" >> $OUTPUT/${EXP}_BACpipeline.log

# Constructing output file name
SLASH=$(echo "/")
OUTTAGGD="$OUTPUT$SLASH$EXP"

# Running taggd
/ahg/regevdata/projects/splotch_spatial_transcriptomics/britta_bac/conda_env/stenv3/bin/taggd_demultiplex.py $IDS $OUTPUT/tmp/R1_TMP.fq ${OUTTAGGD} --metric Hamming --overhang 0

# Take stats regarding TAGGD results and print to log file, output stats-file is $OUTPUT/${EXP}_TAGGD-stats.log
python assist/TAGGD_stats.py ${OUTTAGGD}_results.tsv $OUTPUT/${EXP} 

echo "Stats TAGGD: " >> $OUTPUT/${EXP}_BACpipeline.log
sed -n '1,3 p' $OUTPUT/${EXP}_TAGGD-stats.log >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n " >> $OUTPUT/${EXP}_BACpipeline.log

# To save space, rm R1_TMP
rm $OUTPUT/tmp/R1_TMP.fq

##############################################################
##############################################################
echo "######### BWA #########" >> $OUTPUT/${EXP}_BACpipeline.log

bwa mem $MOUSEREF $OUTPUT/${EXP}_R2_quality_trimmed.fastq > $OUTPUT/${EXP}.sam

# Gettting unmapped reads
samtools view -f 4 $OUTPUT/${EXP}.sam > $OUTPUT/${EXP}_unmapped.sam
samtools flagstat $OUTPUT/${EXP}.sam > $OUTPUT/${EXP}_bwaStats.txt

# Save the amount of reads that mapped to mouse
echo "Stats BWA: " >> $OUTPUT/${EXP}_BACpipeline.log
sed -n '5 p' $OUTPUT/${EXP}_bwaStats.txt >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n " >> $OUTPUT/${EXP}_BACpipeline.log

# From sam to fastq
cat $OUTPUT/${EXP}_unmapped.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $OUTPUT/${EXP}_R2_unaligned.fastq
tr '_' ' ' < $OUTPUT/${EXP}_R2_unaligned.fastq > $OUTPUT/${EXP}_tmp && mv $OUTPUT/${EXP}_tmp $OUTPUT/${EXP}_R2_unaligned.fastq

rm $OUTPUT/${EXP}.sam
rm $OUTPUT/${EXP}_unmapped.sam
rm $OUTPUT/${EXP}_bwaStats.txt
##############################################################
##############################################################
# Collect the R2 reads which belongs to the R1 reads with good barcode and add barcode  to R2
echo -e "######### Pair R1 and R2 ######### \n" >> $OUTPUT/${EXP}_BACpipeline.log

# First only keep R2 which has a corresponding perfectly matched barcode in R1
fastq_pair $OUTPUT/${EXP}_matched.fq $OUTPUT/${EXP}_R2_unaligned.fastq

echo -e "######### Add barcode to R2 ######### \n" >> $OUTPUT/${EXP}_BACpipeline.log

paste -d '~' <(cat $OUTPUT/${EXP}_R2_unaligned.fastq.paired.fq) <(cut -d ' ' -f 2-5 $OUTPUT/${EXP}_matched.fq.paired.fq) | perl -F'~' -lane 'push(@buffer, $F[0]); if($line == 0){@buffer[0] .= "$F[1]"}; if(($line == 3) && @buffer){print join("\n",@buffer); @buffer = ()}; $line = ($line+1) % 4;' > $OUTPUT/tmp/${EXP}_R2_unaligned_matchBarcodeTMP.fastq
sed 's/1:N:0/ 2:N:0/g' $OUTPUT/tmp/${EXP}_R2_unaligned_matchBarcodeTMP.fastq >> $OUTPUT/${EXP}_R2_unaligned_matchBarcodeTMP2.fastq   

# Move the UMI in B3 to last of line
awk 'BEGIN{FS=OFS=" "} {if ($1 ~ /^@/) {a=$2; for (i=2;i<NF; i++) $i=$(i+1); $NF=a; print} else {print $1}}' $OUTPUT/${EXP}_R2_unaligned_matchBarcodeTMP2.fastq > $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq

##############################################################
##############################################################
# Run *unaligned.fastq in Kraken2

echo -e "######### Kraken2 ######### \n" >> $OUTPUT/${EXP}_BACpipeline.log
echo "Number of input reads into Kraken2: " >> $OUTPUT/${EXP}_BACpipeline.log
grep -n '^@' $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq | wc -l >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n" >> $OUTPUT/${EXP}_BACpipeline.log

# Kraken2 custom database
KRAKENDB=files/kraken2_custom_DB

kraken2 \
   --db $KRAKENDB \
   $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq \
   --output $OUTPUT/${EXP}_kraken2_output.txt \
   --classified-out $OUTPUT/${EXP}_R2_TaxaClassified.fastq \
   --report $OUTPUT/${EXP}_kraken2_report.txt \
   --confidence 0.01 \
   --use-names 
 
# Kraken2 stats
echo 'Kraken2 DB: ' $KRAKENDB >> $OUTPUT/${EXP}_BACpipeline.log
echo "Reads assigned as Bacteria (%): " >> $OUTPUT/${EXP}_BACpipeline.log
grep 'Bacteria' $OUTPUT/${EXP}_kraken2_report.txt | head -n 1 | cut -f1 >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n" >> $OUTPUT/${EXP}_BACpipeline.log

# Grep the unaligned fastq headers with a perfectly matched barcode
grep '^@' $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq | cut -c 2- > $OUTPUT/${EXP}_headers_unaligned_matchBarcode.txt
tr ' ' '|' < $OUTPUT/${EXP}_headers_unaligned_matchBarcode.txt > $OUTPUT/${EXP}_headers_unaligned_matchBarcodePIPE.txt 

# Grep the fastq headers 
awk '{ if ($1 == "C") { print } }' $OUTPUT/${EXP}_kraken2_output.txt | awk -F'\t' -v OFS='|' '{print $2, $3}' > $OUTPUT/${EXP}_kraken2_output_fq.txt

# Join the two files to get a file with fastq headers with perfectly matched barcode + taxa classification in one file, unmatched fastq headers are removed
join -j 1 -t"|" <(sort -k1 $OUTPUT/${EXP}_headers_unaligned_matchBarcodePIPE.txt) <(sort -k1 $OUTPUT/${EXP}_kraken2_output_fq.txt) > $OUTPUT/${EXP}_headers_wBarcode_TaxaClass.txt

##############################################################
##############################################################
echo -e "######### DL model, UMI collapsing and printing to output ######### \n" >> $OUTPUT/${EXP}_BACpipeline.log

conda deactivate
source activate /ahg/regevdata/projects/splotch_spatial_transcriptomics/britta_bac/conda_env/keras

DLMODEL=files/model.h5
ENCODER=files/encoder.h5

python assist/DLmodel_and_UMItools_perspot.py \
   $OUTPUT/${EXP}_kraken2_report.txt \
   $OUTPUT/${EXP}_headers_wBarcode_TaxaClass.txt \
   $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq \
   $DLMODEL \
   $ENCODER
   
# UMI collapsing stats
echo 'DL model: ' $DLMODEL >> $OUTPUT/${EXP}_BACpipeline.log
sed -n '1,9 p' $OUTPUT/${EXP}UMIfilt_stats.tsv >> $OUTPUT/${EXP}_BACpipeline.log
echo -e "\n " >> $OUTPUT/${EXP}_BACpipeline.log


##############################################################
##############################################################

echo "######### Finishing up... #########" >> $OUTPUT/${EXP}_BACpipeline.log
gzip $OUTPUT/${EXP}_R2_unaligned.fastq # Unaligned R2 to mouse ref
gzip $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq # unaligned fastq with a perfectly matched barcode
gzip $OUTPUT/${EXP}_R2_TaxaClassified.fastq # fastq file with kraken2 assigned taxa

if [[ $KEEP_ALL_FILES == "no" ]] 
then
   # Compress used files
   rm $OUTPUT/${EXP}_R2_quality_trimmed.log # Remove quality trim log, already in pipeline log
   rm $OUTPUT/${EXP}_R2_quality_trimmed.fastq # Input file to HISAT2
   rm ${OUTTAGGD}_ambiguous.fq # Remove ambigious reads from TAGGD
   rm ${OUTTAGGD}_results.tsv # Results file from TAGGD
   rm ${OUTTAGGD}_unmatched.fq # Unmatched barcode in R1, output from TAGGD
   rm $OUTPUT/${EXP}_TAGGD-stats.log # TAGGD stats file
   rm ${OUTTAGGD}_matched.fq # matched barcode in R1, output from TAGGD
   rm $OUTPUT/${EXP}_*paired* # paired R1 and R2 --> to get the R2 which has a perfectly matched barcode in R1
   rm $OUTPUT/${EXP}_*single* # Singletons after pairing
   rm $OUTPUT/${EXP}_headers_unaligned_matchBarcode.txt # unaligned fastq headers with a perfectly matched barcode
   rm $OUTPUT/tmp/${EXP}_R2_unaligned_matchBarcodeTMP.fastq # Temp file with R2 with perfectly matched barcode
   rm $OUTPUT/${EXP}_headers_unaligned_matchBarcodePIPE.txt # unaligned fastq headers with a perfectly matched barcode, '|' separated
   rm $OUTPUT/${EXP}_kraken2_output.txt # kraken2 output as a list with fastq header and assigned taxa
   rm $OUTPUT/${EXP}_headers_wBarcode_TaxaClass.txt #This is now stored in a pickle file instead ###
   rm $OUTPUT/${EXP}_R2_unaligned_matchBarcodeTMP2.fastq # Temp file 2 with R2 with perfectly matched barcode
   rm $OUTPUT/${EXP}UMIfilt_stats.tsv # UMI collapsing stat, already in Bac log
   rm $OUTPUT/${EXP}_kraken2_output_fq.txt # Fastq headers list
   rm $OUTPUT/${EXP}_kraken2_report.txt

   mv $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq.gz $OUTPUT/tmp/ ###
   mv $OUTPUT/${EXP}_R2_TaxaClassified.fastq.gz $OUTPUT/tmp/
   mv $OUTPUT/${EXP}_R2_unaligned.fastq.gz $OUTPUT/tmp/

   rm -r $OUTPUT/tmp/
else
   mv $OUTPUT/${EXP}_R2_unaligned_matchBarcode.fastq.gz $OUTPUT/tmp/
fi

echo -e "###############################################################"
END=$(date +%s)
DIFF=$(( $END - $START ))
HOURS=$(expr $DIFF / 60)
echo "BAC-pipeline took $HOURS minutes" >> $OUTPUT/${EXP}_BACpipeline.log
done
