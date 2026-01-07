#!/bin/bash

usr="luis"
base_path="/home/${usr}/bioinf_isilon/Research/HALBRITTER/Internal/data/ML2cell/ATACseq_Seruggia_2024_2"
output_file="samplesheet_seruggia_2024-2-step1.csv"

# Create or overwrite the output file and add headers
echo "sample,fastq_1,fastq_2,replicate" > $output_file

# Loop through all combinations of X and Y
for X in {1..8}; do
for Y in {1..2}; do
sample="RBC-step-1-$X-r$Y"
fastq_1="$base_path/${sample}_R1_001.fastq.gz"
fastq_2="$base_path/${sample}_R2_001.fastq.gz"
replicate="1"
echo "$sample,$fastq_1,$fastq_2,$replicate" >> $output_file
done
done

output_file="samplesheet_seruggia_2024-2-step2.csv"
for X in {1..8}; do
for L in {A..H}; do
for Y in {1..2}; do

sample="RBC-step-2--$X$L-r$Y"   
fastq_1="$base_path/${sample}_R1_001.fastq.gz"
fastq_2="$base_path/${sample}_R2_001.fastq.gz"
replicate="1"
echo "$sample,$fastq_1,$fastq_2,$replicate" >> $output_file
done
done
done

cat samplesheet_seruggia_2024-2-step1.csv samplesheet_seruggia_2024-2-step2.csv >samplesheet_seruggia_2024_2_complete.csv
