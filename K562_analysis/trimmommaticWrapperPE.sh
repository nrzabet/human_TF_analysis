#!/bin/bash
# File: trimmommaticPE.sh
#
#$ -wd /home/rtpop/ENCODE_ChIP_data/HepG2/pair_end/
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N TrimPEfiles


java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 1 \
"$1.fastq.gz"  "$2.fastq.gz" "$3.fastq.gz" "$4.fastq.gz" "$5.fastq.gz" "$6.fastq.gz" \
ILLUMINACLIP:/home/rtpop/Trimmomatic-0.39/adapters/adapters.fa:2:3:10:2:keepBothReads \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:"$7"
