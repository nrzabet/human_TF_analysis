#!/bin/bash
# File: trimmommatic_SE.sh
#
#$ -wd /home/rtpop/ENCODE_ChIP_data/HepG2/single_end/
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N trimSE


java -jar /home/rtpop/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 \
"$1.fastq.gz" "$2.fastq.gz" \
ILLUMINACLIP:/home/rtpop/Trimmomatic-0.39/adapters/adapters.fa:2:3:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:"$3"
