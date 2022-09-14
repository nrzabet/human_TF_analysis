#!/bin/bash
# File: bowtieWrapperSE.sh
#
#$ -wd /home/rtpop/ENCODE_ChIP_data/HepG2/single_end/
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N alignSE

bowtie2 -x /home/rtpop/hg38/bwt_index/hg38_bwt_index -U "$1.fastq.gz" -S \
        "$1.sam" 2> "$2.log"
