#!/bin/bash
# File: bowtieWrapperPE.sh
#
#$ -wd /home/rtpop/ENCODE_ChIP_data/HepG2/pair_end/
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N alignPE

bowtie2 -x /home/rtpop/hg38/bwt_index/hg38_bwt_index -1 "$1.fastq.gz"\
        -2 "$2.fastq.gz" -S "$3.sam" 2> "$3.log"
