#!/bin/bash
# File: callpeaksWrapper.sh
#
#$ -wd /home/rtpop/ENCODE_ChIP_data/HepG2/single_end/controls/sam_files
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N callPeaksPE


macs2 callpeak -t "$1.sam" -c "$2.sam" -f SAM -g 2.9e9 \
      -n "$1-peaks" -B -q 0.05 2> "$1.log"
