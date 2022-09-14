#!/bin/bash
# File: fastqc.sh
#
#$ -wd /home/rtpop/ENCODE_ChIP_data/HepG2/pair_end/
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N fastqcSE


fastqc "$1" -o fastqc
