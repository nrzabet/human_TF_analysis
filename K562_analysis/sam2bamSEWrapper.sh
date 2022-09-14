#!/bin/bash
# File: sam2bamSEWrapper.sh
#
#$ -wd /home/rtpop/ENCODE_ChIP_data/HepG2/single_end/controls/
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N Sam2Bam


/home/rtpop/shell_scripts/sam2bamSEcmds.sh $1
