#!/bin/bash

samtools view -q 1 -bt /home/rtpop/hg38/hg38.fa.gz "$1.sam" > "$1.bam"
samtools sort -l 4 -m 20G -O bam -T PREFIX "$1.bam" -o "$1.bam"
samtools index "$1.bam"
samtools idxstats "$1.bam" > "$1.tsv"
