####################################################################
################## Execute ChIPanalyser on IMR90 ###################
####################################################################

# Move 5-ChIPanalyser_general.R script to ChIPanalyser folder and get rid of "5-" to keep its name clean.
screen, cd /, cd storage/st20d/HiC_Kc167/alessandra/IMR90/ChIPanalyser/
chmod 770 ChIPanalyser_general.R

# If I have an old version of ChIPanalyser, download .tar.gz file from https://blog.revolutionanalytics.com/2020/09/mro-402-available.html, 
# move it to wd /ChIPanalyser/ in the cluster, open R3.6.0 and run:
install.packages("ChIPanalyser_1.12.0.tar.gz", repos=NULL) # then quit R and continue from bash. 

# 1st argument: name of tf (e.g. CTCF)
# 2nd argument: QDA ID (e.g. 12)
# 3rd argument: name of cell line (e.g. IMR90)
# 4th argument: regions_to_train_START (these are the regions to train)
# 5th argument: regions_to_train_END
# 6th argument: regions_to_test_START (these are regions to validate)
# 7th argument: regions_to_test_END
# 8th argument: validation_regions
# 9th argument: DNA accessibility method
# EXAMPLE:
qsub -cwd -j y -q all.q \
-o ./CTCF_12.txt \ # name of output file
-b y -N CTCF_12 \ # name of scheduled job
Rscript ./ChIPanalyser_general.R CTCF 12 IMR90 1 10 11 60 topRegions DNase # run only this last line if I want to see errors as they appear

# List of TFs for IMR90: 
ELK1, CTCF, USF2, SMC3, RFX5, RCOR1, RAD21, NFE2L2, MXI1, MAZ, FOS, CHD1, BHLHE40, MAFK, CEBPB

# List of QDAs: 
0.00  0.10  0.20  0.30  0.40  0.50  0.60  0.70  0.80  0.90  0.95  0.99
  1    2      3     4    5      6     7     8     9    10    11    12
  
# IMR90 - ATAC, MNase, DNase, NOMe -  top, middle and bottom regions - DONE