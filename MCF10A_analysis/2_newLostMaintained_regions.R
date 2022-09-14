###### Build GRanges from excel of peaks files to build datasets of new, lost and maintained regions within early and late regions


cell_line <- "MCF10"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/peaks/"))
library(GenomicRanges)

# Ctrl atac
df_1 <- read.csv("MCF10A_CTRL_24H_rep1_peaks.xls", sep ="\t", header=T)
gr_1 <- makeGRangesFromDataFrame(df_1,keep.extra.columns=T) # One grange object per replicate
df_2 <- read.csv("MCF10A_CTRL_24H_rep2_peaks.xls", sep ="\t", header=T)
gr_2 <- makeGRangesFromDataFrame(df_2,keep.extra.columns=T)
atac_ctrl_24 <- intersect(gr_1,gr_2) # Intersect the 2 replicates

df_1 <- read.csv("MCF10A_CTRL_48H_rep1_peaks.xls", sep ="\t", header=T)
gr_1 <- makeGRangesFromDataFrame(df_1,keep.extra.columns=T)
df_2 <- read.csv("MCF10A_CTRL_48H_rep2_peaks.xls", sep ="\t", header=T)
gr_2 <- makeGRangesFromDataFrame(df_2,keep.extra.columns=T)
atac_ctrl_48 <- intersect(gr_1,gr_2)

ctrl_atac <- union(atac_ctrl_24, atac_ctrl_48)

# HER2 atac
df_1 <- read.csv("MCF10A_HER2_24H_rep1_peaks.xls", sep ="\t", header=T)
gr_1 <- makeGRangesFromDataFrame(df_1,keep.extra.columns=T)
df_2 <- read.csv("MCF10A_HER2_24H_rep1_peaks.xls", sep ="\t", header=T)
gr_2 <- makeGRangesFromDataFrame(df_2,keep.extra.columns=T)
atac_her2_24 <- intersect(gr_1,gr_2)

df_1 <- read.csv("MCF10A_HER2_48H_rep1_peaks.xls", sep ="\t", header=T)
gr_1 <- makeGRangesFromDataFrame(df_1,keep.extra.columns=T)
df_2 <- read.csv("MCF10A_HER2_48H_rep1_peaks.xls", sep ="\t", header=T)
gr_2 <- makeGRangesFromDataFrame(df_2,keep.extra.columns=T)
atac_her2_48 <- intersect(gr_1,gr_2)

her2_atac <- union(atac_her2_24, atac_her2_48)



# Find maintained, lost and new regions
lost_strong <- ctrl_atac[!overlapsAny(ctrl_atac, her2_atac, maxgap=1000)]
lost_reduced <- reduce(lost_strong, min.gapwidth=1000)
save(lost_reduced, file = "lost_regions.Rda")
lost_reduced <-  lost_reduced[seqnames(lost_reduced)%in%seqlevels(lost_reduced)[1:24]]
lost_reduced <-  lost_reduced + 2000
seqlevels(lost_reduced) <- seqlevelsInUse(lost_reduced)
names(lost_reduced) <- paste0(seqnames(lost_reduced),":",start(lost_reduced),"..",end(lost_reduced))
save(lost_reduced, file = "lost_regions_withNames.Rda")



new_strong <- her2_atac[!overlapsAny(her2_atac, ctrl_atac, maxgap=1000)]
new_reduced <- reduce(new_strong, min.gapwidth=1000)
save(new_reduced, file = "new_regions.Rda")
new_reduced <-  new_reduced[seqnames(new_reduced)%in%seqlevels(new_reduced)[1:24]]
new_reduced <-  new_reduced + 2000
seqlevels(new_reduced) <- seqlevelsInUse(new_reduced)
names(new_reduced) <- paste0(seqnames(new_reduced),":",start(new_reduced),"..",end(new_reduced))
save(new_reduced, file = "new_regions_withNames.Rda")



maintained <- intersect(ctrl_atac,her2_atac)
maintained_reduced <- reduce(maintained, min.gapwidth=1000) 
save(maintained_reduced, file = "maintained_regions.Rda")
maintained_reduced <-  maintained_reduced[seqnames(maintained_reduced)%in%seqlevels(maintained_reduced)[1:24]]
maintained_reduced <-  maintained_reduced + 2000
seqlevels(maintained_reduced) <- seqlevelsInUse(maintained_reduced)
names(maintained_reduced) <- paste0(seqnames(maintained_reduced),":",start(maintained_reduced),"..",end(maintained_reduced))
save(maintained_reduced, file = "maintained_regions_withNames.Rda")



# Plot histograms
pdf(file="New_regions_histogram.pdf")
hist(log10(width(new_reduced)), main=NULL)
title(main ="Histogram of log10 of gained regions between late and early HER2 overexpression analysis", cex.main = 0.8)
dev.off() # Plot histograms

pdf(file="Lost_regions_histogram.pdf")
hist(log10(width(lost_reduced)), main=NULL)
title(main="Histogram of log10 of lost regions between early and late HER2 overexpression analysis", cex.main = 0.8)
dev.off()

pdf(file="Maintained_regions_histogram.pdf")
hist(log10(width(maintained_reduced)), main=NULL)
title(main="Histogram of log10 of maintained regions between early and late HER2 overexpression analysis", cex.main = 0.8)
dev.off()