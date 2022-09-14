################################################################################
# HER2 PLOTS
################################################################################

# The example script below is just for one TF (ATF1), over all regions (lost, new and maintained).
# But the analysis below was actually run over all TFs.
# Below it's just an example of how I obtained the files used to make the plots.

library(ChIPanalyser)
binsPerSide <- 20
setwd("/storage/st20d/HiC_Kc167/alessandra/MCF10/rdaObjects/")



################################################################################
# ATF1
################################################################################
TF <- "ATF1"



# LOST REGIONS
regions <- "lost"

# Read data lost borders - control
if(file.exists("forMila/chip_ctrl_lostRegions_ATF1_matrix.RData")){
  load("forMila/chip_ctrl_lostRegions_ATF1_matrix.RData")
} else{
  #chip_ctrl_lostRegions_ATF1 <- get(load(paste0(TF,"/",regions,"/predicted_chip_ctrl_lostRegions_MCF10_ATF1_gw.Rda")))
  chip_ctrl_lostRegions_ATF1 <- get(load(paste0("C:/Users/mitica l'ale/Desktop/cluster/predicted_chip_ctrl_lostRegions_MCF10_ATF1_gw.Rda")))
  chip_ctrl_lostRegions_ATF1_matrix <- matrix(0, ncol=binsPerSide*2+1, nrow=length(chip_ctrl_lostRegions_ATF1@profiles[[1]]))
  for(i in 1:length(chip_ctrl_lostRegions_ATF1@profiles[[1]])){
    points_available <- length(chip_ctrl_lostRegions_ATF1@profiles[[1]][[i]]$ChIP)
    points_keep <- binsPerSide*2+1
    points_start <- round((points_available-points_keep)/2)+1
    points_stop <- points_start+points_keep-1
    chip_ctrl_lostRegions_ATF1_matrix[i,] <- chip_ctrl_lostRegions_ATF1@profiles[[1]][[i]]$ChIP[points_start:points_stop]
  }
  save(chip_ctrl_lostRegions_ATF1_matrix, file="forMila/chip_ctrl_lostRegions_ATF1_matrix.RData")
}

# Read data lost borders - her2 over expression
if(file.exists("forMila/chip_her2_lostRegions_ATF1_matrix.RData")){
  load("forMila/chip_her2_lostRegions_ATF1_matrix.RData")
} else{
  #chip_her2_lostRegions_ATF1 <- get(load(paste0(TF,"/",regions,"/predicted_chip_her2_lostRegions_MCF10_ATF1_gw.Rda")))
  chip_her2_lostRegions_ATF1 <- get(load(paste0("C:/Users/mitica l'ale/Desktop/cluster/predicted_chip_her2_lostRegions_MCF10_ATF1_gw.Rda")))
  chip_her2_lostRegions_ATF1_matrix <- matrix(0, ncol=binsPerSide*2+1, nrow=length(chip_her2_lostRegions_ATF1@profiles[[1]]))
  for(i in 1:length(chip_her2_lostRegions_ATF1@profiles[[1]])){
    points_available <- length(chip_her2_lostRegions_ATF1@profiles[[1]][[i]]$ChIP)
    points_keep <- binsPerSide*2+1
    points_start <- round((points_available-points_keep)/2)+1
    points_stop <- points_start+points_keep-1
    chip_her2_lostRegions_ATF1_matrix[i,] <- chip_her2_lostRegions_ATF1@profiles[[1]][[i]]$ChIP[points_start:points_stop]
  }
  save(chip_her2_lostRegions_ATF1_matrix, file="forMila/chip_her2_lostRegions_ATF1_matrix.RData")
}

# Compute difference lost borders between her2 over expression and control
chip_diff_lostRegions_ATF1_matrix <- log2(chip_her2_lostRegions_ATF1_matrix/chip_ctrl_lostRegions_ATF1_matrix)
save(chip_diff_lostRegions_ATF1_matrix, file="forMila/chip_diff_lostRegions_ATF1_matrix.RData")



############################################################################################################



# NEW REGIONS
regions <- "new"

# Read data new borders - control
if(file.exists("forMila/chip_ctrl_newRegions_ATF1_matrix.RData")){
  load("forMila/chip_ctrl_newRegions_ATF1_matrix.RData")
} else{
  #chip_ctrl_newRegions_ATF1 <- get(load(paste0(TF,"/",regions,"/predicted_chip_ctrl_newRegions_MCF10_ATF1_gw.Rda")))
  chip_ctrl_newRegions_ATF1 <- get(load(paste0("C:/Users/mitica l'ale/Desktop/cluster/predicted_chip_ctrl_newRegions_MCF10_ATF1_gw.Rda")))
  chip_ctrl_newRegions_ATF1_matrix <- matrix(0, ncol=binsPerSide*2+1, nrow=length(chip_ctrl_newRegions_ATF1@profiles[[1]]))
  for(i in 1:length(chip_ctrl_newRegions_ATF1@profiles[[1]])){
    points_available <- length(chip_ctrl_newRegions_ATF1@profiles[[1]][[i]]$ChIP)
    points_keep <- binsPerSide*2+1
    points_start <- round((points_available-points_keep)/2)+1
    points_stop <- points_start+points_keep-1
    chip_ctrl_newRegions_ATF1_matrix[i,] <- chip_ctrl_newRegions_ATF1@profiles[[1]][[i]]$ChIP[points_start:points_stop]
  }
  save(chip_ctrl_newRegions_ATF1_matrix, file="forMila/chip_ctrl_newRegions_ATF1_matrix.RData")
}

# Read data new borders - her2 over expression
if(file.exists("forMila/chip_her2_newRegions_ATF1_matrix.RData")){
  load("forMila/chip_her2_newRegions_ATF1_matrix.RData")
} else{
  #chip_her2_newRegions_ATF1 <- get(load(paste0(TF,"/",regions,"predicted_chip_her2_newRegions_MCF10_ATF1_gw.Rda")))
  chip_her2_newRegions_ATF1 <- get(load(paste0("C:/Users/mitica l'ale/Desktop/cluster/predicted_chip_her2_newRegions_MCF10_ATF1_gw.Rda")))
  chip_her2_newRegions_ATF1_matrix <- matrix(0, ncol=binsPerSide*2+1, nrow=length(chip_her2_newRegions_ATF1@profiles[[1]]))
  for(i in 1:length(chip_her2_newRegions_ATF1@profiles[[1]])){
    points_available <- length(chip_her2_newRegions_ATF1@profiles[[1]][[i]]$ChIP)
    points_keep <- binsPerSide*2+1
    points_start <- round((points_available-points_keep)/2)+1
    points_stop <- points_start+points_keep-1
    chip_her2_newRegions_ATF1_matrix[i,] <- chip_her2_newRegions_ATF1@profiles[[1]][[i]]$ChIP[points_start:points_stop]
  }
  save(chip_her2_newRegions_ATF1_matrix, file="forMila/chip_her2_newRegions_ATF1_matrix.RData")
}

# Compute difference new borders between her2 over expression and control
chip_diff_newRegions_ATF1_matrix <- log2(chip_her2_newRegions_ATF1_matrix/chip_ctrl_newRegions_ATF1_matrix)
save(chip_diff_newRegions_ATF1_matrix, file="forMila/chip_diff_newRegions_ATF1_matrix.RData")



############################################################################################################



# MAINTAINED REGIONS
regions <- "maintained"

# Read data maintained borders - control
if(file.exists("forMila/chip_ctrl_maintainedRegions_ATF1_matrix.RData")){
  load("forMila/chip_ctrl_maintainedRegions_ATF1_matrix.RData")
} else{
  #chip_ctrl_maintainedRegions_ATF1 <- get(load(paste0(TF,"/",regions,"/predicted_chip_ctrl_maintainedRegions_MCF10_ATF1_gw.Rda")))
  chip_ctrl_maintainedRegions_ATF1 <- get(load(paste0("C:/Users/mitica l'ale/Desktop/cluster/predicted_chip_ctrl_maintainedRegions_MCF10_ATF1_gw.Rda")))
  chip_ctrl_maintainedRegions_ATF1_matrix <- matrix(0, ncol=binsPerSide*2+1, nrow=length(chip_ctrl_maintainedRegions_ATF1@profiles[[1]]))
  for(i in 1:length(chip_ctrl_maintainedRegions_ATF1@profiles[[1]])){
    points_available <- length(chip_ctrl_maintainedRegions_ATF1@profiles[[1]][[i]]$ChIP)
    points_keep <- binsPerSide*2+1
    points_start <- round((points_available-points_keep)/2)+1
    points_stop <- points_start+points_keep-1
    chip_ctrl_maintainedRegions_ATF1_matrix[i,] <- chip_ctrl_maintainedRegions_ATF1@profiles[[1]][[i]]$ChIP[points_start:points_stop]
  }
  save(chip_ctrl_maintainedRegions_ATF1_matrix, file="forMila/chip_ctrl_maintainedRegions_ATF1_matrix.RData")
}

# Read data maintained borders - her2 over expression
if(file.exists("forMila/chip_her2_maintainedRegions_ATF1_matrix.RData")){
  load("forMila/chip_her2_maintainedRegions_ATF1_matrix.RData")
} else{
  #chip_her2_maintainedRegions_ATF1 <- get(load(paste0(TF,"/",regions,"/predicted_chip_her2_maintainedRegions_MCF10_ATF1_gw.Rda")))
  chip_her2_maintainedRegions_ATF1 <- get(load(paste0("C:/Users/mitica l'ale/Desktop/cluster/predicted_chip_her2_maintainedRegions_MCF10_ATF1_gw.Rda")))
  chip_her2_maintainedRegions_ATF1_matrix <- matrix(0, ncol=binsPerSide*2+1, nrow=length(chip_her2_maintainedRegions_ATF1@profiles[[1]]))
  for(i in 1:length(chip_her2_maintainedRegions_ATF1@profiles[[1]])){
    points_available <- length(chip_her2_maintainedRegions_ATF1@profiles[[1]][[i]]$ChIP)
    points_keep <- binsPerSide*2+1
    points_start <- round((points_available-points_keep)/2)+1
    points_stop <- points_start+points_keep-1
    chip_her2_maintainedRegions_ATF1_matrix[i,] <- chip_her2_maintainedRegions_ATF1@profiles[[1]][[i]]$ChIP[points_start:points_stop]
  }
  save(chip_her2_maintainedRegions_ATF1_matrix, file="forMila/chip_her2_maintainedRegions_ATF1_matrix.RData")
}

# Compute difference maintained borders between her2 over expression and control
chip_diff_maintainedRegions_ATF1_matrix <- log2(chip_her2_maintainedRegions_ATF1_matrix/chip_ctrl_maintainedRegions_ATF1_matrix)
save(chip_diff_maintainedRegions_ATF1_matrix, file="forMila/chip_diff_maintainedRegions_ATF1_matrix.RData")



############################################################################################################


################################################################################
# Now do the same for all the other TFs.
################################################################################