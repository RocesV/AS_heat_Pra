########################################################################################
################################# 5. Meta-Networks  ####################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(seqRFLP))
suppressPackageStartupMessages(library(missForest))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(neat))
suppressPackageStartupMessages(library(fmsb))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(gridExtra))

###### 2. Input Data & Create MetaNetworks Input ######
## Network analysis with all the EVENTS that passed both tests
## Side_note01: If only one of the isoforms from the event passed both tests i will consider whole events isoforms 

Isoforms <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.Intersection/Doublediff/3.Doublediff_1_resShrink_vstrlog_git.xlsx", sheet = "vst_DoubDiff_Passed")
KisSplice <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.Splicing_FinalAnnotation_k51_Complementary_git.xlsx")

# Impute values: Code from preprocess_omics_list (Processomics pkg) -- like MOFA2 input
colnames(vst_MetaNetworks)[2:7] <- c("C1", "C2", "T1.1", "T1.2", "T3.1", "T3.2")
vst_MetaNetworks_C <- vst_MetaNetworks[,c(2,3)]
vst_MetaNetworks_C$C3 <- NA
vst_MetaNetworks_T1 <- vst_MetaNetworks[,c(4,5)]
vst_MetaNetworks_T1$T1.3 <- NA
vst_MetaNetworks_T3 <- vst_MetaNetworks[,c(6,7)]
vst_MetaNetworks_T3$T3.3 <- NA

vst_MetaNetworks_imputed <- list(C = t(vst_MetaNetworks_C), T1 = t(vst_MetaNetworks_T1), T3 = t(vst_MetaNetworks_T3))

vst_MetaNetworks_imputed <- lapply(vst_MetaNetworks_imputed,function(matriz2){
  imputed <- matriz2
  registerDoParallel(cores= detectCores() - 1)
  imputed<-(missForestedit(imputed,parallelize="variable"))$ximp
  matriz2 <- imputed
  matriz2
})

vst_MetaNetworks_Seidr <- t(do.call(rbind, vst_MetaNetworks_imputed))
colnames(vst_MetaNetworks_Seidr) <-  c(rep("C",3), rep("T1",3), rep("T3",3))
rownames(vst_MetaNetworks_Seidr) <- vst_MetaNetworks$Network_names
write.table(vst_MetaNetworks_Seidr, file = "vst_Seidr.txt", sep = " ")

# Add some names
KissIDs <- read.delim("clipboard", header = T) # from 1.Seidr_inputs | AllKissIDs
KissIDs$Symbol <- NA
KissIDs$KissID <- as.character(KissIDs$KissID)
datas <- KissIDs
for(j in 1:nrow(datas)){
  if(!is.na(as.character(Splicing_annotation$DESC_quickGO[which(Splicing_annotation$KissID == datas$KissID[j])]))){
    datas$Symbol[j] <- as.character(Splicing_annotation$DESC_quickGO[which(Splicing_annotation$KissID == datas$KissID[j])])  
  }else if(is.na(as.character(Splicing_annotation$DESC_quickGO[which(Splicing_annotation$KissID == datas$KissID[j])]))){
    if(!is.na(as.character(Splicing_annotation$DESC_InterPro[which(Splicing_annotation$KissID == datas$KissID[j])]))){
      datas$Symbol[j] <- as.character(Splicing_annotation$DESC_InterPro[which(Splicing_annotation$KissID == datas$KissID[j])])
    }else if(is.na(as.character(Splicing_annotation$DESC_InterPro[which(Splicing_annotation$KissID == datas$KissID[j])]))){
      if(!is.na(as.character(Splicing_annotation$DESC_Panther[which(Splicing_annotation$KissID == datas$KissID[j])]))){
        datas$Symbol[j] <- as.character(Splicing_annotation$DESC_Panther[which(Splicing_annotation$KissID == datas$KissID[j])])
      }else if(is.na(as.character(Splicing_annotation$DESC_Panther[which(Splicing_annotation$KissID == datas$KissID[j])]))){
        datas$Symbol[j] <- datas$KissID[j]
      }
    } 
  }
}

datas$Symbol <- paste0(datas$Symbol, ".", rownames(datas))
write.table(datas, file = "AllKissIDs_names.txt", sep = "\t", col.names = T, row.names = F)

datas <- read.delim("clipboard", header = T, dec = ",", row.names = 1)
datas <- t(datas)

# isoforms that do not vary at all create problems, so it´s better to drop them

vars <- apply(datas, 2, var)
filt_id <- which(is.finite(vars))
datas <- datas[,filt_id] # all passed the filter as expected because they are doudblediff

# Center samples around their median, has been shown to improve reconstruction accuracy

medians <- apply(datas, 1, median)
datas <- sweep(datas, MARGIN = 1, FUN = '-', STATS = medians)

MASS::write.matrix(x = unname(datas), sep = "\t", file = "expression.tsv")
write(colnames(datas), file = "genes.txt")

###### 3. Seidr:Meta-Networks Toolkit  ######
## Run Seidr toolkit
## Windows10 Pro- wsl 2 - Ubuntu 20.04


###### 4. Network BackBonning and Infomap community discovering  ######
## Run Seidr toolkit and informap
## Windows10 Pro- wsl 2 - Ubuntu 20.04

###### 5. Network Enrichment Analysis ######
## Input data | from 3.Seidr_MetaNetworks_Infomap & 2.Heatmap
network_128 <- read.delim("clipboard", header = T)
network_1 <- read.delim("clipboard", header = T)
cluster_set_128 <- read.delim("clipboard", header = T)
cluster_set_128 <- cluster_set_128[,c(1,ncol(cluster_set_128))]
cluster_set_1 <- read.delim("clipboard", header = T)
cluster_set_1 <- cluster_set_1[,c(1,ncol(cluster_set_1))]
functional_set <- read.delim("clipboard", header = T)

## Analysis 1: NEAT Topology Enrichment - import final network from 4.Seidr_MetaNetwork_BackBone_nc
# NEAT answer the following question: if one cluster is related to a functional bin, the number of links
# should be larger o smaller than we would expect to observ by chance
# Side_note: bad idea idea if one of your sets contains all nodes of the network

# Define sets as lists: cluster_set (Infomap clusters) 
# Define sets as lists: functional_set (Mercator MapMan functional bins)

cluster_sets_dfs <- list(cluster_set_1 = cluster_set_1, cluster_set_128 = cluster_set_128)

cluster_sets_list <- lapply(cluster_sets_dfs, function(x){
  cluster_sets <- list()
  for(i in levels(as.factor(x$Infomap))){
    cluster_sets[[i]] <- x$KissID[which(x$Infomap == i)]
  }
  cluster_sets
})


functional_sets_list <- list()
for(i in levels(as.factor(functional_set$bincodes))){
  functional_sets_list[[i]] <- functional_set$KissID[which(functional_set$bincodes == i)]
}
desc <- read.delim("clipboard", header = T)
names(functional_sets_list) == rownames(desc)
names(functional_sets_list) <- desc$description

network_list <- list(network_1 = network_1, network_128 = network_128)

# Computation of the test


Network_Enrichment_1 <- neat(alist = cluster_sets_list[[1]], blist = functional_sets_list, network = as.matrix(network_list[[1]]), 
                             nettype = 'directed', nodes = cluster_sets_dfs[[1]]$KissID, alpha = 0.05)
Network_Enrichment_128 <- neat(alist = cluster_sets_list[[2]], blist = functional_sets_list, network = as.matrix(network_list[[2]]), 
                             nettype = 'directed', nodes = cluster_sets_dfs[[2]]$KissID, alpha = 0.05)
write.table(Network_Enrichment_1, file = "NEAT_1.txt", sep = "\t")
write.table(Network_Enrichment_128, file = "NEAT_128.txt", sep = "\t")

# plot the results: Radar plots
# all radars show bins for which some Over/Under enrichment is detected
# pdf with radars for each cluster || only for Network_Enrichment_1 (final network)
coul <- brewer.pal(3, "YlGnBu")
colors_border <- coul
colors_in <- alpha(coul,0.3)

Cluster.Enrichments <- list()

for(j in levels(as.factor(Network_Enrichment_1$A))){
  if(length(which(Network_Enrichment_1$A == j & Network_Enrichment_1$conclusion == "No enrichment")) == 29){
    cat(paste0("\n No enrichments detected for ", j, " Informap Cluster \n"))
    next
  }
  categories <- unique(Network_Enrichment_1$B[which(Network_Enrichment_1$conclusion != "No enrichment")])
  df <- Network_Enrichment_1[which(Network_Enrichment_1$A == j),]
  df <- df[which(df$B %in% categories),]
  df$value <- NA
  df$adjusted_p[which(df$conclusion == "No enrichment")] <- 1
  df$value <- -log10(df$adjusted_p)
  max <- rep(13 , nrow(df))
  min <- rep(0  , nrow(df))
  # prepare rows
  df <- df[order(df$B, decreasing = F),]
  up <- df$value
  down <- df$value
  up[which(df$conclusion == "Underenrichment")] <- 0
  down[which(df$conclusion == "Overenrichment")] <- 0
  # categories in same order
  Rad = rbind(max, min,up , down)
  colnames(Rad) <- df$B
  Cluster.Enrichments[[j]] <- radarchart(as.data.frame(Rad), axistype = 0, 
             # polygon
             pcol = colors_border, pfcol = colors_in, plwd = 4, plty = 1,
             # grid 
             cglcol = "grey", cglty = 1, axislabcol = "grey", cglwd = 0.8,
             # label
             vlcex = 1, title = paste0("Infomap Cluster:", j)
             )}

###### 6. Session info  ######
writeLines(capture.output(sessionInfo()), "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Logs_SessionInfo/MetaNetworks_SelectRegulators_sessionInfo.txt")
