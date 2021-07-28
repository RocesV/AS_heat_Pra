########################################################################################
############################ 1. KISSPLICE ##############################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######
suppressPackageStartupMessages(library(seqRFLP))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(stringr))

###### 2. Define Functions ######

compare.DNA <- function(x,y){
  ########################
  # compare.DNA: search differences between two sequences
  ########################
  # x - sequence of nts
  # y - sequence of nts
  ########################
  sapply(seq(length(x)),
         function(i){
           x[i] == y[i]
         }
  )
}

Variants.seq.diff <- function(matrix1, kmer, export = TRUE, path) {
  ########################
  # Variants.seq.diff: KisSplice extension. Finds the variable sequence (S) between two isoforms. 
  ########################
  # matrix1 - matrix or dataframe (n events, ), first column kissIDs, second column Sequences and third column lengths.
  # kmer - kmer used in KisSplice
  # export - logical. Export S variable seq as fasta to indicated path
  ########################
  Upper <- matrix1[which((as.numeric(rownames(KisSplice)) %% 2) == 0),]
  Lower <- matrix1[which((as.numeric(rownames(KisSplice)) %% 2) == 1),]
  matrix1 <- data.frame(kissID = Upper[,1], Upper.Seq = Upper[,2], Lower.seq = Lower[,2], Length = Upper[,3] - Lower[,3])
  lower <- as.matrix(matrix1[,3])
  upper <- as.matrix(matrix1[,2])
  list.lower <- list()
  list.upper <- list()
  list.length <- list()
  for (i in seq_along(lower)) {
    seq_lower <- unlist(strsplit(lower[i,], split = ""))
    seq_upper <- unlist(strsplit(upper[i,], split = ""))
    list.lower[[i]] <- seq_lower
    list.upper[[i]] <- seq_upper
    list.length[[i]] <- length(seq_upper) - length(seq_lower) 
  }
  
  list.Svariable <- list()
  for (k in 1:length(list.lower)) {
    TrueFalse <- compare.DNA(list.upper[[k]], list.lower[[k]])
    S0 <- grep("FALSE", TrueFalse)[1]
    if(S0 != (kmer + 1)){print(stop("\n Error01: wrong kmer"))}
    Svariable <- list.upper[[k]][S0: (S0+ list.length[[k]] -1)]
    list.Svariable[[k]] <- paste(Svariable,collapse = "")
  }
  df <- data.frame(KissID = matrix1[,1],
                   Lower_seq = lower[,1], 
                   Upper_seq = upper[,1],
                   Length_Variable = unlist(list.length, use.names = FALSE),
                   Variable_seq = unlist(list.Svariable, use.names = FALSE))
  if(export){
    dataframe2fas(df[,c(1,5)], file = paste0(path, "allSparts.fasta"))
  }
  return(df)
}


IR.enrichment <- function(Variants.seq.diff.object, exitron = TRUE){
  ########################
  # IR.enrichment: KisSplice extension. Finds potential retained introns (trying to find exitrons too) in S variables sequences. Based in most frequent intron flanking sites of Arabidopsis
  ########################
  # Variants.seq.diff.object - dataframe (n events, ) containing lower isoforms, upper isoform and S variable sequences
  # exitron - logical. Search potential exitrons based in features described in Staiger and SimpsonGenome Biology (2015) 16:136
  ########################
  returned <- list()
  Filter <- str_sub(as.character(Variants.seq.diff.object[,5]), start = 1, end = 2) == "GT" & str_sub(as.character(Variants.seq.diff.object[,5]), start = -2) == "AG" & nchar(as.character(Variants.seq.diff.object[,5])) > 50
  returned[['introns']] <-  Variants.seq.diff.object[Filter,]
  returned[['Non_introns']] <- Variants.seq.diff.object[!Filter,]
  if(exitron == FALSE) {return(returned)}
  if(exitron == TRUE){
  Filter2 <- str_sub(as.character(Variants.seq.diff.object[,5]), start = 1, end = 2) == "GT" | str_sub(as.character(Variants.seq.diff.object[,5]), start = 1, end = 2) == "GC" | str_sub(as.character(Variants.seq.diff.object[,5]), start = 1, end = 2) == "AT"
  All.pot <- Variants.seq.diff.object[Filter2,]
  Filter3 <- str_sub(as.character(All.pot[,5]), start = -2) == "AG" | str_sub(as.character(All.pot[,5]), start = -2) == "TG" | str_sub(as.character(All.pot[,5]), start = -2) == "AA"
  All.pot2 <- All.pot[Filter3,]
  Filter4 <- nchar(as.character(All.pot2[,5])) > 50
  All.pot <- All.pot2[Filter4,]
  EIx3 <- nchar(as.character(All.pot[,5])) %%3 == 0
  onlysequence <- All.pot[EIx3,5]
  onlysequencedf <- as.data.frame(onlysequence)
  onlysequencedf$onlysequence <- as.character(onlysequencedf$onlysequence)
  StopCodons <- c("TAG", "TAA", "TGA")
  Finally <- as.data.frame(matrix(ncol = 4, nrow = length(which(EIx3 == TRUE))))
  Finally[,1] <- All.pot[EIx3,1]
  colnames(Finally) <- c("ID", "StopCodonds", "%GC", "Sequence")
  Finally[,4] <- onlysequence
    for (i in 1:nrow(onlysequencedf)) {
      bananasplit <- strsplit(onlysequencedf[i,1] , "")[[1]]
      GCporcentaje <- ((sum(str_count(bananasplit, pattern = "G")) + sum(str_count(bananasplit, pattern = "C"))) / length(bananasplit)) *100
      Finally[i,3] <- GCporcentaje
      bananasplitxxi <- paste0(bananasplit[c(TRUE, FALSE, FALSE)], bananasplit[c(FALSE, TRUE, FALSE)], bananasplit[c(FALSE, FALSE, TRUE)])
      Number <- unique(grep( paste(StopCodons, collapse = "|"), bananasplitxxi))
      if (length(Number) != 0){Filling <- 1} else {Filling <- 0}
      Finally[i,2] <- Filling
    } 
  returned[['Exitrons']] <- Finally[which(Finally[,2]  == 0),]
  returned[['introns']][,"%GC"] <- NA
    for (l in 1:nrow(returned[['introns']])) {
      intronesTruesplit <- strsplit(as.character(returned[['introns']][l,2]) , "")[[1]]
      GCporcentaje <- ((sum(str_count(intronesTruesplit, pattern = "G")) + sum(str_count(intronesTruesplit, pattern = "C"))) / length(intronesTruesplit)) *100
      returned[['introns']][l,"%GC"] <- GCporcentaje
    }
  
  cat(paste0("\n Potential Retained Introns Events:", sep = " ", round(nrow(returned$introns)/nrow(Variants.seq.diff.object) * 100, digits = 3), "%"))
  return(returned)
  } }


###### 3. Obtain Results ######
# Import data from KisSplice table. odd rows = lower isoforms. even rows = upper isoforms.

KisSplice <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/3.Variants.diff.xlsx")


# All S variable sequences

allSparts <- Variants.seq.diff(KisSplice, kmer = 51, export = TRUE, path = "D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/") 

# Potential retained introns and exitrons

AS.enrichment <- IR.enrichment(allSparts, exitron = TRUE)

###### 4. Session info ######

writeLines(capture.output(sessionInfo()), "D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Logs_SessionInfo/VariantsSeq_IREnrichment_sessionInfo.txt")
