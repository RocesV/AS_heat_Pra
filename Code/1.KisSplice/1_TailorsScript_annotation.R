########################################################################################
############################## Functional Annotation  ##################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(seqRFLP))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))

###### 2. Fasta construction -- retain KisSplice IDs ######

Trinity.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_Trinity_sorted.bam")
FConsensus.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_FinalConsensus_sorted.bam")
Transsbyss.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_Transabyss_sorted.bam")
NTranssbyss.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_TransabyssNOTNORM_sorted.bam")
Trinity.fas <- read.fasta("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Data/Other assemblies/Trinity.fasta")
FConsensus.fas <- read.fasta("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Data/Other assemblies/Pra_transcripts_db_200524.fa")
Transsbyss.fas <- read.fasta("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Data/Other assemblies/transabyss-merged.fa")
NTranssbyss.fas <- read.fasta("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Data/Other assemblies/transabyss-mergedNOTNORM.fa")
KisSplice <- readDNAStringSet("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.results_k51_type1.fa")
KisSplice.names <- vapply(names(KisSplice), FUN.VALUE = character(1), USE.NAMES = F,FUN = function(x){new <- paste0(strsplit(x, split = "|", fixed = T)[[1]][1], sep = "|", strsplit(x, split = "|", fixed = T)[[1]][2])})
KisSplice.names <- paste0(KisSplice.names, sep = ".", c(1,2))
names(KisSplice) <- KisSplice.names

Chimera_query <- list(Trinity = as.data.frame(matrix(NA, nrow = length(KisSplice.names), ncol = 2)), FConsesus = as.data.frame(matrix(NA, nrow = length(KisSplice.names), ncol = 2)), Transabyss = as.data.frame(matrix(NA, nrow = length(KisSplice.names), ncol = 2)), NTranabyss = as.data.frame(matrix(NA, nrow = length(KisSplice.names), ncol = 2)))
Chimera_query <- lapply(Chimera_query, FUN = function(x){
  colnames(x) <- c("IDs", "sequences")
  x[,1] <- KisSplice.names
  x
}) 


.fas <- list(Trinity = Trinity.fas, FConsensus = FConsensus.fas, Transabyss = Transsbyss.fas, NTranssbyss = NTranssbyss.fas)
.bam <- list(Trinity = Trinity.bam, FConsensus = FConsensus.bam, Transabyss = Transsbyss.bam, NTranssbyss = NTranssbyss.bam)

for(j in 1:length(.fas)){
  for(i in 1:length(KisSplice.names)){
    # dont care about dups because i dont have that much seqs and want final data.frame straigh foward
    condition <- which(.bam[[j]][[1]]$qname == KisSplice.names[i])
    if(length(condition) == 1){# detected or redundant
      Chimera_query[[j]][i,2] <- toupper(paste0(.fas[[j]][[which(names(.fas[[j]]) == as.character(.bam[[j]][[1]]$rname[condition]))]],
                                                collapse = ""))
    } else if(length(condition) < 1){# missing or orphan
      # care conditional preventing orphan isoforms
      if(substr(KisSplice.names[i], nchar(KisSplice.names[i]), nchar(KisSplice.names[i])) == "1"){
        NNames <- gsub(pattern = ".1", replacement = ".2", x = KisSplice.names[i], fixed = T)
      } else if(substr(KisSplice.names[i], nchar(KisSplice.names[i]), nchar(KisSplice.names[i])) == "2"){
        NNames <- gsub(pattern = ".2", replacement = ".1", x = KisSplice.names[i], fixed = T)
      }
      
      # kisSplice >> Transabyss so same KisSplice events it is suposed to be annotated as the same transcripts -- symmetrics 
      condition2 <- which(.bam[[j]][[1]]$qname == NNames)
      if(length(condition2) == 1){ # orphan
        Chimera_query[[j]][i,2] <- toupper(paste0(.fas[[j]][[which(names(.fas[[j]]) == as.character(.bam[[j]][[1]]$rname[condition2]))]],
                                                  collapse = ""))
      } else if(length(condition2) < 1){ # missing -- retain KisSplice seqs
        Chimera_query[[j]][i,2] <-  as.character(KisSplice[[i]])
        
      } else if(length(condition2) > 1){ stop("\n You are doing something wrong boy \n")}
    } else if(length(condition) > 1){ stop("\n You are doing something wrong boy \n")}
  }
}

# Export

dataframe2fas(Chimera_query$Transabyss, file = paste0("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Final query/", names(Chimera_query[3]), ".fa"))
dataframe2fas(Chimera_query$NoTransabyss, file = paste0("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Final query/", names(Chimera_query[4]), ".fa"))
dataframe2fas(Chimera_query$Trinity, file = paste0("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Final query/", names(Chimera_query[1]), ".fa"))
dataframe2fas(Chimera_query$FConsesus, file = paste0("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Final query/", names(Chimera_query[2]), ".fa"))


###### 3. First annotation round -- custom dammit + Mercator ######
# pass to Ubuntu and online plug

###### 4. Complementary identification ######
# I will do separately Mercator and dammit, finally will join them

#### Mercator ####

setwd("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Outputs/1_First_dammitMercator/mercator/")
path <- "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Outputs/1_First_dammitMercator/mercator/"
  
Complementary.splicing <- list()
for(i in list.files(".")){
  Mercator.fasta <- readDNAStringSet(paste0("./", i, "/", i, "_Fasta_mapman.txt"))
  Identifications <- names(Mercator.fasta)
  if(length(Identifications) != 13890){stop("\n You are doing something wrong boy \n")}
  final.table <- as.data.frame(matrix(NA, nrow = length(Identifications),  ncol = 5))
  colnames(final.table) <- c("KissID", "BinCo", "NameCo", "Bin", "Name")
  for(j in 1:length(Identifications)){
    final.table$KissID[j] <- strsplit(Identifications[j], split = " | ", fixed = T)[[1]][1]
    if(length(strsplit(Identifications[j], split = " | ", fixed = T)[[1]]) < 3){
      macro <- strsplit(Identifications[j], split = " | ", fixed = T)[[1]][2]  
    }else if(length(strsplit(Identifications[j], split = " | ", fixed = T)[[1]]) == 3){
      macro <- strsplit(Identifications[j], split = " | ", fixed = T)[[1]][3] 
     }else if(length(strsplit(Identifications[j], split = " | ", fixed = T)[[1]]) > 3){
      pseudo_macro <- strsplit(Identifications[j], split = " | ", fixed = T)[[1]]
      if(length(grep("RNA", pseudo_macro)) > 0){
        macro <- pseudo_macro[grep("RNA", pseudo_macro)[1]]
      }else if(length(grep("RNA", pseudo_macro)) == 0){
        macro <- pseudo_macro[2]}
    }else {stop("Check your conditionals boy")}
    if(macro == "mapman4 classification: 35.2 not assigned.not annotated"){
      final.table$BinCo[j] <- "35.2"
      final.table$NameCo[j] <- "not assigned.not annotated"
      final.table$Bin[j] <- "35"
      final.table$Name[j] <- NA
    }else if(macro != "mapman4 classification: 35.2 not assigned.not annotated"){
      if(str_count(macro, ":") < 4){
        # see strings || conditional between different string beacause some of them same number of : but no complete name so bad row
        if(strsplit(macro, split = " : ")[[1]][4] != " mapman4.2.0"){
          final.table$BinCo[j] <- strsplit(macro, split = ":", fixed = T)[[1]][2]
          final.table$NameCo[j] <- strsplit(macro, split = ":", fixed = T)[[1]][3]
          final.table$Name[j] <- strsplit(macro, split = ":", fixed = T)[[1]][1]
          final.table$Bin[j] <- strsplit(final.table$BinCo[j], split = ".", fixed = T)[[1]][1]  
        }else if(strsplit(macro, split = " : ")[[1]][4] == " mapman4.2.0"){
          final.table$BinCo[j] <- strsplit(macro, split = " : ")[[1]][3]
          final.table$Bin[j] <- strsplit(final.table$BinCo[j], split = ".", fixed = T)[[1]][1]
          final.table$NameCo[j] <- paste0(strsplit(macro, split = " : ")[[1]][1], " : ", strsplit(macro, split = " : ")[[1]][2])
          final.table$Name[j] <- final.table$NameCo[j]
        }
      }else if(str_count(macro, ":") >= 4){
        final.table$BinCo[j] <-  strsplit(macro, split = ":", fixed = T)[[1]][4]
        final.table$Bin[j] <- strsplit(final.table$BinCo[j], split = ".", fixed = T)[[1]][1]
        final.table$NameCo[j] <- paste0(strsplit(macro, split = ":", fixed = T)[[1]][1], " : ", strsplit(macro, split = ":", fixed = T)[[1]][2], " : ", strsplit(macro, split = ":", fixed = T)[[1]][3])
        final.table$Name[j] <- strsplit(macro, split = ":", fixed = T)[[1]][5]
      }
    }
  }
  Complementary.splicing[[i]] <- final.table
  rm(final.table)
  rm(Mercator.fasta)
  rm(i)
  rm(j)
  rm(macro)
  rm(Identifications)
  rm(pseudo_macro)
}

### Complementary identification between same splicing event ###

Complemetary.assembler <- lapply(Complementary.splicing, FUN = function(x){
  Identifications <- substr(x[,1], 1, nchar(x[,1])-2)
  for(i in 1:nrow(x)){
    ide <- x[i,1]
    event.ide <- gsub(".1", "", ide, fixed = T)
    event.ide <- gsub(".2", "", event.ide, fixed = T)
    event.number <- which(Identifications == event.ide)
    if(x[event.number[1],4] != x[event.number[2],4]){    # Bin
     if(x[event.number[1],4] == "35" | x[event.number[2],4] == "35"){
       x[event.number[which(x[event.number,4] == "35")],2] <- x[event.number[which(x[event.number,4] != "35")],2]
       x[event.number[which(x[event.number,4] == "35")],3] <- x[event.number[which(x[event.number,4] != "35")],3]
       x[event.number[which(x[event.number,4] == "35")],5] <- x[event.number[which(x[event.number,4] != "35")],5]
       x[event.number[which(x[event.number,4] == "35")],4] <- x[event.number[which(x[event.number,4] != "35")],4]
       }
    }else if(x[event.number[1],4] == x[event.number[2],4]){    # Name
      if(is.na(x[event.number[1],5]) & !is.na(x[event.number[2],5]) | is.na(x[event.number[2],5]) & !is.na(x[event.number[1],5])){
        x[event.number[is.na(x[event.number,5])],3] <- x[event.number[!is.na(x[event.number,5])],3]
        x[event.number[is.na(x[event.number,5])],5] <- x[event.number[!is.na(x[event.number,5])],5]
      }  
    }
  }
  x
})

### Complementary identification between assemblers ###

colnames(Complemetary.assembler$KisSplice_FConsensus) <- paste0(colnames(Complemetary.assembler$KisSplice_FConsensus), "_FConsensus")
colnames(Complemetary.assembler$KisSplice_NoTransabyss) <- paste0(colnames(Complemetary.assembler$KisSplice_NoTransabyss), "_NoTransabyss")
colnames(Complemetary.assembler$KisSplice_Transabyss) <- paste0(colnames(Complemetary.assembler$KisSplice_Transabyss), "_Transabyss")
colnames(Complemetary.assembler$KisSplice_Trinity) <- paste0(colnames(Complemetary.assembler$KisSplice_Trinity), "_Trinity")
Complementary.final <- do.call(cbind, Complemetary.assembler)
Complementary.final.bin <- Complementary.final[,c(1,2,4,7,9,12,14,17,19)]
Complementary.final.name <- Complementary.final[,c(1,3,5,8,10,13,15,18,20)]

## Priority stablished by bowtie2 mapping results -- NoTransabyss >> Transabyss >> Trinity >> FConsensus || Other criteria may be most common result but it may be 35 so discarded

# Bin: manual revision in each line

Complementary.fina.bin.max <- Complementary.final.bin[,c(1,4,5)] #NoTransabyss columns
colnames(Complementary.fina.bin.max) <- c("KissID", "BinCo", "Bin")
compare.before <- length(which(Complementary.fina.bin.max$Bin == "35"))

Complementary.final.bin$KisSplice_Transabyss.BinCo_Transabyss[c(2939,2940, 8575, 8576, 11003, 11004)] <- c(" 5.2.4.3 ", " 5.2.4.3 ", rep(" 50.1.14 ", 4))
Complementary.final.bin$KisSplice_Transabyss.Bin_Transabyss[c(2939,2940, 8575, 8576, 11003, 11004)] <- c("  5", " 5", rep(" 50", 4))
Complementary.fina.bin.max$BinCo[which(Complementary.final.bin$KisSplice_NoTransabyss.Bin_NoTransabyss == "35" & Complementary.final.bin$KisSplice_Transabyss.Bin_Transabyss != "35")] <- Complementary.final.bin$KisSplice_Transabyss.BinCo_Transabyss[which(Complementary.final.bin$KisSplice_NoTransabyss.Bin_NoTransabyss == "35" & Complementary.final.bin$KisSplice_Transabyss.Bin_Transabyss != "35")]
Complementary.fina.bin.max$Bin[which(Complementary.final.bin$KisSplice_NoTransabyss.Bin_NoTransabyss == "35" & Complementary.final.bin$KisSplice_Transabyss.Bin_Transabyss != "35")] <- Complementary.final.bin$KisSplice_Transabyss.Bin_Transabyss[which(Complementary.final.bin$KisSplice_NoTransabyss.Bin_NoTransabyss == "35" & Complementary.final.bin$KisSplice_Transabyss.Bin_Transabyss != "35")]

Complementary.final.bin$KisSplice_Trinity.BinCo_Trinity[c(1903,1904,2681,2682,7213,7214,10135,10136)] <- c(rep(" 50.1.14 ", 6), rep(" 24.2.8.1", 2))
Complementary.final.bin$KisSplice_Trinity.Bin_Trinity[c(1903,1904,2681,2682,7213,7214,10135,10136)] <- c(rep(" 50 ", 6), rep(" 24", 2))
Complementary.fina.bin.max$BinCo[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_Trinity.Bin_Trinity != "35")] <- Complementary.final.bin$KisSplice_Trinity.BinCo_Trinity[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_Trinity.Bin_Trinity != "35")]
Complementary.fina.bin.max$Bin[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_Trinity.Bin_Trinity != "35")] <- Complementary.final.bin$KisSplice_Trinity.Bin_Trinity[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_Trinity.Bin_Trinity != "35")]

Complementary.fina.bin.max$BinCo[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_FConsensus.Bin_FConsensus != "35")] <- Complementary.final.bin$KisSplice_FConsensus.BinCo_FConsensus[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_FConsensus.Bin_FConsensus != "35")]
Complementary.fina.bin.max$Bin[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_FConsensus.Bin_FConsensus != "35")] <- Complementary.final.bin$KisSplice_FConsensus.Bin_FConsensus[which(Complementary.fina.bin.max$Bin == "35" & Complementary.final.bin$KisSplice_FConsensus.Bin_FConsensus != "35")]
compare.after <- length(which(Complementary.fina.bin.max$Bin == "35"))

compare.before - compare.after # plus classification: 600 | Sequences Classified: 6018
write.table(Complementary.fina.bin.max, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Outputs/2_Complementary identification/mercator/bin/Complementary_final_bins.txt", sep = "\t")

# Names: manual revision in each line

Complementary.fina.name.max <- Complementary.final.name[,c(1,4,5)] #NoTransabyss columns
colnames(Complementary.fina.name.max) <- c("KissID", "NameCo", "Name")
Compare.before <- length(which(is.na(Complementary.fina.name.max$Name) == TRUE))

Complementary.fina.name.max$NameCo[is.na(Complementary.final.name$KisSplice_NoTransabyss.Name_NoTransabyss) & !is.na(Complementary.final.name$KisSplice_Transabyss.Name_Transabyss)] <- Complementary.final.name$KisSplice_Transabyss.NameCo_Transabyss[is.na(Complementary.final.name$KisSplice_NoTransabyss.Name_NoTransabyss) & !is.na(Complementary.final.name$KisSplice_Transabyss.Name_Transabyss)]
Complementary.fina.name.max$Name[is.na(Complementary.final.name$KisSplice_NoTransabyss.Name_NoTransabyss) & !is.na(Complementary.final.name$KisSplice_Transabyss.Name_Transabyss)] <- Complementary.final.name$KisSplice_Transabyss.Name_Transabyss[is.na(Complementary.final.name$KisSplice_NoTransabyss.Name_NoTransabyss) & !is.na(Complementary.final.name$KisSplice_Transabyss.Name_Transabyss)]

Complementary.fina.name.max$NameCo[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_Trinity.Name_Trinity)] <- Complementary.final.name$KisSplice_Trinity.NameCo_Trinity[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_Trinity.Name_Trinity)] 
Complementary.fina.name.max$Name[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_Trinity.Name_Trinity)] <- Complementary.final.name$KisSplice_Trinity.Name_Trinity[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_Trinity.Name_Trinity)]

Complementary.fina.name.max$NameCo[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_FConsensus.Name_FConsensus)] <- Complementary.final.name$KisSplice_FConsensus.NameCo_FConsensus[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_FConsensus.Name_FConsensus)]
Complementary.fina.name.max$Name[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_FConsensus.Name_FConsensus)] <- Complementary.final.name$KisSplice_FConsensus.Name_FConsensus[is.na(Complementary.fina.name.max$Name) & !is.na(Complementary.final.name$KisSplice_FConsensus.Name_FConsensus)]
Compare.after <- length(which(is.na(Complementary.fina.name.max$Name) == TRUE))

Compare.before - Compare.after # plus annotation: 700 | Sequences Annotated: 8270
write.table(Complementary.fina.name.max, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Outputs/2_Complementary identification/mercator/names/Complementary_final_names.txt", sep = "\t")

# Final Excel: manual revision in Excel

Complementary.max <- cbind(Complementary.fina.bin.max, Complementary.fina.name.max[,-1])
write.table(Complementary.max, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Outputs/2_Complementary identification/mercator/final/Complementary_final_woRevision.txt", sep = "\t")

#### Dammit ####

setwd("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Outputs/1_First_dammitMercator/dammit/")
path <- "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional Annotation/Outputs/1_First_dammitMercator/dammit/"

### Custom dammit.parser ###

dammit.parsero.NoTransabyss <- gotubro.parser(dammit.output = "./KisSplice_NoTransabyss.fa.dammit/KisSplice_NoTransabyss.fa.dammit.gff3",
                                              dammit.namemap = "./KisSplice_NoTransabyss.fa.dammit/KisSplice_NoTransabyss.fa.dammit.namemap.csv")

dammit.parsero.Transabyss <- gotubro.parser(dammit.output = "./KisSplice_Transabyss.fa.dammit/KisSplice_Transabyss.fa.dammit.gff3",
                                              dammit.namemap = "./KisSplice_Transabyss.fa.dammit/KisSplice_Transabyss.fa.dammit.namemap.csv")


dammit.parsero.Trinity <- gotubro.parser(dammit.output = "./KisSplice_Trinity.fa.dammit/KisSplice_Trinity.fa.dammit.gff3",
                                            dammit.namemap = "./KisSplice_Trinity.fa.dammit/KisSplice_Trinity.fa.dammit.namemap.csv")

dammit.parsero.FConsensus <- gotubro.parser(dammit.output = "./KisSplice_FConsesus.fa.dammit/KisSplice_FConsesus.fa.dammit.gff3",
                                            dammit.namemap = "./KisSplice_FConsesus.fa.dammit/KisSplice_FConsesus.fa.dammit.namemap.csv")

Complementary.splicing.dammit <- list(NoTransabyss = dammit.parsero.NoTransabyss, Transabyss = dammit.parsero.Transabyss,
                                      Trinity = dammit.parsero.Trinity, FConsensus = dammit.parsero.FConsensus)

### Complementary identification between same splicing event ###

Complemetary.assembler.dammit <- lapply(Complementary.splicing.dammit, FUN = function(x){
  x[,1] <- as.character(x[,1])
  Identifications <- substr(x[,1], 1, nchar(x[,1])-2)
  for(i in 1:nrow(x)){
    ide <- x[i,1]
    event.ide <- gsub(".1", "", ide, fixed = T)
    event.ide <- gsub(".2", "", event.ide, fixed = T)
    event.number <- which(Identifications == event.ide)
    for(j in colnames(x)[-c(1:2)]){
      if(is.na(x[event.number[1],j]) & !is.na(x[event.number[2],j]) | !is.na(x[event.number[1],j]) & is.na(x[event.number[2],j])){ # complementary splicing event annotation
        x[event.number[is.na(x[event.number,j])],j] <- x[event.number[!is.na(x[event.number,j])],j]
      }else if(is.na(x[event.number[1],j]) & is.na(x[event.number[2],j]) | !is.na(x[event.number[1],j]) & !is.na(x[event.number[2],j])){ # no annotatios needed | next iter 
        next 
        }
    }
  }
  x
})

### Complementary identification between assemblers ###

Complementary.final.dammit <- Complemetary.assembler.dammit$NoTransabyss # bowtie2 stablished priority | as Mercator above

for(i in 2:length(Complemetary.assembler.dammit)){
  for(j in colnames(Complementary.final.dammit)[-c(1,2)]){
    if(length(which(is.na(Complementary.final.dammit[,j]) & !is.na(Complemetary.assembler.dammit[[i]][,j]) == TRUE)) == 0){
      cat(paste0("\n Complementary.final.dammit ", j, " is more complete than ", names(Complemetary.assembler.dammit[i]), " ", j, " \n"))
    }else if(length(which(is.na(Complementary.final.dammit[,j]) & !is.na(Complemetary.assembler.dammit[[i]][,j]) == TRUE)) > 0){
      Complementary.final.dammit[is.na(Complementary.final.dammit[,j]) & !is.na(Complemetary.assembler.dammit[[i]][,j]),j] <- Complemetary.assembler.dammit[[i]][is.na(Complementary.final.dammit[,j]) & !is.na(Complemetary.assembler.dammit[[i]][,j]),j]  
    }
  }
}

write.table(Complementary.final.dammit, file = "../../2_Complementary identification/dammit/Complementary_final_dammit.txt", sep = "\t", col.names = T)

#### InterPro Scan ####

setwd("../Outputs/3_Functional_annotation/FunctionalCrossFiles_FCF/InterproScan_PfamRfam")

### Formatting ###

for(i in list.files()){
  file <- read.delim(i, sep = "\t", header = F)
  file$V1 <- as.character(file$V1)
  new.df <- as.data.frame(matrix(NA, nrow = nrow(file), ncol = 3))
  colnames(new.df) <- c("ID", "desc", "GO")
  for(j in 1:nrow(file)){
    if(i == "HAMAP2go"){
      new.df$ID[j] <- gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1])
      new.df$desc[j] <- strsplit(strsplit(file$V1[j], split = " > ", fixed = T)[[1]][2], split = " ; ", fixed = T)[[1]][1]
      new.df$GO[j] <- strsplit(strsplit(file$V1[j], split = " > ", fixed = T)[[1]][2], split = " ; ", fixed = T)[[1]][2]
    }else if(i == "InterPro2go"){ 
      new.df$ID[j] <- substr(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1]), 1, 9)
      new.df$desc[j] <- substr(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1]), 10, nchar(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1])))
      new.df$GO[j] <- strsplit(strsplit(file$V1[j], split = " > ", fixed = T)[[1]][2], split = " ; ", fixed = T)[[1]][2]
    }else if(i == "PIRSF2go"){
      new.df$ID[j] <- substr(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1]), 1, 11)
      new.df$desc[j] <- substr(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1]), 12, nchar(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1])))
      new.df$GO[j] <- strsplit(strsplit(file$V1[j], split = " > ", fixed = T)[[1]][2], split = " ; ", fixed = T)[[1]][2]
    }else if(i == "PRINTS2go" | i == "PROSITE2go" | i == "SMART2go" | i == "Pfam2go" | i == "Rfam2go"){
      new.df$ID[j] <- substr(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1]), 1, 7)
      new.df$desc[j] <- substr(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1]), 8, nchar(gsub(pattern = paste0(substr(i, 1, nchar(i) - 3), ":") , replacement = "", strsplit(file$V1[j], split = " > ", fixed = T)[[1]][1])))
      new.df$GO[j] <- strsplit(strsplit(file$V1[j], split = " > ", fixed = T)[[1]][2], split = " ; ", fixed = T)[[1]][2]
    }else {stop(cat("\n Something go wrong boy \n"))}
  }
  write.table(new.df, file = paste0("./", i, "_formated.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
}

setwd("../Outputs/3_Functional_annotation/FunctionalCrossFiles_FCF/Panther")

PANTHER <- read.delim("panther2go", header = F, sep = "\t")
panther.formatted <- as.data.frame(matrix(NA, nrow(PANTHER), 5))
colnames(panther.formatted) <- c("ID", "desc", "GO_MF", "GO_BP", "GO_CC")
panther.formatted$ID <- PANTHER$V1
panther.formatted$desc <- PANTHER$V2
panther.formatted$GO_MF <- as.character(PANTHER$V3)
panther.formatted$GO_BP <- as.character(PANTHER$V4)
panther.formatted$GO_CC <- as.character(PANTHER$V5)
for(i in 1:nrow(panther.formatted)){
  if(!panther.formatted$GO_MF[i] == ""){
    panther.formatted$GO_MF[i] <-  paste0(do.call(rbind,strsplit(strsplit(panther.formatted$GO_MF[i], split = ";", fixed = T)[[1]], split = "#", fixed = T))[,2], collapse = ";")
  }
  if(!panther.formatted$GO_BP[i] == ""){
    panther.formatted$GO_BP[i] <-  paste0(do.call(rbind,strsplit(strsplit(panther.formatted$GO_BP[i], split = ";", fixed = T)[[1]], split = "#", fixed = T))[,2], collapse = ";")
  }
  if(!panther.formatted$GO_CC[i] == ""){
    panther.formatted$GO_CC[i] <-  paste0(do.call(rbind,strsplit(strsplit(panther.formatted$GO_CC[i], split = ";", fixed = T)[[1]], split = "#", fixed = T))[,2], collapse = ";")
  }
}
write.table(panther.formatted, file = "./panther2go_formated.txt", quote = F, col.names = T, row.names = F, sep = "\t")


### Custom InterPro.parser ###

interpro.parsero.NoTransabyss <- gotubro.parser(interpro.output = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/interpro/KisSplice_NoTransabyss/",
                                              dammit.namemap = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/dammit/KisSplice_NoTransabyss.fa.dammit/KisSplice_NoTransabyss.fa.dammit.namemap.csv",
                                              interpro = TRUE)

interpro.parsero.Transabyss <- gotubro.parser(interpro.output = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/interpro/KisSplice_Transabyss/",
                                            dammit.namemap = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/dammit/KisSplice_Transabyss.fa.dammit/KisSplice_Transabyss.fa.dammit.namemap.csv",
                                            interpro = TRUE)


interpro.parsero.Trinity <- gotubro.parser(interpro.output = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/interpro/KisSplice_Trinity/",
                                         dammit.namemap = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/dammit/KisSplice_Trinity.fa.dammit/KisSplice_Trinity.fa.dammit.namemap.csv",
                                         interpro = TRUE)

interpro.parsero.FConsensus <- gotubro.parser(interpro.output = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/interpro/KisSplice_FConsensus/",
                                            dammit.namemap = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/1_First_dammitMercator_interpro/dammit/KisSplice_FConsesus.fa.dammit/KisSplice_FConsesus.fa.dammit.namemap.csv",
                                            interpro = TRUE)

Complementary.splicing.interpro <- list(NoTransabyss = interpro.parsero.NoTransabyss, Transabyss = interpro.parsero.Transabyss,
                                      Trinity = interpro.parsero.Trinity, FConsensus = interpro.parsero.FConsensus)

### Complementary identification between same splicing event ###

Complemetary.assembler.interpro <- lapply(Complementary.splicing.interpro, FUN = function(x){
  x[,1] <- as.character(x[,1])
  Identifications <- substr(x[,1], 1, nchar(x[,1])-2)
  for(i in 1:nrow(x)){
    ide <- x[i,1]
    event.ide <- gsub(".1", "", ide, fixed = T)
    event.ide <- gsub(".2", "", event.ide, fixed = T)
    event.number <- which(Identifications == event.ide)
    for(j in colnames(x)[-1]){
      if(is.na(x[event.number[1],j]) & !is.na(x[event.number[2],j]) | !is.na(x[event.number[1],j]) & is.na(x[event.number[2],j])){ # complementary splicing event annotation
        x[event.number[is.na(x[event.number,j])],j] <- x[event.number[!is.na(x[event.number,j])],j]
      }else if(is.na(x[event.number[1],j]) & is.na(x[event.number[2],j]) | !is.na(x[event.number[1],j]) & !is.na(x[event.number[2],j])){ # no annotatios needed | next iter 
        next 
      }
    }
  }
  x
})

### Complementary identification between assemblers ###

Complementary.final.interpro <- Complemetary.assembler.interpro$NoTransabyss # bowtie2 stablished priority | as Mercator and dammit above

for(i in 2:length(Complemetary.assembler.interpro)){
  for(j in colnames(Complementary.final.interpro)[-1]){
    if(length(which(is.na(Complementary.final.interpro[,j]) & !is.na(Complemetary.assembler.interpro[[i]][,j]) == TRUE)) == 0){
      cat(paste0("\n Complementary.final.interpro ", j, " is more complete than ", names(Complemetary.assembler.interpro[i]), " ", j, " \n"))
    }else if(length(which(is.na(Complementary.final.interpro[,j]) & !is.na(Complemetary.assembler.interpro[[i]][,j]) == TRUE)) > 0){
      Complementary.final.interpro[is.na(Complementary.final.interpro[,j]) & !is.na(Complemetary.assembler.interpro[[i]][,j]),j] <- Complemetary.assembler.interpro[[i]][is.na(Complementary.final.interpro[,j]) & !is.na(Complemetary.assembler.interpro[[i]][,j]),j]  
    }
  }
}

write.table(Complementary.final.interpro, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/2_Complementary identification/interpro/Complementary_final_interpro.txt", sep = "\t", col.names = T)

###### 5. Pass to functional annotation (gene onthology) ######

Splicing.GOs <- gotubro.crossrefs(gotubro.parser.dammit = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/2_Complementary identification/dammit/Complementary_final_dammit.txt",
                                  gotubro.parser.ip =  "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/2_Complementary identification/interpro/Complementary_final_interpro.txt",
                                  SQLite3db = "F:/SQLite3/gotubro", 
                                  quickGO = TRUE, 
                                  mode = c("reduce", "complete", "normal"), 
                                  return.merged = T)

write.table(Splicing.GOs$reduce, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/3_Functional_annotation/Outputs/GOs_gotubro_reduce.txt", sep = "\t")
write.table(Splicing.GOs$complete, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/3_Functional_annotation/Outputs/GOs_gotubro_complete.txt", sep = "\t")
write.table(Splicing.GOs$normal, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Functional_Annotation/Outputs/3_Functional_annotation/Outputs/GOs_gotubro_normal.txt", sep = "\t")