########################################################################################
################################### 6. Validation  #####################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(seqRFLP))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(d3r))
suppressPackageStartupMessages(library(sunburstR))

###### 2. Data prepare ######
KisSplice <- read.fasta("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.results_k51_type1.fa")
KisSplice.names <- KisSplice[which((1:length(KisSplice) %% 2) == 1)]
KisSplice.names <- vapply(KisSplice.names, FUN.VALUE = character(1), USE.NAMES = F,FUN = function(x){new <- paste0(strsplit(x, split = "|", fixed = T)[[1]][1], sep = "|", strsplit(x, split = "|", fixed = T)[[1]][2])})
KisSplice.names <- paste0(KisSplice.names, sep = ".", c(1,2))
KisSplice.events <- length(KisSplice.names)/2 # 2 isoforms per event

KisSplice.ok <- seqRFLP::rename.fas(KisSplice, names = KisSplice.names)
write.fasta(sequences = KisSplice.ok, file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Data/Other assemblies/KisSplice_ok.fasta")

###### 3. Mapping to different assemblies ######
# Run Script 6_Mapping_Kissplice.sh

###### 4. Splicing event compara ######
## Import Bams and KisSplice events
Trinity.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_Trinity_sorted.bam")
Transsbyss.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_Transabyss_sorted.bam")
Old.Consensus.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_Consensus_sorted.bam")
Final.Consensus.bam <- scanBam(file = "D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Validation/Figures/In silico/KisSplice_FinalConsensus_sorted.bam")
Events <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.k51_Annotation_SYMMETRIC.xlsx")
Events <- levels(as.factor(Events$events.names))

## Event Compara. KisSplice Events are missed, detected, orphan or redundant by other assemblers?
# Missing - no mapped; detected - same events isoforms hit different contigs
# Redundant - same event isoforms hit same contig; orphan - only one isoform (2 per event) mapped
All.assembly <- list(Transabyss = Transsbyss.bam, Trinity = Trinity.bam, Old.Consensus = Old.Consensus.bam, Final.Consensus = Final.Consensus.bam)

All.assembly.2 <- lapply(All.assembly, FUN = function(x){
  x[[1]]$qname <- gsub(".1", "", x[[1]]$qname, fixed = T)
  x[[1]]$qname <- gsub(".2", "", x[[1]]$qname, fixed = T)
  missing <- list()
  detected <- list()
  redundant <- list()
  orphan <- list()
  for(i in Events){
    hit <- x[[1]]$rname[which(x[[1]]$qname == i)]
    if(length(hit) == 0){
      missing[[i]] <- i
    }else if(length(hit) == 1){
      orphan[[i]] <- i
    }else if(length(hit) == 2){
      if(hit[1] == hit[2]){ redundant[[i]] <- i}
      if(hit[1] != hit[2]){ detected[[i]]  <- i}
      # with this is enough because it was performed end to end mapping
    }else if(length(hit) > 2){stop(cat("\n Something go wrong .. \n"))}
  }
  bar <- data.frame(Condition = c("Missing", "Orphan", "Redundant", "Detected"), Value = c(length(missing), length(orphan), length(redundant), length(detected)),
                    Percentage =  c(length(missing)/KisSplice.events, length(orphan)/KisSplice.events, length(redundant)/KisSplice.events, length(detected)/KisSplice.events))
  Event.Compara <- list(missing = unlist(missing, use.names = F), orphan = unlist(orphan, use.names = F), redundant = unlist(redundant, use.names = F), detected = unlist(detected, use.names = F))
  
  x[[1]]$bar <- bar
  x[[1]]$Event.Compara <- Event.Compara
  x
})

###### 5. Plot the results ######
### Stacked bar plot

data <- rbind(All.assembly.2$Transabyss[[1]]$bar, All.assembly.2$Trinity[[1]]$bar, All.assembly.2$Old.Consensus[[1]]$bar, All.assembly.2$Final.Consensus[[1]]$bar)
data$Assembly <- c(rep("Transabyss", 4), rep("Trinity", 4), rep("Old.Consensus", 4), rep("Fin.Consensus", 4))

ggplot(data, aes(fill=Condition, y=Percentage, x=Assembly)) + 
  geom_bar(position="fill", stat="identity") + theme_light() + scale_fill_manual(values = brewer.pal(4, "Pastel1")) +
  ggtitle("Event compara") + xlab(element_blank()) + ylab(element_blank())

### Sunburst plot
# Prepare data
size2 <- list()
names <- list()
All.events <- list()
x <- 0
for(i in 1:length(All.assembly.2$Transabyss[[1]]$Event.Compara)){
  for(j in 1:length(All.assembly.2$Transabyss[[1]]$Event.Compara)){
    for(z in 1:length(All.assembly.2$Transabyss[[1]]$Event.Compara)){
      x <- x +1
      Trans.Trin <- which((All.assembly.2$Transabyss[[1]]$Event.Compara[[j]] %in% All.assembly.2$Trinity[[1]]$Event.Compara[[i]]) == TRUE)
      Trans.Trin.2 <- All.assembly.2$Transabyss[[1]]$Event.Compara[[j]][Trans.Trin]
      Cons.Trans.Trin <- length(which((All.assembly.2$Final.Consensus[[1]]$Event.Compara[[z]] %in% Trans.Trin.2)))
      size2[x] <- Cons.Trans.Trin
      names[[x]] <- c(paste("Trinity", names(All.assembly.2$Transabyss[[1]]$Event.Compara[i]), sep = "."), paste("Transabbys", names(All.assembly.2$Transabyss[[1]]$Event.Compara[j]), sep = "."), paste("F.Consensus", names(All.assembly.2$Transabyss[[1]]$Event.Compara[z]), sep = "."))
      jeje <- paste(names[[x]], collapse = "_")
      Cons.Trans.Trin.2 <- which((All.assembly.2$Final.Consensus[[1]]$Event.Compara[[z]] %in% Trans.Trin.2))
      All.events[[jeje]] <- All.assembly.2$Final.Consensus[[1]]$Event.Compara[[z]][Cons.Trans.Trin.2]
  }
 }
}

dat <- data.frame(
  level1 = rep("KisSplice.Events",64),
  level2 = NA,
  level3 = NA,
  level4 = NA,
  size = unlist(size2, use.names = F)
)

for(i in 1:nrow(dat)){
  dat$level2[i] <- names[[i]][1]
  dat$level3[i] <- names[[i]][2]
  dat$level4[i] <- names[[i]][3]
}

tree <- d3_nest(dat, value_cols = "size")
tree

sb1 <- sunburst(tree, width="100%", height=600,
                colors = c("#74ADD1", # all.detected consesus detected
                           "darksalmon", # inside all circle
                           "white", #trinity.missing
                           "#FEE08B" , # trinity.orphan
                           "#FDAE61", # trinity redundant 
                           "#F46D43", # trinity.detected 
                           "white", # transabyss.missing 
                           "#D9EF8B", # transabyss.orphan 
                           "#A6D96A", # transabyss.redundant
                           "#66BD63", # transabyss.detected
                           "white", # consensus.missing
                           "#E0F3F8", # consensus.orphan
                           "#ABD9E9"),
                legend = list(w = 155)) # 

###### 6. Export some results ######
All.events <- lapply(All.events, function(x){
  paste(x, collapse = ",")
})
export <- data.frame(Compara = names(All.events), Events = unlist(All.events, use.names = F))

write.table(dat, file = "./Validation/Data/v2Event_namesEvent_numberEvent_Compara/Event_number.txt", sep = " ")
write.table(export, file = "./Validation/Data/v2Event_namesEvent_numberEvent_Compara/Event_names.txt", sep = " ")

###### 7. Session info  ######
writeLines(capture.output(sessionInfo()), "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Logs_SessionInfo/ValidationInsilico_sessionInfo.txt")
