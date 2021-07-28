########################################################################################
############################ 3. DESeq2  ################################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(vsn))

###### 2. Input Data ######
## Odd - Upper isoform; Even - Lower isoform
countData <- read.table("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.Raw_counts.txt", sep = " ")
rownames(countData) <- paste0(countData$events.names, sep = ".", c(1,2))
countData <- countData[,-c(1,2)]
colnames(countData) <- paste0(c(rep("C1",2), rep("C2",2), rep("T1_1",2), rep("T1_2",2), rep("T3_1",2), rep("T3_2",2)), sep = ".", c("1", "2"))

## Condition vector
# C = Control; T1 = Short-term stress (1-Day HS); T3 = Long/medium-term stress (3-Day HS)
condition <- factor(c(rep("C",4), rep("T1",4), rep("T3",4)))
# No remove: Vector represents all replicates -- ONLY FOR DATA INSPECTION
replicates1 <- factor(c(1:12))
# Removed Technical Replicates: Vector represents biological replicates -- USED DURING ANALYSIS
replicates2 <- factor(c(rep("1",2), rep("2",2), rep("3",2), rep("4",2), rep("5",2), rep("6",2)))
# Removed All Replicates: Vector represents treatments -- ONLY FOR DATA INSPECTION
replicates3 <- factor(c(rep("1",4), rep("2",4), rep("3",4)))
# Create DataFrame with condition info
colData <-  data.frame(condition = condition, All = replicates1, Biological = replicates2, Treatments = replicates3, row.names = colnames(countData))

## Check countData colnames match condition DataFrame rownames: Warning - DESeq2 will not internally check this
if(all(rownames(colData) %in% colnames(countData))){cat("\n Pass \n")} else{stop("\n Error01: check names pls \n")}

## Construct DESeqDataSet
# Construct
dds <- DESeqDataSetFromMatrix(countData, colData, ~condition)

# Add additional feature data, i.e: Isoform annotation or names (SMA3_Annotation_Symmetric) & functional ontology (Mercator_bin_symmetric). In order to avoid problems -- make unique identifiers
KisSplice <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.k51_Annotation_SYMMETRIC.xlsx")
KisSplice$events.names <- paste0(KisSplice$events.names, sep = ".", c(1,2))
gene <- list()
bin <- list()
for(i in 1:nrow(countData)){
  gene[[i]] <- paste0(KisSplice[which(KisSplice$events.names == rownames(countData)[i]),4], sep = "..", i)
  bin[[i]] <- paste0(KisSplice[which(KisSplice$events.names == rownames(countData)[i]),12], sep = "..", i)
}
AddfeatureData <- data.frame(gene = unlist(gene,use.names = F), bin = unlist(bin,use.names = F)) 
mcols(dds) <- DataFrame(mcols(dds), AddfeatureData)

###### 3. Differential Expression Analysis ######
exported <- list()
## Pre-filtering (not needed because it is automatically applied via independent filtering on the mean of normalized counts within results function)
keep <- rowSums(DESeq2::counts(dds)) >= 10
dds <- dds[keep,]

## Define the reference level to compare: not needed in this case but is a good practice 
dds$condition <- relevel(dds$condition, ref = "C")

##Collapse replicates: in analysis only technical replicates should be collapsed (Second). All replicates possibilites are collapsed for data inspection
dds_Bio <- collapseReplicates(dds, groupby = dds$Biological)

## Check that total counts are maintained post collapse
CHECK <- rowSums( DESeq2::counts(dds)[ , dds$Biological == "1" ])
if(all( CHECK == DESeq2::counts(dds_Bio)[ , "1" ])){cat("\n Pass \n")} else{ stop("\n Error02: Check collapse \n")}
rm(CHECK)

## Differential Expression
dds_Bio_DE <- DESeq(dds_Bio)

# Contrasts
exported[['resCT3']] <- results(dds_Bio_DE, contrast = c("condition", "T3", "C"))
exported[['resCT1']] <- results(dds_Bio_DE, contrast = c("condition", "T1", "C"))
exported[['resT1T3']] <- results(dds_Bio_DE, contrast = c("condition", "T3", "T1"))

# Log fold change shrinkage for visualization and ranking
exported[['resShrinkCT3']] <- lfcShrink(dds_Bio_DE, coef = 3, res = exported[['resCT3']], type = "apeglm")
exported[['resShrinkCT1']] <- lfcShrink(dds_Bio_DE, coef = 2, res = exported[['resCT1']], type = "apeglm")
exported[['resShrinkT1T3']] <- lfcShrink(dds_Bio_DE, contrast = c("condition", "T3", "T1"), res = exported[['resT1T3']], type = "ashr")

# p-values and adj.p-values. order, metrics (DE and strict DE (up-down regulation 3.24 fold)), summary and filter 
exported[['resCT3']] <- exported$resCT3[order(exported$resCT3$padj),]
exported[['resCT1']] <- exported$resCT1[order(exported$resCT1$padj),]
exported[['resT1T3']] <- exported$resT1T3[order(exported$resT1T3$padj),]

exported[['resShrinkCT3']] <- exported$resShrinkCT3[order(exported$resShrinkCT3$padj),]
exported[['resShrinkCT1']] <- exported$resShrinkCT1[order(exported$resShrinkCT1$padj),]
exported[['resShrinkT1T3']] <- exported$resShrinkT1T3[order(exported$resShrinkT1T3$padj),]

exported[['metrics']] <- data.frame(DE_CT1 = sum(exported$resCT1$padj < 0.05, na.rm = TRUE), 
Shrink_DE_CT1 = sum(exported$resShrinkCT1$padj < 0.05, na.rm = TRUE),
Strict_DE_CT1 = sum(exported$resCT1$padj < 0.05 & abs(exported$resCT1$log2FoldChange) > 1.8, na.rm = TRUE),
Strict_Shrink_DE_CT1 = sum(exported$resShrinkCT1$padj < 0.05 & abs(exported$resShrinkCT1$log2FoldChange) > 1.8, na.rm = TRUE),
DE_CT3 = sum(exported$resCT3$padj < 0.05, na.rm = TRUE),
Shrink_DE_CT3 = sum(exported$resShrinkCT3$padj < 0.05, na.rm = TRUE),
Strict_DE_CT3 = sum(exported$resCT3$padj < 0.05 & abs(exported$resCT3$log2FoldChange) > 1.8, na.rm = TRUE),
Strict_Shrink_DE_CT3 = sum(exported$resShrinkCT3$padj < 0.05 & abs(exported$resShrinkCT3$log2FoldChange) > 1.8, na.rm = TRUE),
DE_T1T3 = sum(exported$resT1T3$padj < 0.05, na.rm = TRUE),
Shrink_DE_T1T3 = sum(exported$resShrinkT1T3$padj < 0.05, na.rm = TRUE),
Strict_DE_T1T3 = sum(exported$resT1T3$padj < 0.05 & abs(exported$resT1T3$log2FoldChange) > 1.8, na.rm = TRUE),
Strict_Shrink_DE_T1T3 = sum(exported$resShrinkT1T3$padj < 0.05 & abs(exported$resShrinkT1T3$log2FoldChange) > 1.8, na.rm = TRUE))

exported[['summaryCT3']] <- capture.output(summary(exported$resCT3)) 
exported[['summaryCT1']] <- capture.output(summary(exported$resCT1))
exported[['summaryT1T3']] <- capture.output(summary(exported$resT1T3))

exported[['resCT3_DE']] <- exported$resCT3[which(exported$resCT3$padj < 0.05),]
exported[['resCT1_DE']] <- exported$resCT1[which(exported$resCT1$padj < 0.05),]
exported[['resT1T3_DE']] <- exported$resT1T3[which(exported$resT1T3$padj < 0.05),] 

exported[['resShrinkCT3_DE']] <- exported$resShrinkCT3[which(exported$resShrinkCT3$padj < 0.05),]
exported[['resShrinkCT1_DE']] <- exported$resShrinkCT1[which(exported$resShrinkCT1$padj < 0.05),]
exported[['resShrinkT1T3_DE']] <- exported$resShrinkT1T3[which(exported$resShrinkT1T3$padj < 0.05),]

exported[['resCT3_DE_Strict']] <- exported$resCT3[which(exported$resCT3$padj < 0.05 & abs(exported$resCT3$log2FoldChange) > 1.8),]
exported[['resCT1_DE_Strict']] <- exported$resCT1[which(exported$resCT1$padj < 0.05 & abs(exported$resCT1$log2FoldChange) > 1.8),]
exported[['resT1T3_DE_Strict']] <- exported$resT1T3[which(exported$resT1T3$padj < 0.05 & abs(exported$resT1T3$log2FoldChange) > 1.8),]

exported[['resShrinkCT3_DE_Strict']] <- exported$resShrinkCT3[which(exported$resShrinkCT3$padj < 0.05 & abs(exported$resShrinkCT3$log2FoldChange) > 1.8),]
exported[['resShrinkCT1_DE_Strict']] <- exported$resShrinkCT1[which(exported$resShrinkCT1$padj < 0.05 & abs(exported$resShrinkCT1$log2FoldChange) > 1.8),]
exported[['resShrinkT1T3_DE_Strict']] <- exported$resShrinkT1T3[which(exported$resShrinkT1T3$padj < 0.05 & abs(exported$resShrinkT1T3$log2FoldChange) > 1.8),]

# Diagnostic Dispersion/MA plots
DESeq2::plotDispEsts(dds_Bio_DE)
p1 <- DESeq2::plotMA(exported$resCT3)
p2 <- DESeq2::plotMA(exported$resShrinkCT3)
p3 <- DESeq2::plotMA(exported$resCT1)
p4 <- DESeq2::plotMA(exported$resShrinkCT1)
p5 <- DESeq2::plotMA(exported$resCT1)
p6 <- DESeq2::plotMA(exported$resShrinkCT1)
par(mfrow=c(3,2))
dev.off()

# Data transformation (variance stabilization; remove mean-variance dependence): rlog vs vst. Check which one is more stable. blind = FALSE for downstream analysis (Networks, Integration etc)
# vst does not assume any prior
DESeq2::sizeFactors(dds_Bio_DE) # not much variation in sizeFactors ("rlog tend has better performance than vst when exists 10x fold variation between the smallest and biggest sizeFactor" - M.Love (Developer))
vsd <- DESeq2::vst(dds_Bio_DE, blind = FALSE)
rld <- rlog(dds_Bio_DE, blind = FALSE)
meanSdPlot(assay(vsd)) 
meanSdPlot(assay(rld))
# vsd looks slightly better than rlog: i will use vst for downstream analysis
vsd <- as.data.frame(assay(vsd))
rld <- as.data.frame(assay(rld))
colnames(vsd) <- c(rep("C", 2), rep("T1", 2), rep("T3",2))
colnames(rld) <- c(rep("C", 2), rep("T1", 2), rep("T3",2))
meta <- mcols(dds)
if(all(rownames(meta) == rownames(vsd))){cat("\n Pass Annotation is straight-foward applied \n")} else{stop("\n Error03: Check annotations \n")}
if(all(rownames(meta) == rownames(rld))){cat("\n Pass Annotation is straight-foward applied \n")} else{stop("\n Error03: Check annotations \n")}
vsd <- cbind(gene = meta$gene,bin = meta$bin, kissID = rownames(vsd), vsd)
rld <- cbind(gene = meta$gene,bin = meta$bin, kissID = rownames(rld), rld)

exported[['vst']] <- vsd
exported[['rlog']] <- rld
# vst and rlog filtered by DE, DE_Strict and DE_Shrink_Strict
exported[['vst_DE_CT3']] <- exported$vst[rownames(exported$resCT3_DE),]
exported[['vst_DE_Strict_CT3']] <- exported$vst[rownames(exported$resCT3_DE_Strict),]
exported[['vst_DE_Shrink_Strict_CT3']] <- exported$vst[rownames(exported$resShrinkCT3_DE_Strict),]

exported[['vst_DE_CT1']] <- exported$vst[rownames(exported$resCT1_DE),]
exported[['vst_DE_Strict_CT1']] <- exported$vst[rownames(exported$resCT1_DE_Strict),]
exported[['vst_DE_Shrink_Strict_CT1']] <- exported$vst[rownames(exported$resShrinkCT1_DE_Strict),]

exported[['vst_DE_T1T3']] <- exported$vst[rownames(exported$resT1T3_DE),]
exported[['vst_DE_Strict_T1T3']] <- exported$vst[rownames(exported$resT1T3_DE_Strict),]
exported[['vst_DE_Shrink_Strict_T1T3']] <- exported$vst[rownames(exported$resShrinkT1T3_DE_Strict),]

DE_Passed <- unique(c(rownames(exported$resCT1_DE), rownames(exported$resT1T3_DE), rownames(exported$resCT3_DE)))
DE_Strict_Passed <- unique(c(rownames(exported$resCT1_DE_Strict), rownames(exported$resT1T3_DE_Strict), rownames(exported$resCT3_DE_Strict)))
DE_Shrink_Strict_Passed <- unique(c(rownames(exported$resShrinkCT1_DE_Strict), rownames(exported$resShrinkT1T3_DE_Strict), rownames(exported$resShrinkCT3_DE_Strict)))

exported[['vst_DE_final']] <- exported$vst[DE_Passed,]
exported[['vst_DE_Strict_final']] <- exported$vst[DE_Strict_Passed,]
exported[['vst_DE_Shrink_Strict_final']] <- exported$vst[DE_Shrink_Strict_Passed,]

###### 4. Integration with KissDE results ######
## Match Events (considering all events that passed any KissDE test and filtering low counts)
KissDE_CvsHS <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/1.Annotated_kissDE51_KmeansInput.xlsx", sheet = "kissDE_k51_C_HeatStress")
KissDE_CvsT1 <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/1.Annotated_kissDE51_KmeansInput.xlsx", sheet = "kissDE_k51_C_T1")
KissDE_CvsT3 <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/1.Annotated_kissDE51_KmeansInput.xlsx", sheet = "kissDE_k51_C_T3")
KissDE_T1vsT3 <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/1.Annotated_kissDE51_KmeansInput.xlsx", sheet = "kissDE_k51_T1_T3")

# KissDE events without low counts flag from KissDE
KissDE_Pass <- c(KissDE_T1vsT3$`#1.ID`[!as.logical(KissDE_T1vsT3$`21.lowcounts`)],
                 KissDE_CvsT3$`#1.ID`[!as.logical(KissDE_CvsT3$`21.lowcounts`)],
                 KissDE_CvsT1$`#1.ID`[!as.logical(KissDE_CvsT1$`21.lowcounts`)],
                 KissDE_CvsHS$`#1.ID`[!as.logical(KissDE_CvsHS$`29.lowcounts`)])

KissDE_Pass <- unique(KissDE_Pass)
KissDE_Pass_2 <- c(paste0(KissDE_Pass, sep = ".","1"), paste0(KissDE_Pass, sep = ".","2"))
# Retain only KissDE events that were not discarde in pre-filtering lowcounts DESeq2
KissDE_Pass_3 <- KissDE_Pass_2[which(KissDE_Pass_2 %in% rownames(exported$resCT3))]

## Filter Res, ResShrink, rlog and vsd (only KissDE and Double Dif)
# Only KissDE
exported[['resCT1_KissDE']] <- exported$resCT1[KissDE_Pass_3,]
exported[['resCT3_KissDE']] <- exported$resCT3[KissDE_Pass_3,]
exported[['resT1T3_KissDE']] <- exported$resT1T3[KissDE_Pass_3,]

exported[['resShrinkCT1_KissDE']] <- exported$resShrinkCT1[KissDE_Pass_3,]
exported[['resShrinkCT3_KissDE']] <- exported$resShrinkCT3[KissDE_Pass_3,]
exported[['resShrinkT1T3_KissDE']] <- exported$resShrinkT1T3[KissDE_Pass_3,]

exported[['rlog_KissDE']] <- exported$rlog[KissDE_Pass_3,]
exported[['vst_KissDE']] <- exported$vst[KissDE_Pass_3,]

# Double differential (KissDE-DESeq2)
# Side_Note Explanation 1: Unit: Each EVENT with two isoforms (lower and upper). Interpretation -- Events that passed KissDE test have variants/isoforms with diff expression during treatment. It compares expression between isoforms from the same event during treatment.
# Side_Note Explanation 2: Unit: Each ISOFORM. Interpretation -- Isoforms that passed DESeq2 test have differential expression during treatment. It compares unique isoform expression between two different conditions of the treatment.
# Side_Note Explanation 3: Unit: Integration DESeq2-KissDE. Unique isoforms that passed DESeq2 test may be not pass KissDE test because both isoforms from the same event have same expression pattern during treatment.
# Side_Note Explanation 4: Unit: Integration KissDE_DESeq2. Events (2 isoforms) that passed KissDE test may be not pass DESeq2 test because: 1) Only one isoform has diff expression pattern and the other remain constant or lowly expressed. 2) One isoform is more abundant than the other but both remain constant during treatment. 3) etc  
# x3: DE_Passed, Strict and Strict Shrinkage
KissDE_DESeq2_Passed <-  DE_Passed[DE_Passed %in% KissDE_Pass_3]
KissDE_DESeq2_Strict_Passed <- DE_Strict_Passed[DE_Strict_Passed %in% KissDE_Pass_3]
KissDE_DESeq2_Shrink_Strict_Passed <- DE_Shrink_Strict_Passed[DE_Shrink_Strict_Passed %in% KissDE_Pass_3]

exported[['resCT1_KissDE_DESeq2_Passed']] <- exported$resCT1[KissDE_DESeq2_Passed,]
exported[['resCT3_KissDE_DESeq2_Passed']] <- exported$resCT3[KissDE_DESeq2_Passed,]
exported[['resT1T3_KissDE_DESeq2_Passed']] <- exported$resT1T3[KissDE_DESeq2_Passed,]

exported[['resCT1_KissDE_DESeq2_Strict_Passed']] <- exported$resCT1[KissDE_DESeq2_Strict_Passed,]
exported[['resCT3_KissDE_DESeq2_Strict_Passed']] <- exported$resCT3[KissDE_DESeq2_Strict_Passed,]
exported[['resT1T3_KissDE_DESeq2_Strict_Passed']] <- exported$resT1T3[KissDE_DESeq2_Strict_Passed,]

exported[['resCT1_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$resCT1[KissDE_DESeq2_Shrink_Strict_Passed,]
exported[['resCT3_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$resCT3[KissDE_DESeq2_Shrink_Strict_Passed,]
exported[['resT1T3_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$resT1T3[KissDE_DESeq2_Shrink_Strict_Passed,]

#-----------------------------------------------------------

exported[['resShrinkCT1_KissDE_DESeq2_Passed']] <- exported$resShrinkCT1[KissDE_DESeq2_Passed,]
exported[['resShrinkCT3_KissDE_DESeq2_Passed']] <- exported$resShrinkCT3[KissDE_DESeq2_Passed,]
exported[['resShrinkT1T3_KissDE_DESeq2_Passed']] <- exported$resShrinkT1T3[KissDE_DESeq2_Passed,]

exported[['resShrinkCT1_KissDE_DESeq2_Strict_Passed']] <- exported$resShrinkCT1[KissDE_DESeq2_Strict_Passed,]
exported[['resShrinkCT3_KissDE_DESeq2_Strict_Passed']] <- exported$resShrinkCT3[KissDE_DESeq2_Strict_Passed,]
exported[['resShrinkT1T3_KissDE_DESeq2_Strict_Passed']] <- exported$resShrinkT1T3[KissDE_DESeq2_Strict_Passed,]

exported[['resShrinkCT1_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$resShrinkCT1[KissDE_DESeq2_Shrink_Strict_Passed,]
exported[['resShrinkCT3_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$resShrinkCT3[KissDE_DESeq2_Shrink_Strict_Passed,]
exported[['resShrinkT1T3_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$resShrinkT1T3[KissDE_DESeq2_Shrink_Strict_Passed,]

#-----------------------------------------------------------

exported[['rlog_KissDE_DESeq2_Passed']] <- exported$rlog[KissDE_DESeq2_Passed,]
exported[['vst_KissDE_DESeq2_Passed']] <- exported$vst[KissDE_DESeq2_Passed,]

exported[['rlog_KissDE_DESeq2_Strict_Passed']] <- exported$rlog[KissDE_DESeq2_Strict_Passed,]
exported[['vst_KissDE_DESeq2_Strict_Passed']] <- exported$vst[KissDE_DESeq2_Strict_Passed,]

exported[['rlog_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$rlog[KissDE_DESeq2_Shrink_Strict_Passed,]
exported[['vst_KissDE_DESeq2_Shrink_Strict_Passed']] <- exported$vst[KissDE_DESeq2_Shrink_Strict_Passed,]

###### 5. Export all results ######

for (z in 1:length(exported)) {
  if(class(exported[[z]]) == "DESeqResults"){
    write.table(as.data.frame(exported[[z]]), file = paste0("./", names(exported[z]), ".txt"), sep = " ")
  }else if(class(exported[[z]]) == "data.frame" | class(exported[[z]]) == "character"){
    write.table(exported[[z]], file = paste0("./", names(exported[z]), ".txt"), sep = " ")
  }else{stop("\n Error04: Check data class")}
}
  
###### 6. Session info ######

writeLines(capture.output(sessionInfo()), "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Logs_SessionInfo/DESeq2_DiffIntegration_sessionInfo.txt")
