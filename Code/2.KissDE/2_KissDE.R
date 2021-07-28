########################################################################################
############################ 2. KISSDE  ################################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######
suppressPackageStartupMessages(library(kissDE))

###### 2. Input Data ######
## Condition vectors: 

# (1) Control (C) vs HeatStress (1-Day,T1 & 3-Day,T3)
myCondition_1  <- c(rep("Control_Condition1",4), rep("HeatStress_Condition2",8))
# (2) Control (C) vs Short-term Stress (1-Day,T1)
myCondition_2 <-  c(rep("Control_Condition1",4), rep("HeatStress_T1_Condition2",4), rep("*",4))
# (3) Control (C) vs Long/Medium-term Stress (3-Day,T3)
myCondition_3 <- c(rep("Control_Condition1",4), rep("*",4), rep("HeatStress_T3_Condition2",4))
# (4) Short-term (1-Day, T1) vs Long/Medium-term Stress (3-Day,T3)
myCondition_4 <-  c(rep("*",4), rep("HeatStress_T1_Condition2",4), rep("HeatStress_T3_Condition2",4))

## Obtain raw counts from KisSplice output

fpath <- "D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/results_k51_type1.fa"
myCounts <- kissplice2counts(fpath, counts = 2, pairedEnd = TRUE)
write.table(myCounts$countsEvents, "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/Raw_counts.txt", sep = " ")

###### 2. Analysis ######
## Quality-Control

qualityControl(myCounts, myCondition_1)
qualityControl(myCounts, myCondition_2)
qualityControl(myCounts, myCondition_3)
qualityControl(myCounts, myCondition_4)

## Differential variant analysis

myResults1 <- diffExpressedVariants(countsData = myCounts,
                                    conditions = myCondition_1, pvalue = 0.05,
                                    filterLowCountsVariants = 10, flagLowCountsConditions = 10,
                                    technicalReplicates = FALSE, nbCore = 3)

myResults2 <- diffExpressedVariants(countsData = myCounts,
                                    conditions = myCondition_2, pvalue = 0.05,
                                    filterLowCountsVariants = 10, flagLowCountsConditions = 10,
                                    technicalReplicates = FALSE, nbCore = 3)

myResults3 <- diffExpressedVariants(countsData = myCounts,
                                    conditions = myCondition_3, pvalue = 0.05,
                                    filterLowCountsVariants = 10, flagLowCountsConditions = 10,
                                    technicalReplicates = FALSE, nbCore = 3)

myResults4 <- diffExpressedVariants(countsData = myCounts,
                                    conditions = myCondition_4, pvalue = 0.05,
                                    filterLowCountsVariants = 10, flagLowCountsConditions = 10,
                                    technicalReplicates = FALSE, nbCore = 3)

###### 3. Export ######

writeOutputKissDE(myResults1, output = "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/kissDE_Control_HeatStress_k51.tab", adjPvalMax = 0.05, dPSImin = 0)
writeOutputKissDE(myResults2, output = "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/kissDE_Control_T1_k51.tab", adjPvalMax = 0.05, dPSImin = 0)
writeOutputKissDE(myResults3, output = "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/kissDE_Control_T3_k51.tab", adjPvalMax = 0.05, dPSImin = 0)
writeOutputKissDE(myResults4, output = "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/2.KissDE/kissDE_T1_T3_k51.tab", adjPvalMax = 0.05, dPSImin = 0)

###### 4. Session info ######

writeLines(capture.output(sessionInfo()), "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Logs_SessionInfo/KissDE_sessionInfo.txt")
