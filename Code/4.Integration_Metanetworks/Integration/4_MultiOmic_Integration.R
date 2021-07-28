########################################################################################
############################### 4. MOFA Integration  ###################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######
suppressPackageStartupMessages(library(missForest))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(PCGSE))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(caret))

###### 2. Input Data ######
Isoforms <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.Intersection/Doublediff/3.Doublediff_1_resShrink_vstrlog_git.xlsx", sheet = "vst_DoubDiff_Passed")

###### 3. Impute third replicate with RF ######
## Format 
Isoforms[which(Isoforms$gene != "0"),1] <- paste0(Isoforms$gene[which(Isoforms$gene != "0")], sep = ".", which(Isoforms$gene != "0"))
Isoforms <- Isoforms[,-c(8,9)]
colnames(Isoforms) <- c("KissID","C1", "C2", "T1.1","T1.2", "T3.1", "T3.2" )
rownames(Isoforms) <- Isoforms$KissID
Isoforms_C <- Isoforms[,c(2,3)]
Isoforms_C$C3 <- NA
Isoforms_T1 <- Isoforms[,c(4,5)]
Isoforms_T1$T1.3 <- NA
Isoforms_T3 <- Isoforms[,c(6,7)]
Isoforms_T3$T3.3 <- NA
Isoforms_imputation <- list(C = t(Isoforms_C), T1 = t(Isoforms_T1), T3 = t(Isoforms_T3))

## Code from preprocess_omics_list (Processomics pkg)
Isoforms_imputed <- lapply(Isoforms_imputation,function(matriz2){
  imputed <- matriz2
  registerDoParallel(cores= detectCores() - 1)
  imputed<-(missForestedit(imputed,parallelize="variable"))$ximp
  matriz2 <- imputed
  matriz2
})

## Update25/01/2020
# re-doing cause new annotation and DBs: this time i will use KissID and annotate when its neccesary
Isoforms_MOFA <- read.delim("clipboard", header = T)
Isoforms_MOFA <-  Isoforms_MOFA[which(Isoforms_MOFA$KissID %in% Isoforms_vector$Intersection_doublediff),]

## Re-format and export
Isoforms_MOFA <- t(do.call(rbind, Isoforms_imputed))
rownames(Isoforms_MOFA) <- Isoforms_MOFA$KissID
Isoform.vars <- apply(Isoforms_MOFA[,-1], 1, var)
Isoform.150 <- Isoforms_MOFA[names(sort(Isoform.vars, decreasing = TRUE)[1:150]),] 
Isoform.500 <- Isoforms_MOFA[names(sort(Isoform.vars, decreasing = TRUE)[1:500]),] 
write.table(Isoforms_MOFA, "Isoforms_MOFA.txt", sep = " ")
write.table(Isoform.150, "Isoforms_150.txt", sep = " ")
write.table(Isoform.500, "Isoforms_500.txt", sep = " ")

###### 3. Format and differential test Proteome-Metabolome: Data from Escandon et al., 2021 ######
## Input data
Proteins <- read.delim("clipboard", header = T)
Metabolites <- read.delim("clipboard", header = T)
Escandon.omics <- list(Proteins = Proteins, Metabolites = Metabolites)
Escandon.omics <- lapply(Escandon.omics, function(x){
  rownames(x) <- paste0(x[,1], sep = ".", c(1:nrow(x)))
  x <- x[,-1]
})

## Differential test limma based
Escandon.omics.DE <- lapply(Escandon.omics, function(x){
  # Build models
  phenoData_CT1 <- data.frame(samples = colnames(x)[1:6], Treatment = c(rep("C",3), rep("T1",3)))
  phenoData_CT3 <- data.frame(samples = colnames(x)[c(1:3,7:9)] , Treatment = c(rep("C",3), rep("T3",3)))
  phenoData_T1T3 <- data.frame(samples = colnames(x)[c(4:9)] , Treatment = c(rep("T1",3), rep("T3",3)))
  
  mod_CT1 = model.matrix(~0+as.factor(Treatment), data=phenoData_CT1)
  colnames(mod_CT1) <- c("C", "T1")
  mod_CT3 = model.matrix(~0+as.factor(Treatment), data=phenoData_CT3)
  colnames(mod_CT3) <- c("C", "T3")
  mod_T1T3 = model.matrix(~0+as.factor(Treatment), data=phenoData_T1T3)
  colnames(mod_T1T3) <- c("T1", "T3")
  
  # Fit the model
  fit_CT1<-lmFit(x[,1:6],mod_CT1)
  fit_CT3<-lmFit(x[,c(1:3,7:9)],mod_CT3)
  fit_T1T3<-lmFit(x[,c(4:9)],mod_T1T3)
  
  # Select contrasts
  contrast.matrix_CT1 <- makeContrasts(T1 - C, levels=mod_CT1)
  contrast.matrix_CT3 <- makeContrasts(T3 - C, levels=mod_CT3)
  contrast.matrix_T1T3 <- makeContrasts(T3 - T1, levels=mod_T1T3)
  
  fit2_CT1 <- contrasts.fit(fit_CT1,contrast.matrix_CT1)
  fit2_CT1 <- eBayes(fit2_CT1)
  
  fit2_CT3 <- contrasts.fit(fit_CT3,contrast.matrix_CT3)
  fit2_CT3 <- eBayes(fit2_CT3)
  
  fit2_T1T3 <- contrasts.fit(fit_T1T3,contrast.matrix_T1T3)
  fit2_T1T3 <- eBayes(fit2_T1T3)
  
  # Results and filter by FDR
  DE_CT1_filtered <- limma::topTable(fit2_CT1, number = Inf, coef = 1, p.value = 0.05)
  DE_CT3_filtered <- limma::topTable(fit2_CT3, number = Inf, coef = 1, p.value = 0.05)
  DE_T1T3_filtered <- limma::topTable(fit2_T1T3, number = Inf, coef = 1, p.value = 0.05)
  
  # Select features that passed test
  Feature_names <- Reduce(union, list(rownames(DE_CT1_filtered), 
                                      rownames(DE_CT3_filtered), 
                                      rownames(DE_T1T3_filtered)))
  
  x1 <- x[Feature_names,]
  x1
})

## Split and export
#Only 117 DE prot so top 150 and 500 vars variables from whole dataset
Proteins.vars <- apply(Escandon.omics$Proteins, 1, var)
Proteins.150 <- Escandon.omics$Proteins[names(sort(Proteins.vars, decreasing = TRUE)[1:150]),] 
Proteins.500 <- Escandon.omics$Proteins[names(sort(Proteins.vars, decreasing = TRUE)[1:500]),] 
write.table(Escandon.omics.DE$Proteins, "Proteins_DE.txt", sep = " ")
write.table(Proteins.150, "Proteins_150.txt", sep = " ")
write.table(Proteins.500, "Proteins_500.txt", sep = " ")

Metabolites.vars <- apply(Escandon.omics$Metabolites, 1, var)
Metabolites.150 <- Escandon.omics$Metabolites[names(sort(Metabolites.vars, decreasing = TRUE)[1:150]),] 
Metabolites.500 <- Escandon.omics$Metabolites[names(sort(Metabolites.vars, decreasing = TRUE)[1:500]),] 
write.table(Escandon.omics.DE$Metabolites, "Metabolites_DE.txt", sep = " ")
write.table(Metabolites.150, "Metabolites_150.txt", sep = " ")
write.table(Metabolites.500, "Metabolites_500.txt", sep = " ")

rm(list = ls())
gc()
###### 4. MOFA Integration ######

"C:\Users\bboyl\AppData\Local\Programs\Python\Python38"
reticulate::use_python("C:\\Users\\bboyl\\AppData\\Local\\Programs\\Python\\Python38")

## Load data and create MOFA object | all inputs from MOFA2_inputs

Isoforms_DE <- read.delim("clipboard", header = T, dec = ",")
Proteins_DE <- read.delim("clipboard", header = T, dec = ",")
Metabolites_DE <- read.delim("clipboard", header = T, dec = ",")

Isoforms_150 <- read.delim("clipboard", header = T, dec = ",")
Proteins_150 <- read.delim("clipboard", header = T, dec = ",")
Metabolites_150 <- read.delim("clipboard", header = T, dec = ",")

Isoforms_500 <- read.delim("clipboard", header = T, dec = ",")
Proteins_500 <- read.delim("clipboard", header = T, dec = ",")
Metabolites_500 <- read.delim("clipboard", header = T, dec = ",")

Isoforms_all <- read.delim("clipboard", header = T, dec = ",")
Proteins_all <- read.delim("clipboard", header = T, dec = ",")
Metabolites_all <- read.delim("clipboard", header = T, dec = ",")

MOFA <- list(MOFA_DE = list(Isoforms = Isoforms_DE, Proteins = Proteins_DE, Metabolites = Metabolites_DE),
             MOFA_150 = list(Isoforms = Isoforms_150, Proteins = Proteins_150, Metabolites = Metabolites_150),
             MOFA_500 = list(Isoforms = Isoforms_500, Proteins = Proteins_500, Metabolites = Metabolites_500),
             MOFA_All = list(Isoforms = Isoforms_all, Proteins = Proteins_all, Metabolites = Metabolites_all))

MOFA$MOFA_All <- lapply(MOFA$MOFA_All, function(x){
  x[,1] <- paste0(x[,1], ".", 1:nrow(x))
  x
})

MOFA <- lapply(MOFA, function(y){
  lapply(y, function(x){
    rownames(x) <- x[,1]
    x <- x[,-1]
    as.matrix(x)
  })
})

mymofa <- MOFA$MOFA_500
MOFAobject <- create_mofa(mymofa, groups = NULL)
plot_data_overview(MOFAobject, colors = c("lightsalmon", "darkseagreen", "lightsteelblue2"))

## Define train options and prepare the model

DataOptions <-  get_default_data_options(MOFAobject)
ModelOptions <- get_default_model_options(MOFAobject)
ModelOptions$num_factors <- 3
TrainingOptions <- get_default_training_options(MOFAobject)
TrainingOptions$seed <- 404
TrainingOptions$maxiter <- 60000
TrainingOptions$convergence_mode <- "slow"

MOFAobject <- prepare_mofa(object = MOFAobject, data_options = DataOptions, model_options = ModelOptions,
                          training_options = TrainingOptions)


## Run 25 Random inits during training,compare the models, compare LFs,select best model and export: this dont make sense in Mofa2
n_inits <- 25
MOFAlist <- lapply(1:n_inits, function(it){
  
  MOFAobject <- create_mofa(mymofa, groups = NULL)
  DataOptions <-  get_default_data_options(MOFAobject)
  ModelOptions <- get_default_model_options(MOFAobject)
  ModelOptions$num_factors <- 3
  TrainingOptions <- get_default_training_options(MOFAobject)
  TrainingOptions$seed <- 404 + it
  TrainingOptions$maxiter <- 60000
  TrainingOptions$convergence_mode <- "slow"
  
  MOFAobject <- prepare_mofa(object = MOFAobject, data_options = DataOptions, model_options = ModelOptions,
                             training_options = TrainingOptions)
  
  outfile <- paste0("D:/Splicing_mofa2/500/mofade", it ,".hdf5")
  MOFAobject <- run_mofa(MOFAobject, outfile = outfile)
  MOFAobject
})

compare_elbo(MOFAlist)
compare_factors(MOFAlist)
MOFAobject <- select_model(MOFAlist)

## Save best trained model

#MOFAobject.2 <- MOFA2::load_model("MOFA_EPIOMIDES_DE.hdf5")
MOFAWeights <- get_weights(MOFAobject, as.data.frame = TRUE)
MOFAFactors <- get_factors(MOFAobject, as.data.frame = TRUE)
write.table(MOFAWeights, file = "Weights_MOFA_EPIOMIDES_500.txt", sep = " ")
write.table(MOFAFactors, file = "LFs_MOFA_EPIOMIDES_500.txt", sep = " ")

### DOWNSTREAM ANALYSIS

## add meta-data

sample_metadata <- data.frame(
  sample = samples_names(MOFAobject)[[1]],
  stress = c(rep("C",3), rep("HS",6)),
  treatment = c(rep("C",3), rep("T1",3), rep("T3",3)),
  response = c("H", "H", "H", "R", "R", "R", "H", "H", "H"),
  acclimatation = c("H", "H", "H", "H", "H", "H", "A", "A", "A"),
  sample = c(1:9)
  )

samples_metadata(MOFAobject) <- sample_metadata

## Correlation between factors

plot_factor_cor(MOFAobject)

## Disentangling the heterogeneity. Calculation of variance explained by each factor in each view

MOFAobject@cache$variance_explained$r2_total$group1
MOFAobject@cache$variance_explained$r2_per_factor$group1

a <- plot_variance_explained(MOFAobject, x="view", y="factor")
b <- plot_variance_explained(MOFAobject, x="group", y="factor", plot_total = T)[[2]]
grid.arrange(b, a, ncol = 1, nrow = 2)

## LFs characterization: searching biological meaning of the learned latent factors

# correlation with metadata-covariates: all the coloured ones -- significative association with the meta data

dat <- correlate_factors_with_covariates(MOFAobject, covariates = c("stress", "treatment", "response", "acclimatation", "sample"), alpha = 0.05)

# plot Factors 1 by 1 and matrix 

plot_factor(MOFAobject, 
            factor = 1:2,
            color_by = "treatment",
            shape_by = "stress")

plot_factors(MOFAobject, 
              factor = 1:3,
              color_by = "treatment",
              shape_by = "stress", dot_size = 5)
# LFs of interest: LF1 and LF2

## Weights inspections (LF1-2)

pdf("Topweights_distribution.pdf")
plot_weights(MOFAobject, view = "Isoforms", factors = 1, nfeatures = 10, scale = T)
plot_weights(MOFAobject, view = "Proteins", factors = 1, nfeatures = 10, scale = T)
plot_weights(MOFAobject, view = "Metabolites", factors = 1, nfeatures = 10, scale = T)

plot_weights(MOFAobject, view = "Isoforms", factors = 2, nfeatures = 10, scale = T)
plot_weights(MOFAobject, view = "Proteins", factors = 2, nfeatures = 10, scale = T)
plot_weights(MOFAobject, view = "Metabolites", factors = 2, nfeatures = 10, scale = T)

plot_weights(MOFAobject, view = "Isoforms", factors = 1:2, nfeatures = 10, scale = T)
plot_weights(MOFAobject, view = "Proteins", factors = 1:2, nfeatures = 10, scale = T)
plot_weights(MOFAobject, view = "Metabolites", factors = 1:2, nfeatures = 10, scale = T)
dev.off()

# SideNote LF1 : positive Weights -- molecules more abundance in Control than stress | negative -- stress > control
# SideNote LF2: positive Weights -- molecules more abundance in T1 than C/T3 | negative -- C/T3 > T1

pdf("Topweights.pdf")
plot_top_weights(MOFAobject, 
                 view = c("Proteins", "Isoforms", "Metabolites"),
                 factors = 1,
                 nfeatures = 15,
                 scale = T,
                 abs = T)
plot_top_weights(MOFAobject, 
                 view = c("Proteins", "Isoforms", "Metabolites"),
                 factors = 2,
                 nfeatures = 15,
                 scale = T,
                 abs = T)
plot_top_weights(MOFAobject, 
                 view = c("Proteins", "Isoforms", "Metabolites"),
                 factors = 1:2,
                 nfeatures = 30,
                 scale = T,
                 abs = T)
dev.off()

pdf("DataScatter.pdf")
plot_data_scatter(MOFAobject, view = "Isoforms", features = 6, color_by = "treatment", shape_by = "stress", add_lm = T, factor = 1)
plot_data_scatter(MOFAobject, view = "Proteins", features = 6, color_by = "treatment", shape_by = "stress", add_lm = T, factor = 1)
plot_data_scatter(MOFAobject, view = "Metabolites", features = 6, color_by = "treatment", shape_by = "stress", add_lm = T, factor = 1)

plot_data_scatter(MOFAobject, view = "Isoforms", features = 6, color_by = "treatment", shape_by = "stress", add_lm = T, factor = 2)
plot_data_scatter(MOFAobject, view = "Proteins", features = 6, color_by = "treatment", shape_by = "stress", add_lm = T, factor = 2)
plot_data_scatter(MOFAobject, view = "Metabolites", features = 6, color_by = "treatment", shape_by = "stress", add_lm = T, factor = 2)
dev.off()

df <- MOFAobject@samples_metadata[,c(1,2,3)]
rownames(df) <- df$sample
df$sample <- NULL

pdf("HeatmapWeights.pdf")
plot_data_heatmap(MOFAobject, view = "Isoforms", factor = 1, features = 50, 
                  show_colnames = F, show_rownames = F, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Isoforms", factor = 1, features = 50, 
                  show_colnames = F, show_rownames = T, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Proteins", factor = 1, features = 50, 
                  show_colnames = F, show_rownames = F, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Proteins", factor = 1, features = 50, 
                  show_colnames = F, show_rownames = T, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Metabolites", factor = 1, features = 50, 
                  show_colnames = F, show_rownames = F, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Metabolites", factor = 1, features = 50, 
                  show_colnames = F, show_rownames = T, annotation_samples = df)

plot_data_heatmap(MOFAobject, view = "Isoforms", factor = 2, features = 50, 
                  show_colnames = F, show_rownames = F, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Isoforms", factor = 2, features = 50, 
                  show_colnames = F, show_rownames = T, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Proteins", factor = 2, features = 50, 
                  show_colnames = F, show_rownames = F, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Proteins", factor = 2, features = 50, 
                  show_colnames = F, show_rownames = T, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Metabolites", factor = 2, features = 50, 
                  show_colnames = F, show_rownames = F, annotation_samples = df)
plot_data_heatmap(MOFAobject, view = "Metabolites", factor = 2, features = 50, 
                  show_colnames = F, show_rownames = T, annotation_samples = df)
dev.off()

## IMPORTANT: LF pair inspection, classify samples attending to multi-omic profile

p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "treatment",
                  shape_by = "stress",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed") + theme_bw()
print(p)

### MOFA GSEA
## Mercator bins DB -- Isoforms
bins <- read.delim("clipboard", header = T) #from 2.Heatmap
bin_db <- matrix(0, nrow = 29, ncol = length(bins$bincodes))
bin_names <- read.delim("clipboard", header = T) #from 2.Heatmap
rownames(bin_db) <- bin_names$description
colnames(bin_db) <- bins$KissID


for(i in colnames(bin_db)){
  bn <- as.numeric(bins$bincodes[which(bins$KissID == i)])
  bin_db[bn,i] <- 1
}

# names issue: only for _All
jj <- rownames(MOFAobject@data$Isoforms$group1)
jk <- strsplit(jj, split = ".", fixed = T)
jñ <- lapply(jk, function(x){ paste0(x[-3], collapse = ".")})
jz <- unlist(jñ, use.names = F)
xx <- 0
for(i in 1:ncol(bin_db)){
  if(length(which(jz == colnames(bin_db)[i])) == 0){
    xx <- xx+1
    next
  }
  colnames(bin_db)[i] <- jj[which(jz == colnames(bin_db)[i])]
}

# GSEA on positive weights
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = bin_db, 
                               view = "Isoforms",
                               sign = "positive",
                               factors = c(1,2), min.size = 7
)

# GSEA on negative weights
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = bin_db, 
                               view = "Isoforms",
                               sign = "negative",
                               factors = c(1,2), min.size = 7
)

pdf("GSEAmofa.pdf")
plot_enrichment(res.positive, factor = 1, max.pathways = 10)
plot_enrichment(res.positive, factor = 2, max.pathways = 10)
plot_enrichment(res.negative, factor = 1, max.pathways = 10)
plot_enrichment(res.negative, factor = 2, max.pathways = 10)
plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 1, 
  max.pathways = 3
)
plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 2, 
  max.pathways = 3
)
plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 1, 
  max.pathways = 3
)
plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 2, 
  max.pathways = 3
)
dev.off()

## Mercator bins DB -- Proteins

bins_p <- read.delim("clipboard", header = T) #from 1.MOFA2_inputs
bin_db_p <- matrix(0, nrow = 29, ncol = length(bins_p$bin))
bin_names_p <- read.delim("clipboard", header = T) #from 1.MOFA2_inputs
rownames(bin_db_p) <- bin_names_p$description
colnames(bin_db_p) <- bins_p$ID


for(i in colnames(bin_db_p)){
  bn <- as.numeric(bins$bin[which(bins_p$ID == i)])
  bin_db_p[bn,i] <- 1
}

# names issue: only for _All
jj <- rownames(MOFAobject@data$Proteins$group1)
jk <- strsplit(jj, split = ".", fixed = T)
jñ <- lapply(jk, function(x){x[[1]][1]})
jz <- unlist(jñ, use.names = F)
xx <- 0
for(i in 1:ncol(bin_db_p)){
  if(length(which(jz == colnames(bin_db_p)[i])) == 0){
    xx <- xx+1
    next
  }
  colnames(bin_db_p)[i] <- jj[which(jz == colnames(bin_db_p)[i])]
}


# GSEA on positive weights
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = bin_db_p, 
                               view = "Proteins",
                               sign = "positive",
                               factors = c(1,2), min.size = 7
)

# GSEA on negative weights
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = bin_db_p, 
                               view = "Proteins",
                               sign = "negative",
                               factors = c(1,2), min.size = 7
)

pdf("GSEAmofa_prot_500.pdf")
plot_enrichment(res.positive, factor = 1, max.pathways = 10)
plot_enrichment(res.positive, factor = 2, max.pathways = 10)
plot_enrichment(res.negative, factor = 1, max.pathways = 10)
plot_enrichment(res.negative, factor = 2, max.pathways = 10)
plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 1, 
  max.pathways = 3
)
plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 2, 
  max.pathways = 3
)
plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 1, 
  max.pathways = 3
)
plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 2, 
  max.pathways = 3
)
dev.off()

###### 5. Session info ######

writeLines(capture.output(sessionInfo()), "MOFAIntegration_sessionInfo.txt")

