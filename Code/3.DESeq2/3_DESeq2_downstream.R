########################################################################################
############################ 3. DESeq2  ################################################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 1. Load packages ######
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(cdata))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(nVennR))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(fgsea))

###### 2. Input Data ######

## All: 13080 because the low counts filterig 
DESEq2_CT1 <- readxl::read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_resultsShrink_git.xlsx", sheet = "resShrinkCT1", col_names = T)
DESEq2_CT3 <- readxl::read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_resultsShrink_git.xlsx", sheet = "resShrinkCT3", col_names = T)
DESEq2_T1T3 <- readxl::read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_resultsShrink_git.xlsx", sheet = "resShrinkT1T3", col_names = T)
DESeq2_all <- list(CT1 = DESEq2_CT1, CT3 = DESEq2_CT3, T1T3 = DESEq2_T1T3)

## Differential (< 0.05 padj)
DESEq2_CT1_diff <- readxl::read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_resultsShrink_git.xlsx", sheet = "resShrinkCT1_DE", col_names = T)
DESEq2_CT3_diff <- readxl::read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_resultsShrink_git.xlsx", sheet = "resShrinkCT3_DE", col_names = T)
DESEq2_T1T3_diff <- readxl::read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_resultsShrink_git.xlsx", sheet = "resShrinkT1T3_DE", col_names = T)
DESeq2_diff <- list(CT1 = DESEq2_CT1_diff, CT3 = DESEq2_CT3_diff, T1T3 = DESEq2_T1T3_diff)

## Counts
Kissplice_counts_DESeq2norm <- readxl::read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/2.Heatmap.xlsx", sheet = "NormCounts_DESeq2", col_names = T)

###### 3. Downstream analysis ######

## Remove NA padj isoforms

DESeq2_all <- lapply(DESeq2_all, function(x){
  x <- x[-which(is.na(x$padj)),]
  x
})

## i will re-use the Isoforms_vector list created in 2_KissDE_downstream
## Heatmap with diff expressed isoforms | DESeq2_diff & DESeq2_diff_only

Counts_DESeq2_diff <- Kissplice_counts_DESeq2norm[Kissplice_counts_DESeq2norm$KissID %in% Isoforms_vector$DESeq2_diff,]
Counts_DESeq2_diff_only <- Kissplice_counts_DESeq2norm[Kissplice_counts_DESeq2norm$KissID %in% Isoforms_vector$DESeq2_diff_only,] 

#BINCODES matrix:  col1:"Kiss" col2:"bincodes" | From 2.Heatmap.xlsx
bincodes<-read.delim("clipboard",header=TRUE)
#bindesc matrix structure: rownames=bincode num, col1: description
bindesc<-read.delim("clipboard",header=TRUE)
bindesc$code <- rownames(bindesc)
# matrix
Counts_DESeq2_diff_t <- t(Counts_DESeq2_diff)
Counts_DESeq2_diff_only_t <- t(Counts_DESeq2_diff_only)

Heatmap_matrix <- list(DESeq2_diff = Counts_DESeq2_diff_t, DESeq2_diff_only = Counts_DESeq2_diff_only_t)
Heatmap_matrix <- lapply(Heatmap_matrix, function(matrix){
  colnames(matrix) <- matrix[1,]
  matrix <- matrix[-1,]
  matrix <- apply(X = matrix, MARGIN = 2, FUN = as.numeric)
  C <-apply((matrix[1:4,2:ncol(matrix)]),2,sum)
  T1 <-apply((matrix[5:8,2:ncol(matrix)]),2,sum)
  T3 <-apply((matrix[9:12,2:ncol(matrix)]),2,sum)
  KissID <- colnames(matrix)[-1]
  Treatment_sums<-cbind(C, T1, T3,KissID)
  Treatment_sums <- as.data.frame(Treatment_sums)
  m1 <- merge(Treatment_sums, bincodes, by = "KissID")
  m1[2:(ncol(m1)-1)] <- apply(m1[2:(ncol(m1)-1)], MARGIN = 2, function(x){as.numeric(as.character(x))})
  auxC <-aggregate(C~bincodes,FUN=sum, m1[,-1])
  auxT1<-aggregate(T1~bincodes,FUN=sum, m1[,-1])
  auxT3<-aggregate(T3~bincodes,FUN=sum, m1[,-1])
  aux<-cbind(auxC,auxT1$T1,auxT3$T3)
  colnames(aux)<-c("bin","C","T1","T3")
  auxfsc<-aux[,-1]
  auxsc<-cbind(aux$bin,sweep(auxfsc, 1, rowSums(auxfsc), FUN="/"))
  all<-merge(bindesc, auxsc, by.x=2,by.y=1)
  data<-all[,-1]
  row.names(data)<-(data[,1])
  data<-data[,-1]
  data
})

pdf('normDESeq2_DESeq2_diff_heatmap2.pdf', width = 10, height = 10)
hp<-pheatmap(Heatmap_matrix$DESeq2_diff, color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100),
             clustering_distance_rows=dist(Heatmap_matrix$DESeq2_diff, method = "manhattan"),
             clustering_distance_cols=dist(t(Heatmap_matrix$DESeq2_diff), method = "manhattan"),
             clustering_method="ward.D",display_numbers=TRUE)

dev.off()

pdf('normDESeq2_DESeq2_diff_only_heatmap2.pdf', width = 10, height = 10)
hp<-pheatmap(Heatmap_matrix$DESeq2_diff_only, color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100),
             clustering_distance_rows=dist(Heatmap_matrix$DESeq2_diff_only, method = "manhattan"),
             clustering_distance_cols=dist(t(Heatmap_matrix$DESeq2_diff_only), method = "manhattan"),
             clustering_method="ward.D",display_numbers=TRUE)

dev.off()

## PCA with diff variant isoforms | KissDE_diff & KissDE_diff_only | DESeq2norm & raw

normDESeq2_collapsed <- read.delim("clipboard", header = T, row.names = 1) # only with normDESeq2:collapsed
normDESeq2_collapsed_DESeq2_diff <- normDESeq2_collapsed[Isoforms_vector$DESeq2_diff,] 
normDESeq2_collapsed_DESeq2_diff_only <- normDESeq2_collapsed[Isoforms_vector$DESeq2_diff_only,]

raw <- read.delim("clipboard", header = T, row.names = 1)
raw <- raw[,(ncol(raw)-11):ncol(raw)]
raw_collapsed <- data.frame(C1 = raw$Control_counts1+raw$Control_counts2, 
                            C2 = raw$Control_counts3+raw$Control_counts4,
                            T1.1 = raw$T1_counts5+raw$T1_counts6,
                            T1.2 = raw$T1_counts7+raw$T1_counts8,
                            T3.1 = raw$T3_counts9+raw$T3_counts10,
                            T3.2 = raw$T3_counts11+raw$T3_counts12)
rownames(raw_collapsed) <- rownames(raw)
raw_collapsed_DESeq2_diff <- raw_collapsed[Isoforms_vector$DESeq2_diff,]
raw_collapsed_DESeq2_diff_only <- raw_collapsed[Isoforms_vector$DESeq2_diff_only,]

Stress <- c(rep("Control", 2), rep("Stress", 4))
Treatment <- c(rep("C",2), rep("T1", 2), rep("T3",2))

# dont wanna repeat code just substitute names | too lazy to do looping or apply today
matrix.pca.DESeq2_diff <- prcomp(t(raw_collapsed_DESeq2_diff), center = T, scale. = T, rank. = 5)
fviz_eig(matrix.pca.DESeq2_diff)
matrix.pca.DESeq2_diff_only <- prcomp(t(raw_collapsed_DESeq2_diff_only), center = T, scale. = T, rank. = 5)
fviz_eig(matrix.pca.DESeq2_diff_only)

labelitass <- list(Stress, Treatment)
for(i in 1:length(labelitass)){
  data <- as.data.frame(matrix.pca.DESeq2_diff$x)
  meas_vars <- colnames(data)
  control.Table <- data.frame(expand.grid(meas_vars, meas_vars, stringsAsFactors = FALSE))
  colnames(control.Table) <- c("x", "y")
  control.Table <- cbind(data.frame(pair_key = paste(control.Table[[1]],
                                                     control.Table[[2]], sep = "-"),
                                    stringsAsFactors = FALSE),
                         control.Table)
  data <- cbind(label = labelitass[[i]], data)
  data_aug <- rowrecs_to_blocks(data, control.Table, columnsToCopy = "label")
  splt <- strsplit(data_aug$pair_key, split = "-", fixed = T)
  data_aug$xv <- vapply(splt, function(si) si[[1]], character(1))
  data_aug$yv <- vapply(splt, function(si) si[[2]], character(1))
  data_aug$xv <- factor(as.character(data_aug$xv), meas_vars)
  data_aug$yv <- factor(as.character(data_aug$yv), meas_vars)
  if(length(levels(as.factor(labelitass[[i]]))) >= 6){
    print(ggplot(data_aug, aes(x = x, y = y)) +
            geom_point(aes(color=label, shape ="circle")) +
            facet_grid(yv~xv, labeller = label_both, scales = "free") +
            ggtitle("PCA Biplot Matrix") + 
            ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()))
  }else if(length(levels(as.factor(labelitass[[i]]))) < 6){
    print(ggplot(data_aug, aes(x = x, y = y)) +
            geom_point(aes(color=label, shape =label)) +
            facet_grid(yv~xv, labeller = label_both, scales = "free") +
            ggtitle("PCA Biplot Matrix") + 
            ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()))
  }
  
  a <- readline(prompt = "Enter PCA number: ")
  b <- readline(prompt = "Enter PCA number: ")
  c <- readline(prompt = "Enter PCA number: ")
  LNs <- as.numeric(c(a,b,c))
  par("ask" = TRUE)
  colors <- RColorBrewer::brewer.pal(8, "Set2")
  colors <- colors[as.numeric(as.factor(data$label))]
  print(plot_ly(as.data.frame(data), x = data[,LNs[1]+1], y = data[,LNs[2]+1], z = data[,LNs[3]+1],
                color = data$label, colors = colors) %>%
          add_markers() %>%
          layout(scene = list(xaxis = list(title = paste0("PCA", a)),
                              yaxis = list(title = paste0("PCA", b)),
                              zaxis = list(title = paste0("PCA",c)))))
  
  print(scatterplot3d(x = data[,LNs[1]+1], y = data[,LNs[2]+1], z = data[,LNs[3]+1], color = colors, pch = 16, xlab = paste0("PCA", LNs[1]), ylab = paste0("PCA", LNs[2]), zlab = paste0("PCA", LNs[3]), type = "h"))
  legend("right", legend = data$label[!duplicated(data$label)],
         col =  colors[!duplicated(colors)], pch = 16)
  
}

## Contrasts nVennR/Upset: number of inter-intra contrast specific isoforms

data_DESeq2_diff <- list(CT1 = DESeq2_diff$CT1$KissID, CT3 = DESeq2_diff$CT3$KissID, T1T3 = DESeq2_diff$T1T3$KissID)
data_DESeq2_diff_only <- list(CT1 = DESeq2_diff$CT1$KissID[DESeq2_diff$CT1$KissID %in% Isoforms_vector$DESeq2_diff_only],
                              CT3 = DESeq2_diff$CT3$KissID[DESeq2_diff$CT3$KissID %in% Isoforms_vector$DESeq2_diff_only],
                              T1T3 = DESeq2_diff$T1T3$KissID[DESeq2_diff$T1T3$KissID %in% Isoforms_vector$DESeq2_diff_only])
colors1 <- c("#DFC27D", "#80CDC1", "#D6604D")
colors2 <- c("#BF812D", "#35978F", "#B2182B")

myV1 <- nVennR::plotVenn(data_DESeq2_diff, nCycles = 20000, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile =  "nVennR_3factors_DESeq2_diff1.svg", setColors = colors1)
myV2 <- nVennR::plotVenn(data_DESeq2_diff, nCycles = 20000, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile =  "nVennR_3factors_DESeq2_diff2.svg", setColors = colors2)
myV3 <- nVennR::plotVenn(data_DESeq2_diff, nCycles = 20000, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile =  "nVennR_3factors_DESEq2_diff3.svg")

pdf("UpSet_DESeq2_diff.pdf")
print(upset(fromList(data_DESeq2_diff), sets = names(data_DESeq2_diff), order.by = "freq", 
            keep.order = TRUE, empty.intersections = T, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3"))
dev.off()

tiff("UpSet_DESeq2_diff.tiff")
print(upset(fromList(data_DESeq2_diff), sets = names(data_DESeq2_diff), order.by = "freq", 
            keep.order = TRUE, empty.intersections = T, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3"))
dev.off()

## Kmeans
# DESeq2norm and raw | KissDE_diff and KissDE_diff_only
Kmeans_input <- list(raw_kissde_diff = t(raw_collapsed_DESeq2_diff),
                     raw_kissde_diff_only = t(raw_collapsed_DESeq2_diff_only),
                     norm_kissde_diff = t(normDESeq2_collapsed_DESeq2_diff),
                     norm_kissde_diff_only = t(normDESeq2_collapsed_DESeq2_diff_only))
vectortratamientos <- as.factor(c(rep("C", 2), rep("T1", 2), rep("T3", 2)))
initialcolumn <- 1
initialrow <- 1
data <- lapply(Kmeans_input, function(x){
  datos <- split(as.data.frame(x),vectortratamientos,drop=T)
  datos<-datos[unique(vectortratamientos)]
  datos<-lapply(datos, function(datos) t(apply(datos,1,function(datos) as.numeric(as.character(datos)))))
  treat.means<-lapply(datos, colMeans)
  means.matrix<-do.call("rbind",treat.means)
  colnames(means.matrix)<-colnames(x)
  meansfmatrix<-sweep(means.matrix,2,colSums(means.matrix),"/" )
  meansfmatrix
})

par("ask" = TRUE)
RESULTS.data <-list()
for(i in 1:length(data)){
  clusters <- c(5,25)
  groups<-c()
  groups.ensayo <-c()
  elbow.data<-c()
  lista.kmeans <-list()
  kmeans_matrix <- t(data[[i]])
  z <- 1
  for(j in clusters[1]:clusters[2]){
    set.seed(5881)
    kmeans.result <- kmeans(kmeans_matrix, j,iter.max=30)
    groups <- kmeans.result$cluster #sacamos a que grupo pertenece cada variable
    withinss <- sum(kmeans.result$withinss)
    datosggplot <- as.data.frame(cbind(kmeans_matrix, groups))
    datosggplot <- cbind(datosggplot, rownames(kmeans_matrix))
    colnames(datosggplot)[ncol(datosggplot)] <- "ID"
    result <- list(groups,withinss,datosggplot)
    kmeans.result <- result
    groups.ensayo <- cbind(groups.ensayo,kmeans.result[[1]])
    colnames(groups.ensayo)[z] <- paste(j,"Clusters",sep=" ")
    elbow.data <- c(elbow.data,kmeans.result[[2]])
    names(elbow.data)[z] <-paste(j,"Clusters",sep=" ")
    lista.kmeans[[z]] <-as.data.frame(kmeans.result[[3]])
    names(lista.kmeans)[[z]] <- paste(j,"Clusters",sep=" ")
    z=z+1
  }
  plot(elbow.data, type="o", xlab="Number of clusters, K", ylab="Cumulative within-clusters sum of squares", xaxt = "n")
  axis(1, at=1:length(elbow.data), labels=names(elbow.data))
  resultado <-list(meansfmatrix,lista.kmeans,groups.ensayo,elbow.data)
  names(resultado) <- c("kmeans_matrix","kmeans_list","kmeans_group","elbow.data")
  RESULTS.data[[names(data[i])]] <- resultado
}

fontsizes=c(14,10,16,12)
PLOTS.data <- list()

for(z in 1:length(RESULTS.data)){
  myplotlist <- c()
  for(i in c(12,16)){
    datosggplot<-RESULTS.data[[z]]$kmeans_list[[i]]
    datosggplot <- melt(datosggplot, (ncol(datosggplot)-1):(ncol(datosggplot))) 
    datosggplot <- datosggplot[order(datosggplot$groups),]
    datosggplot2 <- datosggplot
    datosggplot2[,1]<-paste("Cluster",datosggplot2[,1])
    a<-factor(datosggplot2[,1], levels=c(paste("Cluster",unique(datosggplot[,1]),sep=" ")))
    myplot <- ggplot(datosggplot2, aes(x=variable, y=value, group=ID, colour=as.factor(groups), text=datosggplot2$ID)) +
      facet_wrap(facet=a ,scales="free") + 
      geom_line(linetype= "solid",alpha = 0.04) +
      stat_summary(aes(group = groups), fun.y = mean, geom = "line") + 
      theme_minimal()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    if(i == 12){
      myplotlist[[1]] <- myplot  
    }else if(i == 16){
      myplotlist[[2]] <- myplot
    }
  }
  names(myplotlist)<-names(RESULTS.data[[z]]$kmeans_list[c(12,16)])
  PLOTS.data[[names(RESULTS.data[z])]] <- myplotlist   
}

write.table(RESULTS.data$raw_kissde_diff$kmeans_list$`16 Clusters`, file = "kmeans_deseq2_diff_raw_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$raw_kissde_diff$kmeans_list$`20 Clusters`, file = "kmeans_deseq2_diff_raw_20.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$raw_kissde_diff_only$kmeans_list$`16 Clusters`, file = "kmeans_deseq2_diff_only_raw_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$raw_kissde_diff_only$kmeans_list$`20 Clusters`, file = "kmeans_deseq2_diff_only_raw_20.txt", sep = "\t", row.names = F, col.names = T)

write.table(RESULTS.data$norm_kissde_diff$kmeans_list$`16 Clusters`, file = "kmeans_deseq2_diff_norm_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$norm_kissde_diff$kmeans_list$`20 Clusters`, file = "kmeans_deseq2_diff_norm_20.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$norm_kissde_diff_only$kmeans_list$`16 Clusters`, file = "kmeans_deseq2_diff_only_norm_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$norm_kissde_diff_only$kmeans_list$`20 Clusters`, file = "kmeans_deseq2_diff_only_norm_20.txt", sep = "\t", row.names = F, col.names = T)

## GO enrichments and REVIGO
# several advices: fill NAs with generic GOs; use all dfs; use isoforms instead of events

Splicing_annotation <- read.delim("clipboard", header = T)
GO_enrichments <- list(CT1 = list(), CT3 = list(), T1T3 = list())

for(i in 1:length(DESeq2_all)){
  p.value.column <- grep("padj", colnames(DESeq2_all[[i]]))
  p.value <- c(as.data.frame(DESeq2_all[[i]])[,p.value.column])
  isoforms <- DESeq2_all[[i]]$KissID
  # Parse GOs
  GOs <- list()
  for(j in 1:nrow(DESeq2_all[[i]])){
    GOs[[j]] <- Splicing_annotation$GO_gotubro[which(Splicing_annotation$KissID == isoforms[j])]
  }
  GO <- as.character(unlist(GOs, use.names = F))
  df <- tibble(KissID = isoforms, GO = GO, p.value = p.value)
  df2 <- df %>%
    mutate(GO=strsplit(GO, " | ", fixed = T)) %>% unnest(GO)
  # Correct errors -- several GOs together (DBsyntaxerrors)
  errors <- df2[which(nchar(df2$GO) == 20 & substr(df2$GO, 1,3) == "GO:"),]
  df2 <- df2[-which(nchar(df2$GO) == 20 & substr(df2$GO, 1,3) == "GO:"),]
  errors2 <- errors %>% 
    mutate(GO=strsplit(GO, "GO:", fixed = T)) %>% unnest(GO)
  errors2 <- errors2[which(errors2$GO != ""),]
  errors2$GO <- paste0("GO:", errors2$GO)
  df2 <- rbind(df2, errors2)
  # delete all not NA;GO values; i.e: "", NOT etc
  checks <- df2[which(!is.na(df2$GO) & substr(df2$GO, 1,3) != "GO:"),]
  df2 <- df2[-which(!is.na(df2$GO) & substr(df2$GO, 1,3) != "GO:"),]
  checks2 <- checks %>%
    mutate(GO=strsplit(GO, "| ", fixed = T)) %>% unnest(GO)
  checks2 <- checks2[which(substr(checks2$GO,1,3) == "GO:"),]
  df2 <- rbind(df2, checks2)
  # loop filling all NAs with generic GOs
  Class <- list(MF = NA, BP = NA, CC = NA)
  for(z in 1:length(Class)){
    df3 <- df2
    if(names(Class[z]) == "MF"){
      df3$GO[which(is.na(df3$GO))] <- "GO:0003674"
      Class[[z]] <- df3
    }else if(names(Class[z]) == "BP"){
      df3$GO[which(is.na(df3$GO))] <- "GO:0008150"
      Class[[z]] <- df3
    }else if(names(Class[z]) == "CC"){
      df3$GO[which(is.na(df3$GO))] <- "GO:0005575"
      Class[[z]] <- df3
    }
  }
  # Breaker check
  Breaker.check <- lapply(Class, function(x){
    if(length(which((substr(x$GO,1,3) == "GO:"))) != nrow(x)){ stop(" You are doing something wrong boy")} 
  })
  # export to list
  GO_enrichments[[names(DESeq2_all[i])]] <- Class
}

# export-filter tables and  do TreeMap with REVIGO - http://revigo.irb.hr/

for(i in 1:length(GO_enrichments)){
  for(j in 1:length(GO_enrichments[[i]])){
    write.table(GO_enrichments[[i]][[j]], file = paste0("GO_enrichments_", names(GO_enrichments[i]), "_",names(GO_enrichments[[i]][j]), ".txt"), sep = "\t", col.names = T, row.names = F)
  }
}

GO_enrichments_DESeq2_diff <- GO_enrichments
GO_enrichments_DESeq2_diff_only <- GO_enrichments

for(i in 1:length(GO_enrichments_DESeq2_diff)){
  for(j in 1:length(GO_enrichments_DESeq2_diff[[i]])){
    GO_enrichments_DESeq2_diff[[i]][[j]] <- GO_enrichments_DESeq2_diff[[i]][[j]][which(GO_enrichments_DESeq2_diff[[i]][[j]]$KissID %in% Isoforms_vector$DESeq2_diff),]
    write.table(GO_enrichments_DESeq2_diff[[i]][[j]], file = paste0("GO_enrichments_DESeq2d_", names(GO_enrichments_DESeq2_diff[i]), "_",names(GO_enrichments_DESeq2_diff[[i]][j]), ".txt"), sep = "\t", col.names = T, row.names = F)
  }
}

for(i in 1:length(GO_enrichments_DESeq2_diff_only)){
  for(j in 1:length(GO_enrichments_DESeq2_diff_only[[i]])){
    GO_enrichments_DESeq2_diff_only[[i]][[j]] <- GO_enrichments_DESeq2_diff_only[[i]][[j]][which(GO_enrichments_DESeq2_diff_only[[i]][[j]]$KissID %in% Isoforms_vector$DESeq2_diff_only),]
    write.table(GO_enrichments_DESeq2_diff_only[[i]][[j]], file = paste0("GO_enrichments_DESeq2do_", names(GO_enrichments_DESeq2_diff_only[i]), "_",names(GO_enrichments_DESeq2_diff_only[[i]][j]), ".txt"), sep = "\t", col.names = T, row.names = F)
  }
}

# fgsea bin functional enrichments 

bins <- read.delim("clipboard", header = T) #vst_KissDE
bin.names <- read.delim("clipboard", header = T)#2.Heatmap-Bin Description

DESeq2_diff_only <- lapply(DESeq2_diff, function(x){
  do <- x[which(x$KissID %in% Isoforms_vector$DESeq2_diff_only),]
  do
})
FGSEA <- DESeq2_all
Bin_db <- list()
for(i in factor(rownames(bin.names))){
  Bin_db[[i]] <- bins$KissID[which(bins$bincodes == i)]
}
names(Bin_db) <- as.character(bin.names$description)
FGSEA_res <- FGSEA
for(i in 1:length(FGSEA)){
  for(z in 1:length(FGSEA[[i]])){
    datas <- FGSEA[[i]][[z]]
    datas$bin <- NA
    datas$stat <- NA
    for(j in 1:nrow(datas)){
      datas$bin[j] <- bins$bincodes[which(bins$KissID == datas$KissID[j])]
      datas$stat[j] <- datas$log2FoldChange[j]/datas$lfcSE[j] 
    }
    if(length(which(is.na(datas$padj))) != 0){ stop("Stop boy, you doing something wrong")}
    FGSEA[[i]][[z]] <- datas
    datas <- datas[,c(1,8)]
    ranks <- deframe(datas)
    fGseaRES <- fgsea(pathways = Bin_db, stats = ranks, nperm = 100000, minSize = 7, maxSize = 1000)
    fgseaResTidy <- fGseaRES %>%
      as_tibble() %>%
      arrange(desc(NES))
    fgseaResTidy %>% 
      dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
      arrange(padj) %>% 
      DT::datatable()
    FGSEA_res[[i]][[z]] <- fgseaResTidy
    write.table(fgseaResTidy[,-8], file = paste0("binfgsea_DESeq2_", names(FGSEA[i]), "_", names(FGSEA[[i]][z]), ".txt"), sep = "\t", col.names = T)
    pdf(file = paste0("binfgsea_DESeq2_", names(FGSEA[i]), "_", names(FGSEA[[i]][z]), ".pdf"))
    
    ei1 <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj <0.05)) + coord_flip() +
      labs(x="Mercator Bins", y="Normalized Enrichment Stat", title=paste0("Preranked GSEA (",  names(FGSEA[[i]][z]), ")")) + theme_minimal()
    
    ei2 <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=NES)) + coord_flip() +
      labs(x="Mercator Bins", y="Normalized Enrichment Stat", title=paste0("Preranked GSEA (",  names(FGSEA[[i]][z]), ")")) + theme_minimal()
    
    print(ei2 + scale_fill_gradient2(low = "slateblue2", high = "darkred")) 
    print(ei1 + scale_fill_manual(values = c("lightgoldenrod3", "lightskyblue4")))
    
    dev.off()
  }  
}

# Volcano plots + new annotation

pal <- c( "lightcoral", "lightgoldenrod1", "lightsteelblue", "lightskyblue4")
for(i in 1:length(DESeq2_all)){
  datas <- DESeq2_all[[i]]
  datas$Symbol <- NA
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
 datas$logpvalue <- -log10(datas$padj)
 datas$group <- "NotSignificant"
 datas$group[which(datas$padj < 0.05 & abs(datas$log2FoldChange) < 1.8)] <- "Significant"
 datas$group[which(datas$padj > 0.05 & abs(datas$log2FoldChange) > 1.8)] <- "FoldChange"
 datas$group[which(datas$padj < 0.05 & abs(datas$log2FoldChange) > 1.8)] <- "Significant&FoldChange"
 top_peaks_down <-  datas[with(datas, order(log2FoldChange, padj, decreasing = FALSE)),][1:5,]
 top_peaks_up <-  datas[with(datas, order(log2FoldChange, padj, decreasing = TRUE)),][1:5,]
 top_peaks <- rbind(top_peaks_up, top_peaks_down)
 write.table(top_peaks, file = paste0("volcano_DESeq2_all_", names(DESeq2_all[i]), ".txt"), sep = "\t", col.names = T)
 a <- list()
 for(j in seq_len(nrow(top_peaks))){
   m <- top_peaks[j,]
   a[[j]] <- list(
     x = m[["log2FoldChange"]],
     y = m[["logpvalue"]],
     text = m[["Symbol"]],
     xref = "x",
     yref = "y",
     showarrow = TRUE,
     arrowhead = 0.3,
     ax = -10,
     ay = 50,
     size = 0.5 
   )
 }
 p1 <- plot_ly(data = datas, x = datas$log2FoldChange, y = datas$logpvalue, text = datas$Symbol, mode = "markers", color = datas$group, colors = pal, xaxis = "FoldChange[log2]",  yaxis = "p-value[-log10]") %>% 
   layout(title =names(DESeq2_all[i])) %>% 
   layout(annotations = a)
 
 export(p = p1, file = paste0("volcano_DESeq2_all_", names(DESeq2_all[i]), "_woCoherency.pdf"))
}
