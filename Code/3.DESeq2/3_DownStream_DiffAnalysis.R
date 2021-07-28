########################################################################################
########################### 3. Down-Stream Diff Analysis  ##############################
########################################################################################

setwd("D:\\AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/")

###### 0. SideNote Input ######
# I will try to do all the DoubleDiff analysis with DoubleDiff_Passed isoforms (944)

###### 1. Load packages ######
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cdata))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(DiagrammeR))
suppressPackageStartupMessages(library(nVennR))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(fgsea))

###### 2. Principal Component Analysis  ######
# Input data: DoubDiff_SS_Passed vst data | DoubleDiff_Passed vst data | and the rest types of data with vst  
DoubDiff_vst <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.Intersection/Doublediff/3.Doublediff_1_resShrink_vstrlog_git.xlsx",
                       sheet = "vst_DoubDiff_Passed") 
Global_vst <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_rlog_vst_git.xlsx",
                         sheet = "vst") 
KissDE_diff_vst <- Global_vst[which(Global_vst$KissID %in% Isoforms_vector$KissDE_diff),]
KissDE_diff_only_vst <- Global_vst[which(Global_vst$KissID %in% Isoforms_vector$KissDE_diff_only),]
DESeq2_diff_vst <- Global_vst[which(Global_vst$KissID %in% Isoforms_vector$DESeq2_diff),]
DESeq2_diff_only_vst <- Global_vst[which(Global_vst$KissID %in% Isoforms_vector$DESeq2_diff_only),]

pc <- prcomp(t(DESeq2_diff_only_vst[,-c(1,8,9,10,11)]))
percent <- round(summary(pc)$importance[2,]*100, 1)
samples <- c(rep("C",2), rep("T1",2), rep("T3", 2))

# Color by factor
treatment <- c(rep("Control", 2), rep("HeatStress", 4))
treatment.cols <- c("pink", "lightblue")
condition.cols <- c("darkseagreen1", "darksalmon", "darkred")

#  Beeswarm PCbiplot matrix
data <- as.data.frame(pc$x)
meas_vars <- colnames(data)
control.Table <- data.frame(expand.grid(meas_vars, meas_vars, stringsAsFactors = FALSE))
colnames(control.Table) <- c("x", "y")
control.Table <- cbind(data.frame(pair_key = paste(control.Table[[1]],control.Table[[2]], sep = "-"),
                                  stringsAsFactors = FALSE),control.Table)
data1 <- cbind(samples,data)
data2 <- cbind(treatment,data)

data1_aug <- rowrecs_to_blocks(data1, control.Table, columnsToCopy = "samples")
data2_aug <- rowrecs_to_blocks(data2, control.Table, columnsToCopy = "treatment")

splt1 <- strsplit(data1_aug$pair_key, split = "-", fixed = T)
splt2 <- strsplit(data1_aug$pair_key, split = "-", fixed = T)

data1_aug$xv <- vapply(splt1, function(si) si[[1]], character(1))
data1_aug$yv <- vapply(splt1, function(si) si[[2]], character(1))
data1_aug$xv <- factor(as.character(data1_aug$xv), meas_vars)
data1_aug$yv <- factor(as.character(data1_aug$yv), meas_vars)

data2_aug$xv <- vapply(splt2, function(si) si[[1]], character(1))
data2_aug$yv <- vapply(splt2, function(si) si[[2]], character(1))
data2_aug$xv <- factor(as.character(data2_aug$xv), meas_vars)
data2_aug$yv <- factor(as.character(data2_aug$yv), meas_vars)

samples_Beeswarm <-  ggplot(data1_aug, aes(x = x, y = y)) +
  geom_point(aes(color=samples, shape =samples)) +
  facet_grid(yv~xv, labeller = label_both, scales = "free") +
  ggtitle("PCA Biplot Matrix") + 
  ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

treatment_Beeswarm <-  ggplot(data2_aug, aes(x = x, y = y)) +
  geom_point(aes(color=treatment, shape =treatment)) +
  facet_grid(yv~xv, labeller = label_both, scales = "free") +
  ggtitle("PCA Biplot Matrix") + 
  ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

scatterplot3d(data[,1],
              data[,2],
              data[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=condition.cols[as.integer(factor(samples))],
              pch=19)
legend("bottomright",pch=c(15,15,15),col=c(NA,condition.cols[1:3]),
       legend=c("Color:", levels(factor(samples))))

rm(list = ls())
gc()

###### 3. Heatmaps ######
# Input data: DoubleDiff_Passed normDESeq2 and vst

Kissplice_counts_DESeq2norm_collapsed <- read.delim("clipboard", header = T) # from 1.Kissplice/2.Heatmap
doublediff_vst <- DoubDiff_vst[,-c(8,9,10,11)]
doublediff_norm <-Kissplice_counts_DESeq2norm_collapsed[Kissplice_counts_DESeq2norm_collapsed$KissID %in% Isoforms_vector$Intersection_doublediff,] 

#BINCODES matrix:  col1:"Kiss" col2:"bincodes" | From 2.Heatmap.xlsx
bincodes<-read.delim("clipboard",header=TRUE)
#bindesc matrix structure: rownames=bincode num, col1: description
bindesc<-read.delim("clipboard",header=TRUE)
bindesc$code <- rownames(bindesc)
# matrix
doublediff_vst_t <- t(doublediff_vst)
doublediff_norm_t <- t(doublediff_norm)

Heatmap_matrix_dd <- list(norm = doublediff_norm_t)
Heatmap_matrix_dd <- lapply(Heatmap_matrix_dd, function(matrix){
  colnames(matrix) <- matrix[1,]
  matrix <- matrix[-1,]
  matrix <- apply(X = matrix, MARGIN = 2, FUN = as.numeric)
  C <-apply((matrix[1:2,2:ncol(matrix)]),2,sum)
  T1 <-apply((matrix[3:4,2:ncol(matrix)]),2,sum)
  T3 <-apply((matrix[5:6,2:ncol(matrix)]),2,sum)
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

pdf('vst_doublediff_heatmap.pdf', width = 10, height = 10)
hp<-pheatmap(Heatmap_matrix$vst, color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100),
             clustering_distance_rows=dist(Heatmap_matrix$vst, method = "manhattan"),
             clustering_distance_cols=dist(t(Heatmap_matrix$vst), method = "manhattan"),
             clustering_method="ward.D",display_numbers=TRUE)

dev.off()

pdf('normDESeq2_doublediff_heatmap.pdf', width = 10, height = 10)
hp<-pheatmap(Heatmap_matrix$norm, color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100),
             clustering_distance_rows=dist(Heatmap_matrix$norm, method = "manhattan"),
             clustering_distance_cols=dist(t(Heatmap_matrix$norm), method = "manhattan"),
             clustering_method="ward.D",display_numbers=TRUE)

dev.off()

###### 4. Ven diff & Upset ######
# Input data: DoubleDiff_Passed and all diff ones (8) global Upset

data_doublediff_passed <- list(CT1 = DESeq2_diff$CT1$KissID[DESeq2_diff$CT1$KissID %in% Isoforms_vector$Intersection_doublediff], 
                               CT3 = DESeq2_diff$CT3$KissID[DESeq2_diff$CT3$KissID %in% Isoforms_vector$Intersection_doublediff], 
                               T1T3 = DESeq2_diff$T1T3$KissID[DESeq2_diff$T1T3$KissID %in% Isoforms_vector$Intersection_doublediff])
data_all_diff <- list(CHS_DV = c(paste0(KissDE_diff$CHS$`1.ID`, ".1"), paste0(KissDE_diff$CHS$`1.ID`, ".2")),
                              CT1_DV = c(paste0(KissDE_diff$CT1$`1.ID`, ".1"), paste0(KissDE_diff$CT1$`1.ID`, ".2")),
                              CT3_DV = c(paste0(KissDE_diff$CT3$`1.ID`, ".1"), paste0(KissDE_diff$CT3$`1.ID`, ".2")),
                              T1T3_DV = c(paste0(KissDE_diff$T1T3$`1.ID`, ".1"), paste0(KissDE_diff$T1T3$`1.ID`, ".2")),
                     CT1_DE = DESeq2_diff$CT1$KissID ,  CT3_DE = DESeq2_diff$CT3$KissID, T1T3_DE = DESeq2_diff$T1T3$KissID)
colors1 <- c("#DFC27D", "#80CDC1", "#D6604D")
colors2 <- c("#BF812D", "#35978F", "#B2182B")

myV1 <- nVennR::plotVenn(data_all_diff[-1], nCycles = 60000, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile =  "nVennR_6factors_alldiff1.svg", setColors = colors1)
myV2 <- nVennR::plotVenn(data_all_diff[-1], nCycles = 60000, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile =  "nVennR_6factors_alldif2.svg", setColors = colors2)
myV3 <- nVennR::plotVenn(data_all_diff[-1], nCycles = 60000, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile =  "nVennR_6factors_alldif3.svg")

pdf("UpSet_data_all_diff6.pdf")
print(upset(fromList(data_all_diff[-1]), sets = names(data_all_diff[-1]), order.by = "freq", 
            keep.order = TRUE, empty.intersections = T, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3"))
dev.off()

tiff("UpSet_data_all_diff6.tiff")
print(upset(fromList(data_all_diff[-1]), sets = names(data_all_diff[-1]), order.by = "freq", 
            keep.order = TRUE, empty.intersections = T, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3"))
dev.off()

###### 5. K-means ######
# DESeq2norm and raw | doublediff
Kmeans_input <- list(raw_doublediff = t(raw_collapsed[Isoforms_vector$Intersection_doublediff,]),
                     norm_doublediff = t(normDESeq2_collapsed[Isoforms_vector$Intersection_doublediff,]))
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
  for(i in c(8,12,16)){
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
      myplotlist[[2]] <- myplot  
    }else if(i == 16){
      myplotlist[[3]] <- myplot
    }else if(i == 8){
      myplotlist[[1]] <- myplot
    }
  }
  names(myplotlist)<-names(RESULTS.data[[z]]$kmeans_list[c(8,12,16)])
  PLOTS.data[[names(RESULTS.data[z])]] <- myplotlist   
}

write.table(RESULTS.data$raw_doublediff$kmeans_list$`12 Clusters`, file = "kmeans_doublediff_raw_12.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$raw_doublediff$kmeans_list$`16 Clusters`, file = "kmeans_doublediff_raw_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$raw_doublediff$kmeans_list$`20 Clusters`, file = "kmeans_doublediff_raw_20.txt", sep = "\t", row.names = F, col.names = T)

write.table(RESULTS.data$norm_doublediff$kmeans_list$`12 Clusters`, file = "kmeans_doublediff_norm_12.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$norm_doublediff$kmeans_list$`16 Clusters`, file = "kmeans_doublediff_norm_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$norm_doublediff$kmeans_list$`20 Clusters`, file = "kmeans_doublediff_norm_20.txt", sep = "\t", row.names = F, col.names = T)

###### 6.  GO enrichments and REVIGO ######
# padj | KissDE DV isoforms with DESeq2 stats | This is because very low N in doublediff (944; 93 only T1vsT3) and in order to have some agreement between bins;GOs enrichments and volcanos

Splicing_annotation <- read.delim("clipboard", header = T)
GO_enrichments <- list(CT1 = list(), CT3 = list(), T1T3 = list())
doublediff_passed <- list(CT1 = DESeq2_diff$CT1[DESeq2_diff$CT1$KissID %in% Isoforms_vector$Intersection_doublediff,], 
                          CT3 = DESeq2_diff$CT3[DESeq2_diff$CT3$KissID %in% Isoforms_vector$Intersection_doublediff,], 
                          T1T3 = DESeq2_diff$T1T3[DESeq2_diff$T1T3$KissID %in% Isoforms_vector$Intersection_doublediff,])
intersection <- list(CT1 = DESeq2_all$CT1[DESeq2_all$CT1$KissID %in% Isoforms_vector$KissDE_diff,],
                     CT3 = DESeq2_all$CT3[DESeq2_all$CT3$KissID %in% Isoforms_vector$KissDE_diff,],
                     T1T3 = DESeq2_all$T1T3[DESeq2_all$T1T3$KissID %in% Isoforms_vector$KissDE_diff,])

for(i in 1:length(intersection)){
  p.value.column <- grep("padj", colnames(intersection[[i]]))
  p.value <- c(as.data.frame(intersection[[i]])[,p.value.column])
  isoforms <- intersection[[i]]$KissID
  # Parse GOs
  GOs <- list()
  for(j in 1:nrow(intersection[[i]])){
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
  GO_enrichments[[names(intersection[i])]] <- Class
}

# export-filter tables and  do TreeMap with REVIGO - http://revigo.irb.hr/

for(i in 1:length(GO_enrichments)){
  for(j in 1:length(GO_enrichments[[i]])){
    write.table(GO_enrichments[[i]][[j]], file = paste0("Intersection_GO_enrichments_", names(GO_enrichments[i]), "_",names(GO_enrichments[[i]][j]), ".txt"), sep = "\t", col.names = T, row.names = F)
  }
}


###### 7. FGSEA bin enrichments ###### 
# stat | KissDE DV isoforms with DESeq2 stats | This is because very low N in doublediff (944; 93 only T1vsT3) and in order to have some agreement between bins;GOs enrichments and volcanos

bins <- read.delim("clipboard", header = T) #vst_KissDE
bin.names <- read.delim("clipboard", header = T)#2.Heatmap-Bin Description

FGSEA <- intersection
Bin_db <- list()
for(i in factor(rownames(bin.names))){
  Bin_db[[i]] <- bins$KissID[which(bins$bincodes == i)]
}
names(Bin_db) <- as.character(bin.names$description)
FGSEA_res <- FGSEA
for(i in 1:length(FGSEA)){
  for(z in 1:length(FGSEA[[i]])){
    datas <- FGSEA[[i]]
    datas$bin <- NA
    datas$stat <- NA
    for(j in 1:nrow(datas)){
      datas$bin[j] <- bins$bincodes[which(bins$KissID == datas$KissID[j])]
      datas$stat[j] <- datas$log2FoldChange[j]/datas$lfcSE[j] 
    }
    if(length(which(is.na(datas$padj))) != 0){ stop("Stop boy, you doing something wrong")}
    FGSEA[[i]] <- datas
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
    FGSEA_res[[i]] <- fgseaResTidy
    write.table(fgseaResTidy[,-8], file = paste0("binfgsea_intersection_", names(FGSEA[i]), ".txt"), sep = "\t", col.names = T)
    pdf(file = paste0("binfgsea_intersection_", names(FGSEA[i]), ".pdf"))
    
    ei1 <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj <0.05)) + coord_flip() +
      labs(x="Mercator Bins", y="Normalized Enrichment Stat", title=paste0("Preranked GSEA (",  names(FGSEA[i]), ")")) + theme_minimal()
    
    ei2 <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=NES)) + coord_flip() +
      labs(x="Mercator Bins", y="Normalized Enrichment Stat", title=paste0("Preranked GSEA (",  names(FGSEA[i]), ")")) + theme_minimal()
    
    print(ei2 + scale_fill_gradient2(low = "slateblue2", high = "darkred")) 
    print(ei1 + scale_fill_manual(values = c("lightgoldenrod3", "lightskyblue4")))
    
    dev.off()
  }  
}

###### 8. Volcano Plots ######
# logFC and padj | KissDE DV isoforms with DESeq2 stats | This is because very low N in doublediff (944; 93 only T1vsT3) and in order to have some agreement between bins;GOs enrichments and volcanos

pal <- c( "lightcoral", "lightgoldenrod1", "lightsteelblue", "lightskyblue4")
for(i in 1:length(intersection)){
  datas <- intersection[[i]]
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
  write.table(top_peaks, file = paste0("volcano_intersection_", names(intersection[i]), ".txt"), sep = "\t", col.names = T)
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
    layout(title =names(intersection[i])) %>% 
    layout(annotations = a)
  
  export(p = p1, file = paste0("volcano_intersection_", names(intersection[i]), "_woCoherency.pdf"))
}

###### 9. KissDE, DESeq2 and DoubDiff VennDiagrams | FlowChart ######
## FlowChart
#1 Obtain desired numbers
# KISSPLICE
KisSplice <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.k51_Annotation_SYMMETRIC.xlsx")
Isoforms_KisSplice <- nrow(KisSplice)
Events_KisSplice <- Isoforms_KisSplice/2
Isoform_names_KisSplice <- paste0(KisSplice$events.names[order(KisSplice$events.names)], sep = ".", c(1,2))
Events_names_KisSplice <- unique(KisSplice$events.names)

# KISSDE
KissDE <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/2.KissDE_DESeq2_all.xlsx")
Isoforms_KissDE <- nrow(KissDE)
Isoform_names_KissDE <- KissDE$KissID
Events_KissDE_c <- list()
for(i in 1:nrow(KissDE)){
  Events_KissDE_c[[i]] <-  strsplit(KissDE$KissID[i], split = ".", fixed = T)[[1]][1]
}
Events_KissDE <- length(unique(unlist(Events_KissDE_c, use.names = F)))
Events_names_KissDE <- unique(unlist(Events_KissDE_c, use.names = F))

# NOT-KISSDE
Isoform_names_notKissDE <- Isoform_names_KisSplice[! Isoform_names_KisSplice %in% Isoform_names_KissDE]  
Isoform_notKissDE <- length(Isoform_names_notKissDE)
Events_names_notKissDE <- Events_names_KisSplice[! Events_names_KisSplice %in% Events_names_KissDE]
Events_notKissDE <- length(Events_names_notKissDE)

#DE
DE <- read_excel("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/3.DESeq2/1.DESeq2_rlog_vst.xlsx", sheet = "vst_DE_final")
Isoform_names_DE <- DE$KissID
Isoform_DE <- length(Isoform_names_DE)
Events_DE_c <- list()
for(i in 1:nrow(DE)){
  Events_DE_c[[i]] <-  strsplit(DE$KissID[i], split = ".", fixed = T)[[1]][1]
}
Events_names_DE <-  unique(unlist(Events_DE_c, use.names = F))
Events_DE <- length(Events_names_DE)

#DE-KissDE
Events_names_DE_KissDE <- Events_names_KissDE[Events_names_KissDE %in% Events_names_DE]
Events_DE_KissDE <- length(Events_names_DE_KissDE)
Isoform_names_DE_KissDE <- Isoform_names_KissDE[Isoform_names_KissDE %in% Isoform_names_DE]
Isoform_DE_KissDE <- length(Isoform_names_DE_KissDE)

# KissDE only
Events_names_KissDEonly <- Events_names_KissDE[! Events_names_KissDE %in% Events_names_DE]
Events_KissDEonly <- length(Events_names_KissDEonly)
Isoform_names_KissDEonly <- Isoform_names_KissDE[! Isoform_names_KissDE %in% Isoform_names_DE]
Isoform_KissDEonly <- length(Isoform_names_KissDEonly)

# DE-only
Events_names_DEonly <- Events_names_notKissDE[Events_names_notKissDE %in% Events_names_DE]
Events_DEonly <- length(Events_names_DEonly)
Isoform_names_DEonly <- Isoform_names_notKissDE[Isoform_names_notKissDE %in% Isoform_names_DE]
Isoform_DEonly <- length(Isoform_names_DEonly)

# No regulation
Events_names_noreg <- Events_names_notKissDE[! Events_names_notKissDE %in% Events_names_DE]
Events_noreg <- length(Events_names_noreg)
Isoform_names_noreg <- Isoform_names_notKissDE[! Isoform_names_notKissDE %in% Isoform_names_DE]
Isoform_noreg <- length(Isoform_names_noreg)

# 2 Do flowchart
DiagrammeR::grViz("digraph {

graph [layout = dot, rankdir = LR]

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = Linen]

KisSplice [label = 'KisSplice \n n = 13890 Isoforms \n m = 6945 AS Events', shape = folder, fillcolor = Beige]

KissDE [label = 'KissDE \n n = 5105 Isoforms \n m = 2575 AS Events', shape = box3d]
NotKissDE [label = 'Not KissDE \n n = 8785 Isoforms \n m = 4370 AS Events', shape = box3d]

DEKissDE [label = 'KissDE-DE \n n = 944 Isoforms \n m = 784 AS Events', shape = note, fillcolor = Gainsboro]
KissDEonly [label = 'KissDE-only \n n = 4161 Isoforms \n m = 1791 AS Events', shape = note, fillcolor = Gainsboro]

DEonly [label = 'DE-only \n n = 624 Isoforms \n m = 506 AS Events', shape = note, fillcolor = Gainsboro]
noreg [label = 'No regulation \n n = 8161 Isoforms \n m = 3864 AS Events', shape = note, fillcolor = Gainsboro]

# edge definitions with the node IDs
KisSplice -> {KissDE NotKissDE} 

KissDE -> DEKissDE [color = Sienna, penwidth = 3]
KissDE -> KissDEonly [color = coral, penwidth = 1.5]

NotKissDE -> DEonly [color = coral, penwidth = 1.5]
NotKissDE -> noreg [color = gray, penwidth = 1]
}")

rm(list = ls())
gc()


###### 10. FunctionalChaos between DV,DE and DE-DV ######
# data from normCounts Heatmaps -- DE-only, DV-only, DVDE

FunctionalDischo <- list(dd = Heatmap_matrix_dd$norm, dv = Heatmap_matrix_dv$KissDE_diff_only, de = Heatmap_matrix_de$DESeq2_diff_only)
missbinde <- rownames(FunctionalDischo$dd)[-which(rownames(FunctionalDischo$dd) %in% rownames(FunctionalDischo$de))]

df <- data.frame(C = c(0,0,0,0), T1 = c(0,0,0,0), T3 = c(0,0,0,0), row.names = missbinde)
FunctionalDischo$de <- rbind(FunctionalDischo$de, df)

yep <- list()
plots <- list()
for(i in rownames(FunctionalDischo$dd)){
  dat <- as.factor(rep(c("C", "T1", "T3"),3))
  val <- c(as.numeric(FunctionalDischo$dd[i,]),as.numeric(FunctionalDischo$de[i,]), as.numeric(FunctionalDischo$dv[i,]))
  gro <- c(rep("DD",3), rep("DE",3), rep("DV", 3))
  df <- data.frame(variable = dat, value = val,  ID = gro, groups = rep(i, 9))
  yep[[i]] <- df
  p<-ggplot(df, aes(x=variable, y=value, group=gro)) +
    geom_line(aes(color=gro, size=.05))+
    geom_point(aes(color=gro, shape=gro, size=0.1))
  p <- p + scale_color_brewer(palette="Reds")+
    theme_classic() + ylim(0,1) + ylab("") + xlab("") +  theme(legend.position = "none")
  plots[[i]] <- p
}

ggarrange(plots$PS, plots$`redox homeostasis`, plots$`phytormone action`, plots$`chromatin organisation`, plots$`cell cycle organisation`,
          plots$`DNA damage response`, plots$`RNA biosynthesis`, plots$`RNA processing`, plots$`protein biosynthesis`, plots$`protein modification`,
          plots$`protein homeostasis`, plots$`cellular respiration`, plots$`cytoskeleton organisation`, plots$`cell wall organisation`, plots$`vesicle trafficking`,
          plots$`protein translocation`, plots$`solute transport`, plots$`nutrient uptake`, plots$`external stimuli response`, plots$`multi-process regulation`,
          plots$`not assigned`, plots$`enzyme classification`, plots$`CHO metabolism`, plots$`AA metabolism`, plots$`lipid metabolism`,
          plots$`nucleotide metabolism`, plots$`coenzyme metabolism`, plots$`polyamine metabolism`, plots$`secondary metabolism`,
          labels = names(yep),
          font.label = list(size = 5, color = "black", face =
                              "bold", family = NULL),
          label.x = 0.15,
          ncol = 5, nrow = 6, common.legend = T)

###### 11. Session info ######

writeLines(capture.output(sessionInfo()), "D:\\AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Logs_SessionInfo/DownStreamAnalysis_sessionInfo.txt")
