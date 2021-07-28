########################################################################################
############################ 1. KISSPLICE ##############################################
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

###### 2. Heatmap ######
# I will test this with several transformations: raw counts (not advisable), scaled to sample size (not advisable), tpm (it is okay
# but gene length normalization in this bio context dont make sense), DESeq2 normalization (kinda like it)

#BINCODES matrix:  col1:"Kiss" col2:"bincodes" | From 2.Heatmap.xlsx
bincodes<-read.delim("clipboard",header=TRUE)
bincodes$KissID <- gsub(pattern = "|", replacement = ".", x = bincodes$KissID, fixed = T)
#bindesc matrix structure: rownames=bincode num, col1: description
bindesc<-read.delim("clipboard",header=TRUE)
bindesc$code <- rownames(bindesc)
# one of the previously described
matrix <- read.delim("clipboard", header = T)

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
pdf('normDESeq22_collapsed_heatmap.pdf', width = 10, height = 10)
hp<-pheatmap(data, color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100),
             clustering_distance_rows=dist(data, method = "manhattan"),
             clustering_distance_cols=dist(t(data), method = "manhattan"),
             clustering_method="ward.D",display_numbers=TRUE)

dev.off()

###### 3. PCA ######
# normDESeq2:collapsed & raw
matrix <- read.delim("clipboard", header = T)
matrix.pca <- matrix[,-1]

raw <- read.delim("clipboard", header = T, row.names = 1)
raw <- raw[,(ncol(raw)-11):ncol(raw)]
raw_collapsed <- data.frame(C1 = raw$Control_counts1+raw$Control_counts2, 
                            C2 = raw$Control_counts3+raw$Control_counts4,
                            T1.1 = raw$T1_counts5+raw$T1_counts6,
                            T1.2 = raw$T1_counts7+raw$T1_counts8,
                            T3.1 = raw$T3_counts9+raw$T3_counts10,
                            T3.2 = raw$T3_counts11+raw$T3_counts12)

Stress <- c(rep("Control", 2), rep("Stress", 4))
Treatment <- c(rep("C",2), rep("T1", 2), rep("T3",2))
matrix.pca.done <- prcomp(t(raw_collapsed), center = F, scale. = F, rank. = 5)
fviz_eig(matrix.pca.done)

labelitass <- list(Stress, Treatment)
for(i in 1:length(labelitass)){
  data <- as.data.frame(matrix.pca.done$x)
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

###### 4. Venn & Upset presence/absence ######

matrixC <- matrix[,c(2,3)]
matrixT1 <-  matrix[,c(4,5)]
matrixT3 <-  matrix[,c(6,7)]

C <- matrix$KissID[which(rowSums(matrixC) > 0)]
T1 <- matrix$KissID[which(rowSums(matrixT1) > 0)]
T3 <- matrix$KissID[which(rowSums(matrixT3) > 0)]
data <- list(Control = C, T1 = T1, T3 = T3)
myV <- nVennR::plotVenn(data, nCycles = 20000, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile =  "nVennR.svg")

pdf("UpSet_k51.pdf")    
print(upset(fromList(data), sets = names(data), order.by = "freq", 
            keep.order = TRUE, empty.intersections = T, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3"))
dev.off()

###### 5. Kmeans ######
## normDESeq2 and raw

Kmeans_input0 <- list(raw = t(raw_collapsed), norm = t(matrix.pca))
vectortratamientos <- as.factor(c(rep("C", 2), rep("T1", 2), rep("T3", 2)))
initialcolumn <- 1
initialrow <- 1

data <- lapply(Kmeans_input0, function(x){
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

write.table(RESULTS.data$raw$kmeans_list$`16 Clusters`, file = "kmeans_all_raw_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$raw$kmeans_list$`20 Clusters`, file = "kmeans_all_raw_20.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$norm$kmeans_list$`16 Clusters`, file = "kmeans_all_norm_16.txt", sep = "\t", row.names = F, col.names = T)
write.table(RESULTS.data$norm$kmeans_list$`20 Clusters`, file = "kmeans_all_norm_20.txt", sep = "\t", row.names = F, col.names = T)

## % of events with isoforms in different kmeans clusters for each final kmeans

inputs <- list(raw_16 = RESULTS.data$raw$kmeans_list$`16 Clusters`, 
               raw_20 = RESULTS.data$raw$kmeans_list$`20 Clusters`,
               norm_16 = RESULTS.data$norm$kmeans_list$`16 Clusters`,
               norm_20 = RESULTS.data$norm$kmeans_list$`20 Clusters`)
            
resssperonoesvaca <- lapply(inputs, function(x){
  events <- unique(gsub(".2", "", gsub(pattern = ".1", replacement = "", x = x$ID, fixed = T), fixed = T))
  jj <- 0
  ss <- 0
  for(i in events){
    if(length(x$groups[which(x$ID == paste0(i, ".2"))]) == 0 | length(x$groups[which(x$ID == paste0(i, ".1"))]) == 0){
      ss <- ss + 1
      next
    }
    if(x$groups[which(x$ID == paste0(i, ".1"))] != x$groups[which(x$ID == paste0(i, ".2"))]){
      jj <- jj + 1
    }else if(x$groups[which(x$ID == paste0(i, ".1"))] == x$groups[which(x$ID == paste0(i, ".2"))]){
      jj <- jj + 0
    }
  }
  res <- jj/(length(events)-ss) * 100
  res
})

# raw_16  [1] 85.28438 | raw_20 [1] 86.98344 | norm_16 [1] 86.2635 | norm_20 [1] 88.30814
