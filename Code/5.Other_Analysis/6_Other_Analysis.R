################################################################################
############################# 5.  Other Analysis ###############################
################################################################################

library(plotly)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(compute.es)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(seqRFLP)
library(stringr)
library(uwot)

## Plot intron-length theory (first) # retry intronIC without i

arth <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Arabidopsis_thaliana/arth_chromsizes.txt")
amtr <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Amborella_trichopoda/amtr_chromsizes.txt")
zema <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Zea_mays/zema_chromsizes.txt")
orsa <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Oryza_sativa/orsa_chromsizes.txt")

gibi <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Gibi/gibi_chromsizes.txt")
gnmo <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Gnmo/gnmo_chromsizes.txt")
pita <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Pita/pita_chromsizes.txt")
pila <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Pila/pila_chromsizes.txt")
psme <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Psme/psme_chromsizes.txt")
abal <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Abal/abal_chromsizes.txt")
piab <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Piab/piab_chromsizes.txt")

# ------


arth_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Arabidopsis_thaliana/arabidopsis_thaliana.meta.iic", header = F)
amtr_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Amborella_trichopoda/amborella_trichopoda.meta.iic")
zema_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Zea_mays/zea_mays.meta.iic")
orsa_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Angiosperms/Oryza_sativa/oryza_sativa.meta.iic")

gibi_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Gibi/ginkgo_biloba.meta.iic")
gnmo_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Gnmo/gnetum_montanum.meta.iic")
pita_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Pita/pinus_taeda.meta.iic")
pila_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Pila/pinus_lambertiana.meta.iic")
psme_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Psme/pseudotsuga_menziesii.meta.iic")
abal_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Abal/abies_alba.meta.iic")
piab_introns <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/send_introns/send/Gymnosperms/Piab/picea_abies.meta.iic")



# 1) Bubble plot of genome sizes

data <- data.frame(Species = c("Pita", "Pila", "Psme", "Gnmo", "Gibi", "Amtr", "Zema", "Orsa", "Arth"),
                   GenomeSize = c(sum(pita$V2), sum(pila$V2), sum(psme$V2), sum(gnmo$V2), sum(gibi$V2), sum(amtr$V2), sum(zema$V2), sum(orsa$V2), sum(arth$V2)))
data$GenomeSize <- data$GenomeSize/max(data$GenomeSize) * 100
fig2 <- plot_ly(data, x = ~GenomeSize, y = ~Species, text = ~Species, type = 'scatter', mode = 'markers',
               marker = list(size = ~GenomeSize, opacity = 0.5, color = 'gray'))
fig2 <- fig2 %>% layout(xaxis = list(showgrid = FALSE),
                      yaxis = list(showgrid = FALSE, autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila","Psme","Gnmo", "Gibi","Amtr", "Zema", "Orsa", "Arth"), showline = F))

fig2
plotly::export(fig2, file = "IntronComparaGenomeSize_bubble.pdf")

# 2) Boxplots all introns (U2 majority)

data2 <- data.frame("Size" = c(arth_introns$V6, amtr_introns$V6, zema_introns$V6, orsa_introns$V6,
                               gibi_introns$V6, gnmo_introns$V6, pita_introns$V6, pila_introns$V6, psme_introns$V6, abal_introns$V6, piab_introns$V6),
                    "Species" = c(rep("Arth", nrow(arth_introns)), rep("Amtr", nrow(amtr_introns)), rep("Zema", nrow(zema_introns)), rep("Orsa", nrow(orsa_introns)),
                                  rep("Gibi", nrow(gibi_introns)), rep("Gnmo", nrow(gnmo_introns)), rep("Pita", nrow(pita_introns)), rep("Pila", nrow(pila_introns)),
                                  rep("Psme", nrow(psme_introns)), rep("Abal", nrow(abal_introns)), rep("Piab", nrow(piab_introns))),
                    "Clade" = c(rep("Angiosperms", nrow(arth_introns)+nrow(amtr_introns)+nrow(zema_introns)+nrow(orsa_introns)),
                                rep("Gymnosperms",  nrow(gibi_introns)+nrow(gnmo_introns)+nrow(pita_introns)+nrow(pila_introns)+nrow(psme_introns)+nrow(abal_introns)+nrow(piab_introns)))
                    )


data3 <- list()
for(i in levels(data2$Species)){
  species.filter <- data2[which(data2$Species == i),]
  species.filter <- species.filter[order(species.filter$Size, decreasing = T),]
  Thres <- sum(species.filter$Size)*0.9
  N <- 0
  for(j in 1:nrow(species.filter)){
    N <- N + species.filter$Size[j]
    if(N < Thres){ 
      next
    }else if(N >= Thres){
        N90 <- species.filter$Size[j]
      }
  }
  data3[[i]] <- species.filter[which(species.filter$Size > N90),]
}

data4 <- data2
Thres <- sum(data4$Size)*0.9
N <- 0
for(i in 1:nrow(data4)){
  N <- N + data4$Size[i]
  if(N < Thres){ 
    next
  }else if(N >= Thres){
    N90 <- data4$Size[i]
  }
}
data4 <- data4[which(data4$Size > N90),]

fig <- plot_ly(orientation = "h")
fig <- fig %>% add_boxplot(x=log10(data3$Arth$Size), y=data3$Arth$Species, name='Arth', boxpoints = T, fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Amtr$Size), data3$Amtr$Species, name='Amtr', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Zema$Size), data3$Zema$Species, name='Zema', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Orsa$Size), data3$Orsa$Species, name='Orsa', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Gibi$Size), data3$Gibi$Species, name='Gibi', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Gnmo$Size), data3$Gnmo$Species, name='Gnmo', boxpoints = T, fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Pita$Size), data3$Pita$Species, name='Pita', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Pila$Size), data3$Pila$Species, name='Pila', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Psme$Size), data3$Psme$Species, name='Psme', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Abal$Size), data3$Abal$Species, name='Abal', boxpoints = T, fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% add_boxplot(x=log10(data3$Piab$Size), data3$Piab$Species, name='Piab', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig <- fig %>% layout(
  title = 'Intron length (log10)',
  yaxis = list( autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila", "Piab",  "Psme",  "Abal", "Gnmo", "Gibi","Orsa", "Zema", "Amtr", "Arth"), showline = F),
  margin = list(r = 10, t = 25, b = 40, l = 110),
  showlegend = F
)
fig

data4 <- data4[-which(data4$Species == "Abal"),]
data4 <- data4[-which(data4$Species == "Piab"),]
fig_0 <- plot_ly(orientation = "h")
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Arth")]), data4$Species[which(data4$Species == "Arth")], name='Arth', boxpoints = F, fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Amtr")]), data4$Species[which(data4$Species == "Amtr")], name='Amtr', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Zema")]), data4$Species[which(data4$Species == "Zema")], name='Zema', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Orsa")]), data4$Species[which(data4$Species == "Orsa")], name='Orsa', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Gibi")]), data4$Species[which(data4$Species == "Gibi")], name='Gibi', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Gnmo")]), data4$Species[which(data4$Species == "Gnmo")], name='Gnmo', boxpoints = F, fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Pita")]), data4$Species[which(data4$Species == "Pita")], name='Pita', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Pila")]), data4$Species[which(data4$Species == "Pila")], name='Pila', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% add_boxplot(x=log10(data4$Size[which(data4$Species == "Psme")]), data4$Species[which(data4$Species == "Psme")], name='Psme', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig_0 <- fig_0 %>% layout(
  yaxis = list( autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila", "Psme", "Gnmo", "Gibi","Amtr", "Zema", "Orsa", "Arth"), showline = F),
  margin = list(r = 10, t = 25, b = 40, l = 110),
  showlegend = F
)
fig_0
plotly::export(fig, file = "IntronCompara_intronICwithI_boxplot_F.pdf")

# 3) Boxplot U12 introns

arth_introns_U12 <- arth_introns[which(as.numeric(as.character(arth_introns$V2)) > 0),]
amtr_introns_U12 <- amtr_introns[which(as.numeric(as.character(amtr_introns$V2)) > 0),]
zema_introns_U12 <- zema_introns[which(as.numeric(as.character(zema_introns$V2)) > 0),]
orsa_introns_U12 <- orsa_introns[which(as.numeric(as.character(orsa_introns$V2)) > 0),]
gibi_introns_U12 <- gibi_introns[which(as.numeric(as.character(gibi_introns$V2)) > 0),]
gnmo_introns_U12 <- gnmo_introns[which(as.numeric(as.character(gnmo_introns$V2)) > 0),]
psme_introns_U12 <- psme_introns[which(as.numeric(as.character(psme_introns$V2)) > 0),]
pila_introns_U12 <- pila_introns[which(as.numeric(as.character(pila_introns$V2)) > 0),]
pita_introns_U12 <- pita_introns[which(as.numeric(as.character(pita_introns$V2)) > 0),]
piab_introns_U12 <- piab_introns[which(as.numeric(as.character(piab_introns$V2)) > 0),]
abal_introns_U12 <- abal_introns[which(as.numeric(as.character(abal_introns$V2)) > 0),]



data5 <- data.frame("Size" = c(arth_introns_U12$V6, amtr_introns_U12$V6, zema_introns_U12$V6, orsa_introns_U12$V6,
                               gibi_introns_U12$V6, gnmo_introns_U12$V6, pita_introns_U12$V6, pila_introns_U12$V6, psme_introns_U12$V6, abal_introns_U12$V6, piab_introns_U12$V6),
                    "Species" = c(rep("Arth", nrow(arth_introns_U12)), rep("Amtr", nrow(amtr_introns_U12)), rep("Zema", nrow(zema_introns_U12)), rep("Orsa", nrow(orsa_introns_U12)),
                                  rep("Gibi", nrow(gibi_introns_U12)), rep("Gnmo", nrow(gnmo_introns_U12)), rep("Pita", nrow(pita_introns_U12)), rep("Pila", nrow(pila_introns_U12)),
                                  rep("Psme", nrow(psme_introns_U12)), rep("Abal", nrow(abal_introns_U12)), rep("Piab", nrow(piab_introns_U12))),
                    "Clade" = c(rep("Angiosperms", nrow(arth_introns_U12)+nrow(amtr_introns_U12)+nrow(zema_introns_U12)+nrow(orsa_introns_U12)),
                                rep("Gymnosperms",  nrow(gibi_introns_U12)+nrow(gnmo_introns_U12)+nrow(pita_introns_U12)+nrow(pila_introns_U12)+nrow(psme_introns_U12)+nrow(abal_introns_U12)+nrow(piab_introns_U12)))
)


fig1 <- plot_ly(orientation = "h")
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Arth")]), data5$Species[which(data5$Species == "Arth")], name='Arth', boxpoints = T, fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Amtr")]), data5$Species[which(data5$Species == "Amtr")], name='Amtr', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Zema")]), data5$Species[which(data5$Species == "Zema")], name='Zema', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Orsa")]), data5$Species[which(data5$Species == "Orsa")], name='Orsa', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Gibi")]), data5$Species[which(data5$Species == "Gibi")], name='Gibi', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Gnmo")]), data5$Species[which(data5$Species == "Gnmo")], name='Gnmo', boxpoints = T, fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Pita")]), data5$Species[which(data5$Species == "Pita")], name='Pita', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Pila")]), data5$Species[which(data5$Species == "Pila")], name='Pila', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Psme")]), data5$Species[which(data5$Species == "Psme")], name='Psme', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Piab")]), data5$Species[which(data5$Species == "Piab")], name='Piab', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))
fig1 <- fig1 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Abal")]), data5$Species[which(data5$Species == "Abal")], name='Abal', boxpoints = T,  fillcolor = "white", line = list(color="black", width=2))

fig1 <- fig1 %>% layout(
  yaxis = list( autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila", "Piab","Psme","Abal", "Gnmo", "Gibi","Amtr", "Zema", "Orsa", "Arth"), showline = F),
  margin = list(r = 10, t = 25, b = 40, l = 110),
  showlegend = F
)
fig1


data6 <- data5
Thres <- sum(data6$Size)*0.9
N <- 0
for(i in 1:nrow(data6)){
  N <- N + data6$Size[i]
  if(N < Thres){ 
    next
  }else if(N >= Thres){
    N90 <- data6$Size[i]
  }
}
data6 <- data6[which(data6$Size > N90),]

fig3 <- plot_ly(orientation = "h")
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Arth")]), data5$Species[which(data5$Species == "Arth")], name='Arth', boxpoints = F, fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Amtr")]), data5$Species[which(data5$Species == "Amtr")], name='Amtr', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Zema")]), data5$Species[which(data5$Species == "Zema")], name='Zema', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Orsa")]), data5$Species[which(data5$Species == "Orsa")], name='Orsa', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Gibi")]), data5$Species[which(data5$Species == "Gibi")], name='Gibi', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Gnmo")]), data5$Species[which(data5$Species == "Gnmo")], name='Gnmo', boxpoints = F, fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Pita")]), data5$Species[which(data5$Species == "Pita")], name='Pita', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Pila")]), data5$Species[which(data5$Species == "Pila")], name='Pila', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))
fig3 <- fig3 %>% add_boxplot(x=log10(data5$Size[which(data5$Species == "Psme")]), data5$Species[which(data5$Species == "Psme")], name='Psme', boxpoints = F,  fillcolor = "white", line = list(color="black", width=2))

fig3 <- fig3 %>% layout(
  yaxis = list( autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila","Psme","Gnmo", "Gibi","Amtr", "Zema", "Orsa", "Arth"), showline = F),
  margin = list(r = 10, t = 25, b = 40, l = 110),
  showlegend = F
)
fig3

# 5) Barplot U12 dinucleotides

pal_annotation<-c("#104E8B","#7A378B","#8B3A3A","#778899","#548B54","#CD853F","#8B5F65","#CD9B1D", "#54C7FC", "#F0E68C", "#FF83FA", "#FFB6C1")

U12_flanks <- list(arth = arth_introns_U12, orsa = orsa_introns_U12, zema = zema_introns_U12, amtr = amtr_introns_U12,
                   gibi = gibi_introns_U12, gnmo = gnmo_introns_U12, psme = psme_introns_U12, pila = pila_introns_U12,
                   pita = pita_introns_U12)


U12_flanks_percents <- lapply(U12_flanks, function(x){
  percents <- list()
  total <- nrow(x)
  for(j in levels(as.factor(as.character(x[,3])))){
    percentage <- length(which(x[,3] == j))/total * 100
    percents[j] <- percentage
  }
  if(length(percents) > 3){
    percents <- c(sort(unlist(percents), decreasing = T)[1:3], Others = 100-(sum(sort(unlist(percents), decreasing = T)[1:3])))
  }else if(length(percents) <= 3){
    percents["Others"] <- 0
  }
  data.frame(Percents = sort(unlist(percents), decreasing = T))
})

U12_barplot = data.frame(Species = c("Arth", "Orsa", "Zema", "Amtr", "Gibi", "Gnmo", "Psme", "Pila", "Pita"),
                         'GT-AG' = c(72.81, 71.5, 69.2, 100, 73.9, 82.6, 80.5, 78.85, 79.25),
                         'AT-AC' = c(25.89, 22.4, 15.3, 0, 16.4, 16.9, 16.8, 18.3, 19.1),
                         'AT-AA' = c(1.3,1.7,0,0,3.7,0,1.7,0.85,0.85),
                         'GT-AT' = c(0,0,2.1,0,0,0,0,0,0),
                         'GC-AG' = c(0,0,0,0,0,0.5,0,0,0),
                         Others = c(0,4.4,13.4,0,6,0,1,2,0.8))

fig4 <- plot_ly(U12_barplot, x = ~GT.AG, y = ~Species, type = 'bar', orientation = 'h', name='GT-AG',
               marker = list(color = pal_annotation[1],
                             line = list(color = "black",
                                         width = 1)))
fig4 <- fig4 %>% add_trace(x = ~AT.AC, name = 'AT-AC',
                         marker = list(color = pal_annotation[6],
                                       line = list(color = "black",
                                                   width = 1)))
fig4 <- fig4 %>% add_trace(x = ~AT.AA, name = 'AT-AA',
                         marker = list(color = pal_annotation[3],
                                       line = list(color = "black",
                                                   width = 1)))
fig4 <- fig4 %>% add_trace(x = ~GT.AT, name = 'GT-AT',
                         marker = list(color = pal_annotation[4],
                                       line = list(color = "black",
                                                   width = 1)))
fig4 <- fig4 %>% add_trace(x = ~GC.AG, name = 'GC-AG',
                         marker = list(color = pal_annotation[5],
                                       line = list(color = "black",
                                                   width = 1)))
fig4 <- fig4 %>% add_trace(x = ~Others, name = 'Others',
                         marker = list(color = pal_annotation[2],
                                       line = list(color = "black",
                                                   width = 1)))
fig4 <- fig4 %>% layout(barmode = 'stack',
                      xaxis = list(title = ""),
                      yaxis = list( autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila","Psme","Gnmo", "Gibi","Amtr", "Zema", "Orsa", "Arth"), showline = F))

fig4

# 6) Barplot all dinucleotides (U2 majority)

intron_flanks <- list(arth = arth_introns, orsa = orsa_introns, zema = zema_introns, amtr = amtr_introns,
                   gibi = gibi_introns, gnmo = gnmo_introns, psme = psme_introns, pila = pila_introns,
                   pita = pita_introns)

intron_flanks_percents <- lapply(intron_flanks, function(x){
  percents <- list()
  total <- nrow(x)
  for(j in levels(as.factor(as.character(x[,3])))){
    percentage <- length(which(x[,3] == j))/total * 100
    percents[j] <- percentage
  }
  if(length(percents) > 3){
    percents <- c(sort(unlist(percents), decreasing = T)[1:3], Others = 100-(sum(sort(unlist(percents), decreasing = T)[1:3])))
  }else if(length(percents) <= 3){
    percents["Others"] <- 0
  }
  data.frame(Percents = sort(unlist(percents), decreasing = T))
})

introns_barplot = data.frame(Species = c("Arth", "Orsa", "Zema", "Amtr", "Gibi", "Gnmo", "Psme", "Pila", "Pita"),
                         'GT-AG' = c(98.36, 93.8, 86, 98.4, 90, 98.4, 98.9, 98.01, 98.7),
                         'GC-AG' = c(1.1,0.32,1.4,1.58,0.8,1.5,0.9,1.12,1),
                         'AT-AC' = c(0.12, 0, 0, 0, 0, 0.05, 0.06, 0, 0.05),
                         'CT-AC' = c(0,1.08,0,0,4.8,0,0,0,0),
                         'AA-CA' = c(0,0,0,0.01,0,0,0,0,0),
                         'NN-AG' = c(0,0,0,0,0,0,0,0.11,0),
                         'GT-TG' = c(0,0.8,0.9,0,0,0,0,0,0),
                         Others = c(0.42,4,11.7,0.01,4.4,0.05,0.14,0.76,0.25))

fig5 <- plot_ly(introns_barplot, x = ~GT.AG, y = ~Species, type = 'bar', orientation = 'h', name='GT-AG',
                marker = list(color = pal_annotation[1],
                              line = list(color = "black",
                                          width = 1)))
fig5 <- fig5 %>% add_trace(x = ~GC.AG, name = 'GC-AG',
                           marker = list(color = pal_annotation[5],
                                         line = list(color = "black",
                                                     width = 1)))
fig5 <- fig5 %>% add_trace(x = ~AT.AC, name = 'AT-AC',
                           marker = list(color = pal_annotation[6],
                                         line = list(color = "black",
                                                     width = 1)))
fig5 <- fig5 %>% add_trace(x = ~CT.AC, name = 'CT-AC',
                           marker = list(color = pal_annotation[10],
                                         line = list(color = "black",
                                                     width = 1)))
fig5 <- fig5 %>% add_trace(x = ~AA.CA, name = 'AA-CA',
                           marker = list(color = pal_annotation[9],
                                         line = list(color = "black",
                                                     width = 1)))
fig5 <- fig5 %>% add_trace(x = ~NN.AG, name = 'NN-AG',
                           marker = list(color = pal_annotation[8],
                                         line = list(color = "black",
                                                     width = 1)))
fig5 <- fig5 %>% add_trace(x = ~GT.TG, name = 'GT-TG',
                           marker = list(color = pal_annotation[7],
                                         line = list(color = "black",
                                                     width = 1)))
fig5 <- fig5 %>% add_trace(x = ~Others, name = 'Others',
                           marker = list(color = pal_annotation[2],
                                         line = list(color = "black",
                                                     width = 1)))
fig5 <- fig5 %>% layout(barmode = 'stack',
                        xaxis = list(title = ""),
                        yaxis = list( autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila","Psme","Gnmo", "Gibi","Amtr", "Zema", "Orsa", "Arth"), showline = F))

fig5

# 7) Line plot effect sizes 

Arth_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$arth$V6[-which(intron_flanks$arth$V2 == ".")]), alternative = "two.sided")
Arth_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$arth$V6), alternative = "two.sided")

Orsa_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$orsa$V6[-which(intron_flanks$orsa$V2 == ".")]), alternative = "two.sided")
Orsa_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$orsa$V6), alternative = "two.sided")

Zema_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$zema$V6[-which(intron_flanks$zema$V2 == ".")]), alternative = "two.sided")
Zema_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$zema$V6), alternative = "two.sided")

Amtr_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$amtr$V6[-which(intron_flanks$amtr$V2 == ".")]), alternative = "two.sided")
Amtr_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$amtr$V6), alternative = "two.sided")

Gibi_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$gibi$V6[-which(intron_flanks$gibi$V2 == ".")]), alternative = "two.sided")
Gibi_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$gibi$V6), alternative = "two.sided")

Gnmo_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$gnmo$V6[-which(intron_flanks$gnmo$V2 == ".")]), alternative = "two.sided")
Gnmo_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$gnmo$V6), alternative = "two.sided")

Psme_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$psme$V6[-which(intron_flanks$psme$V2 == ".")]), alternative = "two.sided")
Psme_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$psme$V6), alternative = "two.sided")

Pila_All <- t.test(as.numeric(intron_flanks$pita$V6[-which(intron_flanks$pita$V2 == ".")]), as.numeric(intron_flanks$pila$V6[-which(intron_flanks$pila$V2 == ".")]), alternative = "two.sided")
Pila_U12 <- t.test(as.numeric(U12_flanks$pita$V6), as.numeric(U12_flanks$pila$V6), alternative = "two.sided")

arth_all_tes <- tes(t=Arth_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$arth[-which(intron_flanks$arth$V2 == ".",)])))
arth_u12_tes <- tes(t=Arth_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$arth))

orsa_all_tes <- tes(t=Orsa_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$orsa[-which(intron_flanks$orsa$V2 == ".",)])))
orsa_u12_tes <- tes(t=Orsa_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$orsa))

zema_all_tes <- tes(t=Zema_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$zema[-which(intron_flanks$zema$V2 == ".",)])))
zema_u12_tes <- tes(t=Zema_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$zema))

amtr_all_tes <- tes(t=Amtr_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$amtr[-which(intron_flanks$amtr$V2 == ".",)])))
amtr_u12_tes <- tes(t=Amtr_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$amtr))

gibi_all_tes <- tes(t=Gibi_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$gibi[-which(intron_flanks$gibi$V2 == ".",)])))
gibi_u12_tes <- tes(t=Gibi_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$gibi))

gnmo_all_tes <- tes(t=Gnmo_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$gnmo[-which(intron_flanks$gnmo$V2 == ".",)])))
gnmo_u12_tes <- tes(t=Gnmo_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$gnmo))

psme_all_tes <- tes(t=Psme_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$psme[-which(intron_flanks$psme$V2 == ".",)])))
psme_u12_tes <- tes(t=Psme_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$psme))

pila_all_tes <- tes(t=Pila_All$statistic, as.numeric(nrow(intron_flanks$pita[-which(intron_flanks$pita$V2 == "."),])), as.numeric(nrow(intron_flanks$pila[-which(intron_flanks$pila$V2 == ".",)])))
pila_u12_tes <- tes(t=Pila_U12$statistic, nrow(U12_flanks$pita),nrow(U12_flanks$pila))

EffectSize <- data.frame(Species = c("Arth", "Orsa", "Zema", "Amtr", "Gibi", "Gnmo", "Psme", "Pila", "Pita"),
           All_d = abs(c(arth_all_tes$d, orsa_all_tes$d, zema_all_tes$d, amtr_all_tes$d, gibi_all_tes$d, gnmo_all_tes$d, psme_all_tes$d, pila_all_tes$d, 0)),
           U12_d = abs(c(arth_u12_tes$d, orsa_u12_tes$d, zema_u12_tes$d, amtr_u12_tes$d, gibi_u12_tes$d, gnmo_u12_tes$d, psme_u12_tes$d, pila_u12_tes$d, 0)))

fig6 <- plot_ly(EffectSize, y = ~Species, orientation = "h") %>% layout(yaxis = list( autorange = T, categoryorder = "array", categoryarray = c("Pita", "Pila","Psme","Gnmo", "Gibi","Amtr", "Zema", "Orsa", "Arth"), showline = F))
fig6 <- fig6 %>% add_trace(x = EffectSize$All_d, y= EffectSize$Species,  type = "scatter", mode = "line+markers", line = list(color = "darkgoldenrod", dash = "dot", width = 4), name = "All introns")
fig6 <- fig6 %>% add_trace(x = EffectSize$U12_d, y= EffectSize$Species,  type = "scatter", mode = "line+markers", line = list(color = "purple", dash = "dot", width = 4), name = "U12 introns")
fig6

# 8) Intron compara final fig

fig_final <- subplot(fig2, fig_0, fig3, fig5, fig4, fig6)
fig_final %>% layout(title = "Intron comparative analysis")
plotly::orca(fig_final %>% layout(title = "Intron comparative analysis"), file = "IntronCompara_finalfig5.pdf", format = "pdf", width = 2000, height = 900)

## Barplot expression switch 40ºC and large histo DS vs DE (fisher text) 

# Barplot

Candidates <- read.delim("clipboard", header = T) # /Data/5.Validation/Bio
Events.data <- read.delim("clipboard", header = T) # /Data/2.KissDE
pal_bars <- c("lightgoldenrod", "lightsteelblue3")

data.list <- list()
plot.list <- list()
for(i in levels(Candidates$KissID)){
  custom <- Events.data[which(Events.data$X1.ID == i),]
  total.C <- sum(custom[,5:8], custom[,17:20])
  total.T1 <- sum(custom[,9:12], custom[,21:24])
  total.T3 <- sum(custom[,13:16], custom[,25:28])
  custom.df <- data.frame(Treatment = rep(c("C", "T1", "T3"),2),
                          Variants  = c(rep("Variant1",3), rep("Variant2", 3)),
                          value     = c(sum(custom[,5:8])/total.C, sum(custom[,9:12])/total.T1, sum(custom[,13:16])/total.T3,
                                        sum(custom[,17:20])/total.C, sum(custom[,21:24])/total.T1, sum(custom[,25:28])/total.T3),
                          kissID = c(rep(i, 6)))
  anno<-ggplot(custom.df,aes(x=Treatment,y=value,fill=Variants,color=variants,width=0.75),size=0.01)+
    geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Ratio")+
    scale_fill_manual(values=rev(pal_bars))+
    ggtitle(Candidates$v2Anno[which(Candidates$KissID == i)]) + theme(plot.title = element_text(size = 5))
  data.list[[i]] <- custom.df
  plot.list[[i]] <- anno
  }
pdf("IsformsSwitch_2.pdf")
grid.arrange(plot.list$`bcc_12324|Cycle_0`, plot.list$`bcc_15225|Cycle_0`,
             plot.list$`bcc_17045|Cycle_0`, plot.list$`bcc_29432|Cycle_2`,
             plot.list$`bcc_3998|Cycle_0`, plot.list$`bcc_40022|Cycle_1`,
             plot.list$`bcc_9268|Cycle_9`, plot.list$`bcc_9458|Cycle_1`,
             nrow = 4, ncol = 2)
dev.off()
data.list <- do.call(rbind, data.list)
write.table(data.list, file = "IsoformsSwitch.txt", sep = "\t", col.names = T, row.names = F)

# Large Histo DS vs DE: df- value=%DS, Treatment=CT1, CT3, T1T3, DE=DEornonDE (Fisher test)

total <- c(c(paste0(KissDE_diff$CT1$`1.ID`[which(KissDE_diff$CT1$`21.lowcounts`== "FALSE")], ".1"), paste0(KissDE_diff$CT1$`1.ID`[which(KissDE_diff$CT1$`21.lowcounts`  == "FALSE")], ".2")),
                c(paste0(KissDE_diff$CT3$`1.ID`[which(KissDE_diff$CT3$`21.lowcounts`== "FALSE")], ".1"), paste0(KissDE_diff$CT3$`1.ID`[which(KissDE_diff$CT3$`21.lowcounts`  == "FALSE")], ".2")),
                c(paste0(KissDE_diff$T1T3$`1.ID`[which(KissDE_diff$T1T3$`21.lowcounts`== "FALSE")], ".1"), paste0(KissDE_diff$T1T3$`1.ID`[which(KissDE_diff$T1T3$`21.lowcounts`  == "FALSE")], ".2")))
total <- length(unique(total))

DS.CT1 <- c(paste0(KissDE_diff$CT1$`1.ID`[which(KissDE_diff$CT1$`21.lowcounts`== "FALSE")], ".1"), paste0(KissDE_diff$CT1$`1.ID`[which(KissDE_diff$CT1$`21.lowcounts`  == "FALSE")], ".2"))
NonDE.CT1 <- length(which((DS.CT1 %in% DESeq2_diff$CT1$KissID) == FALSE))/total * 100
DE.CT1 <- length(which((DS.CT1 %in% DESeq2_diff$CT1$KissID) == TRUE))/total * 100

DS.CT3 <- c(paste0(KissDE_diff$CT3$`1.ID`[which(KissDE_diff$CT3$`21.lowcounts`== "FALSE")], ".1"), paste0(KissDE_diff$CT3$`1.ID`[which(KissDE_diff$CT3$`21.lowcounts`  == "FALSE")], ".2"))
NonDE.CT3 <- length(which((DS.CT3 %in% DESeq2_diff$CT3$KissID) == FALSE))/total * 100
DE.CT3 <- length(which((DS.CT3 %in% DESeq2_diff$CT3$KissID) == TRUE))/total * 100

DS.T1T3 <- c(paste0(KissDE_diff$T1T3$`1.ID`[which(KissDE_diff$T1T3$`21.lowcounts`== "FALSE")], ".1"), paste0(KissDE_diff$T1T3$`1.ID`[which(KissDE_diff$T1T3$`21.lowcounts`  == "FALSE")], ".2"))
NonDE.T1T3 <- length(which((DS.T1T3 %in% DESeq2_diff$T1T3$KissID) == FALSE))/total * 100
DE.T1T3 <- length(which((DS.T1T3 %in% DESeq2_diff$T1T3$KissID) == TRUE))/total * 100

large.histo <- data.frame(DS= c(DE.CT1, NonDE.CT1, DE.CT3, NonDE.CT3, DE.T1T3, NonDE.T1T3),
                          Treatment = c(rep("CT1",2), rep("CT3", 2), rep("T1T3", 2)),
                          DE = c("DE", "Non-DE", "DE", "Non-DE", "DE", "Non-DE"))
write.table(large.histo, "DEvsDS.txt", row.names = F, col.names = T)

bp <- ggbarplot(
  large.histo, x = "Treatment", y = "DS",  
  color= "DE", palette = c("#00AFBB", "#E7B800"),
  position = position_dodge(0.8), size = 2
)

# add p-values with affinity

xtab <- as.table(rbind(
  c(11.89*total/100, 61.01*total/100), c(7.73*total/100, 47.46*total/100),
  c(1.36*total/100, 27.98*total/100)
))
dimnames(xtab) <- list(
  Treatment = c("CT1", "CT3", "T1T3"),
  DE = c("DE", "Non-DE")
)
xtab
# Compare the proportion of DE and Non-DE in each category
stat.test <- row_wise_fisher_test(xtab, p.adjust.method = "bonferroni")

# export and add with affinity

write.table(stat.test, file = "DEvsDS_conditions_pvaluesADDAFFINITY.txt", col.names = T, row.names = F)

## Descriptive: %PTCs, %CDS, %UTRs (third)

# From .bam for each kissID extract mapping coordinates | with CodAn GTF annotation match UTRs or CDS | in case of CDS -- ORFs sequences divide by codons and see if longest ORF have PTC associated with the event

# inputs

kissIDs <- Splicing_annotation$KissID
bam.coords <- Transsbyss.bam
CodAn.gtf <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Transabyss_Full//annotation.gtf")
ORFS.seq <- seqinr::read.fasta("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Transabyss_Full//ORF_sequences.fasta") # do with full Transabyss
fiveUTR.seq <- seqinr::read.fasta("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Transabyss_Full//5utr_sequences.fasta") # do with full Transabyss
threeUTR.seq <- seqinr::read.fasta("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Transabyss_Full//3utr_sequences.fasta") # do with full Transabyss
KissSeq <- seqinr::read.fasta("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.results_k51_type1.fa")
Transabyss <- seqinr::read.fasta("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/1.Final query/KisSplice_Transabyss.fa")
names(KissSeq) <- Splicing_annotation$KissID
VariableSeq <- seqinr::read.fasta("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/1.KisSplice/kmer51/3.allSparts.fasta")
for(i in 1:length(VariableSeq)){
  names(VariableSeq)[i] <- paste0(strsplit(names(VariableSeq)[i], "|", fixed = T)[[1]][1], sep = "|", strsplit(names(VariableSeq)[i], "|", fixed = T)[[1]][2])
}
# AS transcript structure classification

Model.results <- data.frame(Isoforms = names(KissSeq), Model = NA, PTC = NA)
for(i in seq(1, length(KissSeq))){
  
  # ID 
  
  ID <- as.character(kissIDs[i])
  
  # Bowtie2 hit Check 
  
  if(length(which(bam.coords[[1]]$qname == ID)) == 0){
      Model.results[i,2] <- "NoMapping"
      next
  }
  
  # Coherency Checks
  
  if(bam.coords[[1]]$qwidth[which(bam.coords[[1]]$qname == ID)] != length(KissSeq[[i]])){ stop ("\n Length Coherency between BAM and KissSeq error \n")}
  
  # Mapping coordinate + kmer (Variable part beginning)
  
  S.Start <- bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == ID)] + 52
  
  # Check available predictions: 
  
  if(length(which(CodAn.gtf$V1 == ID)) == 0){
    Model.results[i,2] <- "NoPrediction"
    next
  }
  
  # Extract prediction information
  
  CodAn.isoform <- CodAn.gtf[which(CodAn.gtf$V1 == ID),]
  
  if(CodAn.isoform$V7[1] == "+"){
    if(length(which(CodAn.isoform$V3 == "start_codon")) != 0){
      if(S.Start < CodAn.isoform$V4[which(CodAn.isoform$V3 == "start_codon")]){
        Model.results[i,2] <- "5UTR"
      }
      StartCodon <- CodAn.isoform[which(CodAn.isoform$V3 == "start_codon"),]
      if(S.Start %in% StartCodon$V4:StartCodon$V5){
        Model.results[i,2] <- "StartCodon_CDS"
      }
    }
    
    if(length(which(CodAn.isoform$V3 == "stop_codon")) != 0){
      if(S.Start > CodAn.isoform$V5[which(CodAn.isoform$V3 == "stop_codon")]){
        Model.results[i,2] <- "3UTR"
      }
      StopCodon <- CodAn.isoform[which(CodAn.isoform$V3 == "stop_codon"),]
      if(S.Start %in% StopCodon$V4:StopCodon$V5){
        Model.results[i,2] <- "StopCodon_CDS"
      }
    }
    
    CDS <- CodAn.isoform[which(CodAn.isoform$V3 == "CDS"),]
    if(S.Start %in% CDS$V4:CDS$V5){
      Model.results[i,2] <- "CDS"
    }

    if(is.na(Model.results[i,2])){ stop( "\n Something go wrong, conditional optional requires \n")}
  }else if(CodAn.isoform$V7[1] == "-"){
    if(length(which(CodAn.isoform$V3 == "start_codon")) != 0){
      if(S.Start > CodAn.isoform$V5[which(CodAn.isoform$V3 == "start_codon")]){
        Model.results[i,2] <- "5UTR"
      }
      StartCodon <- CodAn.isoform[which(CodAn.isoform$V3 == "start_codon"),]
      if(S.Start %in% StartCodon$V4:StartCodon$V5){
        Model.results[i,2] <- "StartCodon_CDS"
      }
    }
    
    if(length(which(CodAn.isoform$V3 == "stop_codon")) != 0){
      if(S.Start < CodAn.isoform$V4[which(CodAn.isoform$V3 == "stop_codon")]){
        Model.results[i,2] <- "3UTR"
      }
      StopCodon <- CodAn.isoform[which(CodAn.isoform$V3 == "stop_codon"),]
      if(S.Start %in% StopCodon$V4:StopCodon$V5){
        Model.results[i,2] <- "StopCodon_CDS"
      }
    }
    
    CDS <- CodAn.isoform[which(CodAn.isoform$V3 == "CDS"),]
    if(S.Start %in% CDS$V4:CDS$V5){
      Model.results[i,2] <- "CDS"
    }
    

    
    if(is.na(Model.results[i,2])){ stop( "\n Something go wrong, conditional optional requires \n")}
  } 
  # CDS-PTC Search
  
  if(Model.results[i,2] == "CDS" | Model.results[i,2] == "StopCodon_CDS" | Model.results[i,2] == "StartCodon_CDS"){
    
    # all coordinates are refered based in + strand
    
    if(CodAn.isoform$V7[1] == "+"){ # if strand + -- 5'-CDS-3'
      
      Variable.i <- VariableSeq[[which(names(VariableSeq) == stri_sub(ID, 1, nchar(ID)-2))]]
      
      # nested conditionals just in case dont have UTRS
      
      if(length(which(names(fiveUTR.seq) == ID)) == 0 & length(which(names(threeUTR.seq) == ID)) == 0){
        Transcript.mod <- paste0(paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""))
      }else if(length(which(names(fiveUTR.seq) == ID)) == 0 & length(which(names(threeUTR.seq) == ID)) != 0){
        Transcript.mod <- paste0(paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(threeUTR.seq[[ which(names(threeUTR.seq) == ID)]]), collapse = ""))
      }else if(length(which(names(fiveUTR.seq) == ID)) != 0 & length(which(names(threeUTR.seq) == ID)) == 0){
        Transcript.mod <- paste0(paste(as.character(fiveUTR.seq[[ which(names(fiveUTR.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""))
      }else if(length(which(names(fiveUTR.seq) == ID)) != 0 & length(which(names(threeUTR.seq) == ID)) != 0){
        Transcript.mod <- paste0(paste(as.character(fiveUTR.seq[[ which(names(fiveUTR.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(threeUTR.seq[[ which(names(threeUTR.seq) == ID)]]), collapse = "")
        )
      }
      
      
      if(!all(compare.DNA((strsplit(Transcript.mod, split = ""))[[1]], Transabyss[[which(names(Transabyss) == ID)]]))){
        # checking coherency of pasted transcripts parts
        stop(" \n Not coherency pasting all transcript parts \n")
      }
      
      if(stri_sub(ID, -2) == ".1"){ # 1) Search S part in upper (.1). If match search PTCs in frame ORF. If not match -- 2) info in .2? pass to .2; if not info in .2 paste S seq to .1 and search PTCs in frame
        
        S.mod <- stri_sub(Transcript.mod, S.Start-1, S.Start-1+length(Variable.i)-1)
        S.mod <- (strsplit(x = S.mod, split = ""))[[1]]
        
        if(length(S.mod) != length(Variable.i)){ stop("\n Different events variable seqs length \n")}
        
        Condition <- compare.DNA(S.mod, as.character(Variable.i))
        
        if(length(which(Condition == TRUE)) >= 0.95*length(Condition)){ # mapped transcript have this splicing event variation into account
          
          StopCodons <- c("tag", "taa", "tga")
          
          # search PTCs in frame
          
          ORF.i <- paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = "")
          
          # Because of Full Transcript Models only one option of frame
          
          bananasplit <- (strsplit(ORF.i, split = ""))[[1]]
          bananasplitxxi <- paste0(bananasplit[c(TRUE, FALSE, FALSE)], bananasplit[c(FALSE, TRUE, FALSE)], bananasplit[c(FALSE, FALSE, TRUE)])
          if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) > 1){
            Model.results[i,3] <- "PTC"
          }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 1){
            if(length(which(CodAn.isoform$V3 == "stop_codon")) == 1){
              Model.results[i,3] <- "CDS_change"  
            }else if(length(which(CodAn.isoform$V3 == "stop_codon")) == 0){
              Model.results[i,3] <- "PTC"
            }else {stop("\n Wait what \n")}
          }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 0){
            Model.results[i,3] <- "CDS_change"
          }
          
        }else if(length(which(Condition == TRUE)) < 0.95*length(Condition)){ # mapped transcript does not have this splicing events variation into account
            
            # paste S seq to ORF seq
            
            S.Start.ORF.i <- S.Start - min(CodAn.isoform$V4)
            ORF.i <- paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = "")
            S.ORF.i <- paste0(str_sub(ORF.i, start = 1, end = S.Start.ORF.i-1), paste(as.character(Variable.i), collapse = ""), str_sub(ORF.i, start = S.Start.ORF.i, end = nchar(ORF.i)))
            
            # search PTCs in frame
            
            StopCodons <- c("tag", "taa", "tga")
            
            # Because of Full Transcript Models only one option of frame
            
            bananasplit <- (strsplit(S.ORF.i, split = ""))[[1]]
            bananasplitxxi <- paste0(bananasplit[c(TRUE, FALSE, FALSE)], bananasplit[c(FALSE, TRUE, FALSE)], bananasplit[c(FALSE, FALSE, TRUE)])
            if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) > 1){
              Model.results[i,3] <- "PTC"
            }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 1){
              if(length(which(CodAn.isoform$V3 == "stop_codon")) == 1){
                Model.results[i,3] <- "CDS_change"  
              }else if(length(which(CodAn.isoform$V3 == "stop_codon")) == 0){
                Model.results[i,3] <- "PTC"
              }else {stop("\n Wait what \n")}
            }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 0){
              Model.results[i,3] <- "CDS_change"
            }
          }
          
      }else if(stri_sub(ID, -2) == ".2"){ 
        
        # paste S seq to ORF seq
        
        S.Start.ORF.i <- S.Start - min(CodAn.isoform$V4)
        ORF.i <- paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = "")
        S.ORF.i <- paste0(str_sub(ORF.i, start = 1, end = S.Start.ORF.i-1), paste(as.character(Variable.i), collapse = ""), str_sub(ORF.i, start = S.Start.ORF.i, end = nchar(ORF.i)))
        
        # search PTCs in frame
        
        StopCodons <- c("tag", "taa", "tga")
        
        # Because of Full Transcript Models only one option of frame
        
        bananasplit <- (strsplit(S.ORF.i, split = ""))[[1]]
        bananasplitxxi <- paste0(bananasplit[c(TRUE, FALSE, FALSE)], bananasplit[c(FALSE, TRUE, FALSE)], bananasplit[c(FALSE, FALSE, TRUE)])
        if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) > 1){
          Model.results[i,3] <- "PTC"
        }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 1){
          if(length(which(CodAn.isoform$V3 == "stop_codon")) == 1){
            Model.results[i,3] <- "CDS_change"  
          }else if(length(which(CodAn.isoform$V3 == "stop_codon")) == 0){
            Model.results[i,3] <- "PTC"
          }else {stop("\n Wait what \n")}
        }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 0){
          Model.results[i,3] <- "CDS_change"
        }
          
      }
      
    }else if(CodAn.isoform$V7[1] == "-"){ # if strand - -- 3'-CDS-5'
      
      Variable.i <- strsplit(tolower(as.character(reverseComplement(DNAString(paste(as.character(VariableSeq[[which(names(VariableSeq) == stri_sub(ID, 1, nchar(ID)-2))]]), collapse = ""))))), split = "")[[1]]
      
      # nested conditional just in case dont have UTRs
      # check changes from Transabyss transcript to transcript mod -- Transcriptmod == reverse complement (Transabyss transcrip) 
      
      if(length(which(names(fiveUTR.seq) == ID)) == 0 & length(which(names(threeUTR.seq) == ID)) == 0){
        Transcript.mod <- paste0(paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""))
      }else if(length(which(names(fiveUTR.seq) == ID)) == 0 & length(which(names(threeUTR.seq) == ID)) != 0){
        Transcript.mod <- paste0(paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(threeUTR.seq[[ which(names(threeUTR.seq) == ID)]]), collapse = ""))
      }else if(length(which(names(fiveUTR.seq) == ID)) != 0 & length(which(names(threeUTR.seq) == ID)) == 0){
        Transcript.mod <- paste0(paste(as.character(fiveUTR.seq[[ which(names(fiveUTR.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""))
      }else if(length(which(names(fiveUTR.seq) == ID)) != 0 & length(which(names(threeUTR.seq) == ID)) != 0){
        Transcript.mod <- paste0(paste(as.character(fiveUTR.seq[[ which(names(fiveUTR.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = ""), 
                                 paste(as.character(threeUTR.seq[[ which(names(threeUTR.seq) == ID)]]), collapse = "")
        )
      }                      
      
      if(!all(compare.DNA((strsplit(Transcript.mod, split = ""))[[1]], strsplit(tolower(as.character(reverseComplement(DNAString(paste(Transabyss[[which(names(Transabyss) == ID)]], collapse = ""))))), split = "")[[1]]))){
        # checking coherency of pasted transcripts parts
        stop(" \n Not coherency pasting all transcript parts \n")
      }
      
      if(stri_sub(ID, -2) == ".1"){ # 1) Search S part in upper (.1). If match search PTCs in frame ORF. If not match -- 2) info in .2? pass to .2; if not info in .2 paste S seq to .1 and search PTCs in frame
        
        S.mod <- stri_sub(Transcript.mod, S.Start-1, S.Start-1+length(Variable.i)-1)
        S.mod <- (strsplit(x = S.mod, split = ""))[[1]]
        
        if(length(S.mod) != length(Variable.i)){ stop("\n Different events variable seqs length \n")}
        
        Condition <- compare.DNA(S.mod, as.character(Variable.i))
        
        if(length(which(Condition == TRUE)) >= 0.95*length(Condition)){ # mapped transcript have this splicing event variation into account
          
          StopCodons <- c("tag", "taa", "tga")
          
          # search PTCs in frame
          
          ORF.i <- paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = "")
          
          # Because of Full Transcript Models only one option of frame
          
          bananasplit <- (strsplit(ORF.i, split = ""))[[1]]
          bananasplitxxi <- paste0(bananasplit[c(TRUE, FALSE, FALSE)], bananasplit[c(FALSE, TRUE, FALSE)], bananasplit[c(FALSE, FALSE, TRUE)])
          if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) > 1){
            Model.results[i,3] <- "PTC"
          }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 1){
            if(length(which(CodAn.isoform$V3 == "stop_codon")) == 1){
              Model.results[i,3] <- "CDS_change"  
            }else if(length(which(CodAn.isoform$V3 == "stop_codon")) == 0){
              Model.results[i,3] <- "PTC"
            }else {stop("\n Wait what \n")}
          }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 0){
            Model.results[i,3] <- "CDS_change"
          }
          
        }else if(length(which(Condition == TRUE)) < 0.95*length(Condition)){ # mapped transcript does not have this splicing events variation into account
            
            # paste S seq to ORF seq
            
            S.Start.ORF.i <- S.Start - min(CodAn.isoform$V4) # test with others
            ORF.i <- paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = "")
            S.ORF.i <- paste0(str_sub(ORF.i, start = 1, end = S.Start.ORF.i-1), paste(as.character(Variable.i), collapse = ""), str_sub(ORF.i, start = S.Start.ORF.i, end = nchar(ORF.i)))
            
            # search PTCs in frame
            
            StopCodons <- c("tag", "taa", "tga")
            
            # Because of Full Transcript Models only one option of frame
            
            bananasplit <- (strsplit(S.ORF.i, split = ""))[[1]]
            bananasplitxxi <- paste0(bananasplit[c(TRUE, FALSE, FALSE)], bananasplit[c(FALSE, TRUE, FALSE)], bananasplit[c(FALSE, FALSE, TRUE)])
            if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) > 1){
              Model.results[i,3] <- "PTC"
            }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 1){
              if(length(which(CodAn.isoform$V3 == "stop_codon")) == 1){
                Model.results[i,3] <- "CDS_change"  
              }else if(length(which(CodAn.isoform$V3 == "stop_codon")) == 0){
                Model.results[i,3] <- "PTC"
              }else {stop("\n Wait what \n")}
            }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 0){
              Model.results[i,3] <- "CDS_change"
            } # PTC conditional
          } # lower isoform contain also CDS info or not
          
      }else if(stri_sub(ID, -2) == ".2"){ 
        
        # paste S seq to ORF seq
        
        S.Start.ORF.i <- S.Start - min(CodAn.isoform$V4)
        ORF.i <- paste(as.character(ORFS.seq[[ which(names(ORFS.seq) == ID)]]), collapse = "")
        S.ORF.i <- paste0(str_sub(ORF.i, start = 1, end = S.Start.ORF.i-1), paste(as.character(Variable.i), collapse = ""), str_sub(ORF.i, start = S.Start.ORF.i, end = nchar(ORF.i)))
        
        # search PTCs in frame
        
        StopCodons <- c("tag", "taa", "tga")
        
        # Because of Full Transcript Models only one option of frame
        
        bananasplit <- (strsplit(S.ORF.i, split = ""))[[1]]
        bananasplitxxi <- paste0(bananasplit[c(TRUE, FALSE, FALSE)], bananasplit[c(FALSE, TRUE, FALSE)], bananasplit[c(FALSE, FALSE, TRUE)])
        if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) > 1){
          Model.results[i,3] <- "PTC"
        }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 1){
          if(length(which(CodAn.isoform$V3 == "stop_codon")) == 1){
            Model.results[i,3] <- "CDS_change"  
          }else if(length(which(CodAn.isoform$V3 == "stop_codon")) == 0){
            Model.results[i,3] <- "PTC"
          }else {stop("\n Wait what \n")}
          
        }else if(length(which(bananasplitxxi == StopCodons[1])) + length(which(bananasplitxxi == StopCodons[2])) + length(which(bananasplitxxi == StopCodons[3])) == 0){
          Model.results[i,3] <- "CDS_change"
        } # PTC conditional
        
      } # isoform conditional
      
      
    } # strand conditional

  } # CDS-PTC search 
} # for

## Event Coherency compara : fill Noprediction and NoMapping with the other isoform from same splicing event

x <- Model.results
  x[,1] <- as.character(x[,1])
  Identifications <- substr(x[,1], 1, nchar(x[,1])-2)
  x[,2] <- gsub("NoPrediction", NA, x[,2], fixed = T)
  x[,2] <- gsub("NoMapping", NA, x[,2], fixed = T)
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

NOCDSrows <- which(x$Model == "5UTR" | x$Model == "3UTR")
x[NOCDSrows[which(x[NOCDSrows,3] == "PTC" | x[NOCDSrows,3] == "CDS_change")],3] <- NA
write.table(x, file = "GlobalDesSplicingEvents.txt", row.names = F, col.names = T, sep = "\t")

# Global: Barplot, all, dv, de and dd (MODEL/w and wo NoRes ones and PTCs predict/CDS)

x <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/GlobalDesSplicingEvents.txt",  sep = "\t")
colnames(x) <- c("Isoforms", "Model", "PTC")
x <- x[-1,]
x$Model <- as.character(x$Model)
x$PTC <- as.character(x$PTC)
x[is.na(x$Model),2] <- "NoRes"
x[is.na(x$PTC),3] <- "NoRes"

# WITH/WO noRES and All, DV, DE, DD

All_CDS <- length(which(x$Model == "CDS" | x$Model == "StopCodon_CDS"))/length(x$Isoforms) * 100
All_3UTR <- length(which(x$Model == "3UTR"))/length(x$Isoforms) * 100
All_5UTR <- length(which(x$Model == "5UTR"))/length(x$Isoforms) * 100
All_NoRes <- length(which(x$Model == "NoRes"))/length(x$Isoforms) * 100

x.DVonly <- x[which(x$Isoforms %in% Isoforms_vector$KissDE_diff_only),]

DV_CDS <- length(which(x.DVonly$Model == "CDS" | x.DVonly$Model == "StopCodon_CDS"))/length(x.DVonly$Isoforms) * 100
DV_3UTR <- length(which(x.DVonly$Model == "3UTR"))/length(x.DVonly$Isoforms) * 100
DV_5UTR <- length(which(x.DVonly$Model == "5UTR"))/length(x.DVonly$Isoforms) * 100
DV_NoRes <- length(which(x.DVonly$Model == "NoRes"))/length(x.DVonly$Isoforms) * 100

x.DEonly <- x[which(x$Isoforms %in% Isoforms_vector$DESeq2_diff_only),]

DE_CDS <- length(which(x.DEonly$Model == "CDS" | x.DEonly$Model == "StopCodon_CDS"))/length(x.DEonly$Isoforms) * 100
DE_3UTR <- length(which(x.DEonly$Model == "3UTR"))/length(x.DEonly$Isoforms) * 100
DE_5UTR <- length(which(x.DEonly$Model == "5UTR"))/length(x.DEonly$Isoforms) * 100
DE_NoRes <- length(which(x.DEonly$Model == "NoRes"))/length(x.DEonly$Isoforms) * 100

x.DDonly <- x[which(x$Isoforms %in% Isoforms_vector$Intersection_doublediff),]

DD_CDS <- length(which(x.DDonly$Model == "CDS" | x.DDonly$Model == "StopCodon_CDS"))/length(x.DDonly$Isoforms) * 100
DD_3UTR <- length(which(x.DDonly$Model == "3UTR"))/length(x.DDonly$Isoforms) * 100
DD_5UTR <- length(which(x.DDonly$Model == "5UTR"))/length(x.DDonly$Isoforms) * 100
DD_NoRes <- length(which(x.DDonly$Model == "NoRes"))/length(x.DDonly$Isoforms) * 100

All.Mod <- data.frame(level = c(rep("All",4), rep("DS", 4), rep("DE", 4), rep("DD", 4)), variation = rep(c("CDS", "3UTR", "5UTR", "NoRes"),4), 
                      value = c(All_CDS, All_3UTR, All_5UTR, All_NoRes,
                                DV_CDS, DV_3UTR, DV_5UTR, DV_NoRes,
                                DE_CDS, DE_3UTR, DE_5UTR, DE_NoRes,
                                DD_CDS, DD_3UTR, DD_5UTR, DD_NoRes))

All.wnoRes<- ggplot(All.Mod, aes(fill=variation, y=value, x=level, color=variation, width=0.75), size=0.01) + 
  geom_bar(stat="identity", color="black", size=0.1)+theme_classic()+labs(x=NULL, y="Ratio")+
  scale_fill_npg()+scale_colour_manual(values=rep("#000000",11))+
  theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
        legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
        axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
        axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
        legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))

# wo NoRes

x.wiNoRes <- x
x.woNoRes <- x[which(x$Model != "NoRes"),]

All_CDS <- length(which(x.woNoRes$Model == "CDS" | x.woNoRes$Model == "StopCodon_CDS"))/length(x.woNoRes$Isoforms) * 100
All_3UTR <- length(which(x.woNoRes$Model == "3UTR"))/length(x.woNoRes$Isoforms) * 100
All_5UTR <- length(which(x.woNoRes$Model == "5UTR"))/length(x.woNoRes$Isoforms) * 100

x.woNoRes.DVonly <- x.woNoRes[which(x.woNoRes$Isoforms %in% Isoforms_vector$KissDE_diff_only),]

DV_CDS <- length(which(x.woNoRes.DVonly$Model == "CDS" | x.woNoRes.DVonly$Model == "StopCodon_CDS"))/length(x.woNoRes.DVonly$Isoforms) * 100
DV_3UTR <- length(which(x.woNoRes.DVonly$Model == "3UTR"))/length(x.woNoRes.DVonly$Isoforms) * 100
DV_5UTR <- length(which(x.woNoRes.DVonly$Model == "5UTR"))/length(x.woNoRes.DVonly$Isoforms) * 100

x.woNoRes.DEonly <- x.woNoRes[which(x.woNoRes$Isoforms %in% Isoforms_vector$DESeq2_diff_only),]

DE_CDS <- length(which(x.woNoRes.DEonly$Model == "CDS" | x.woNoRes.DEonly$Model == "StopCodon_CDS"))/length(x.woNoRes.DEonly$Isoforms) * 100
DE_3UTR <- length(which(x.woNoRes.DEonly$Model == "3UTR"))/length(x.woNoRes.DEonly$Isoforms) * 100
DE_5UTR <- length(which(x.woNoRes.DEonly$Model == "5UTR"))/length(x.woNoRes.DEonly$Isoforms) * 100

x.woNoRes.DDonly <- x.woNoRes[which(x.woNoRes$Isoforms %in% Isoforms_vector$Intersection_doublediff),]

DD_CDS <- length(which(x.woNoRes.DDonly$Model == "CDS" | x.woNoRes.DDonly$Model == "StopCodon_CDS"))/length(x.woNoRes.DDonly$Isoforms) * 100
DD_3UTR <- length(which(x.woNoRes.DDonly$Model == "3UTR"))/length(x.woNoRes.DDonly$Isoforms) * 100
DD_5UTR <- length(which(x.woNoRes.DDonly$Model == "5UTR"))/length(x.woNoRes.DDonly$Isoforms) * 100

All.Mod <- data.frame(level = c(rep("All",3), rep("DS", 3), rep("DE", 3), rep("DD", 3)), variation = rep(c("CDS", "3UTR", "5UTR"),4), 
                      value = c(All_CDS, All_3UTR, All_5UTR, 
                                DV_CDS, DV_3UTR, DV_5UTR, 
                                DE_CDS, DE_3UTR, DE_5UTR, 
                                DD_CDS, DD_3UTR, DD_5UTR))

All.wonoRes<- ggplot(All.Mod, aes(fill=variation, y=value, x=level, color=variation, width=0.75), size=0.01) + 
  geom_bar(stat="identity", color="black", size=0.1)+theme_classic()+labs(x=NULL, y="Ratio")+
  scale_fill_npg()+scale_colour_manual(values=rep("#000000",11))+
  theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
        legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
        axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
        axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
        legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))

# PTCs CDS Change

pal_bar <- brewer.pal(3, "Dark2")

x.PTC <- x[which(x$Model == "CDS" | x$Model == "StopCodon_CDS"),]

All_PTC <- length(which(x.PTC$PTC == "PTC"))/length(x.PTC$Isoforms) * 100
All_CDS_change <-length(which(x.PTC$PTC == "CDS_change"))/length(x.PTC$Isoforms) * 100

x.PTC.DVonly <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$KissDE_diff_only),] 

DV_PTC <- length(which(x.PTC.DVonly$PTC == "PTC"))/length(x.PTC.DVonly$Isoforms) * 100
DV_CDS_chage <- length(which(x.PTC.DVonly$PTC == "CDS_change"))/length(x.PTC.DVonly$Isoforms) * 100

x.PTC.DEonly <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$DESeq2_diff_only),] 

DE_PTC <- length(which(x.PTC.DEonly$PTC == "PTC"))/length(x.PTC.DEonly$Isoforms) * 100
DE_CDS_chage <- length(which(x.PTC.DEonly$PTC == "CDS_change"))/length(x.PTC.DEonly$Isoforms) * 100

x.PTC.DDonly <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$Intersection_doublediff),] 

DD_PTC <- length(which(x.PTC.DDonly$PTC == "PTC"))/length(x.PTC.DDonly$Isoforms) * 100
DD_CDS_chage <- length(which(x.PTC.DDonly$PTC == "CDS_change"))/length(x.PTC.DDonly$Isoforms) * 100

All.PTC <- data.frame(level = c(rep("All",2), rep("DS", 2), rep("DE", 2), rep("DD", 2)), variation = rep(c("PTC", "Change"),4), 
                      value = c(All_PTC, All_CDS_change,
                                DV_PTC, DV_CDS_chage,
                                DE_PTC, DE_CDS_chage,
                                DD_PTC, DD_CDS_chage))

All.PTC.plot <- ggplot(All.PTC, aes(fill=variation, y=value, x=level, color=variation, width=0.75), size=0.01) + 
  geom_bar(stat="identity", color="black", size=0.1)+theme_classic()+labs(x=NULL, y="Ratio")+
  scale_fill_manual(values=rev(pal_bar))+scale_colour_manual(values=rep("#000000",11))+
  theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
        legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
        axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
        axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
        legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))

pdf(file = "GlobalDescriptive_AllDVDEDD.pdf", paper = "special", height = 1.1, width = 2.2, onefile = T)

All.wnoRes
All.wonoRes
All.PTC.plot

dev.off()


# WITH/WO noRES ad All, DV (CT1-CT3), DE (CT1-CT3), DD (CT1-CT3) and do with T1T3 for supps

All_CDS <- length(which(x$Model == "CDS" | x$Model == "StopCodon_CDS"))/length(x$Isoforms) * 100
All_3UTR <- length(which(x$Model == "3UTR"))/length(x$Isoforms) * 100
All_5UTR <- length(which(x$Model == "5UTR"))/length(x$Isoforms) * 100
All_NoRes <- length(which(x$Model == "NoRes"))/length(x$Isoforms) * 100

x.DVonly.CT1 <- x[which(x$Isoforms %in% Isoforms_vector$KissDE_diff_only[Isoforms_vector$KissDE_diff_only %in% c(paste0(KissDE_diff$CT1$`1.ID`,".1"), paste0(KissDE_diff$CT1$`1.ID`,".2"))]),]

DV_CDS_CT1 <- length(which(x.DVonly.CT1$Model == "CDS" | x.DVonly.CT1$Model == "StopCodon_CDS"))/length(x.DVonly.CT1$Isoforms) * 100
DV_3UTR_CT1 <- length(which(x.DVonly.CT1$Model == "3UTR"))/length(x.DVonly.CT1$Isoforms) * 100
DV_5UTR_CT1 <- length(which(x.DVonly.CT1$Model == "5UTR"))/length(x.DVonly.CT1$Isoforms) * 100
DV_NoRes_CT1 <- length(which(x.DVonly.CT1$Model == "NoRes"))/length(x.DVonly.CT1$Isoforms) * 100

x.DVonly.CT3 <- x[which(x$Isoforms %in% Isoforms_vector$KissDE_diff_only[Isoforms_vector$KissDE_diff_only %in% c(paste0(KissDE_diff$CT3$`1.ID`,".1"), paste0(KissDE_diff$CT3$`1.ID`,".2"))]),]

DV_CDS_CT3 <- length(which(x.DVonly.CT3$Model == "CDS" | x.DVonly.CT3$Model == "StopCodon_CDS"))/length(x.DVonly.CT3$Isoforms) * 100
DV_3UTR_CT3 <- length(which(x.DVonly.CT3$Model == "3UTR"))/length(x.DVonly.CT3$Isoforms) * 100
DV_5UTR_CT3 <- length(which(x.DVonly.CT3$Model == "5UTR"))/length(x.DVonly.CT3$Isoforms) * 100
DV_NoRes_CT3 <- length(which(x.DVonly.CT3$Model == "NoRes"))/length(x.DVonly.CT3$Isoforms) * 100

x.DVonly.T1T3 <- x[which(x$Isoforms %in% Isoforms_vector$KissDE_diff_only[Isoforms_vector$KissDE_diff_only %in% c(paste0(KissDE_diff$T1T3$`1.ID`,".1"), paste0(KissDE_diff$T1T3$`1.ID`,".2"))]),]

DV_CDS_T1T3 <- length(which(x.DVonly.T1T3$Model == "CDS" | x.DVonly.T1T3$Model == "StopCodon_CDS"))/length(x.DVonly.T1T3$Isoforms) * 100
DV_3UTR_T1T3 <- length(which(x.DVonly.T1T3$Model == "3UTR"))/length(x.DVonly.T1T3$Isoforms) * 100
DV_5UTR_T1T3 <- length(which(x.DVonly.T1T3$Model == "5UTR"))/length(x.DVonly.T1T3$Isoforms) * 100
DV_NoRes_T1T3 <- length(which(x.DVonly.T1T3$Model == "NoRes"))/length(x.DVonly.T1T3$Isoforms) * 100

##

x.DEonly.CT1 <- x[which(x$Isoforms %in% Isoforms_vector$DESeq2_diff_only[Isoforms_vector$DESeq2_diff_only %in% DESeq2_diff_only$CT1$KissID]),]

DE_CDS_CT1 <- length(which(x.DEonly.CT1$Model == "CDS" | x.DEonly.CT1$Model == "StopCodon_CDS"))/length(x.DEonly.CT1$Isoforms) * 100
DE_3UTR_CT1 <- length(which(x.DEonly.CT1$Model == "3UTR"))/length(x.DEonly.CT1$Isoforms) * 100
DE_5UTR_CT1 <- length(which(x.DEonly.CT1$Model == "5UTR"))/length(x.DEonly.CT1$Isoforms) * 100
DE_NoRes_CT1 <- length(which(x.DEonly.CT1$Model == "NoRes"))/length(x.DEonly.CT1$Isoforms) * 100

x.DEonly.CT3 <- x[which(x$Isoforms %in% Isoforms_vector$DESeq2_diff_only[Isoforms_vector$DESeq2_diff_only %in% DESeq2_diff_only$CT3$KissID]),]

DE_CDS_CT3 <- length(which(x.DEonly.CT3$Model == "CDS" | x.DEonly.CT3$Model == "StopCodon_CDS"))/length(x.DEonly.CT3$Isoforms) * 100
DE_3UTR_CT3 <- length(which(x.DEonly.CT3$Model == "3UTR"))/length(x.DEonly.CT3$Isoforms) * 100
DE_5UTR_CT3 <- length(which(x.DEonly.CT3$Model == "5UTR"))/length(x.DEonly.CT3$Isoforms) * 100
DE_NoRes_CT3 <- length(which(x.DEonly.CT3$Model == "NoRes"))/length(x.DEonly.CT3$Isoforms) * 100

x.DEonly.T1T3 <- x[which(x$Isoforms %in% Isoforms_vector$DESeq2_diff_only[Isoforms_vector$DESeq2_diff_only %in% DESeq2_diff_only$T1T3$KissID]),]

DE_CDS.T1T3 <- length(which(x.DEonly.T1T3$Model == "CDS" | x.DEonly.T1T3$Model == "StopCodon_CDS"))/length(x.DEonly.T1T3$Isoforms) * 100
DE_3UTR.T1T3 <- length(which(x.DEonly.T1T3$Model == "3UTR"))/length(x.DEonly.T1T3$Isoforms) * 100
DE_5UTR.T1T3 <- length(which(x.DEonly.T1T3$Model == "5UTR"))/length(x.DEonly.T1T3$Isoforms) * 100
DE_NoRes.T1T3 <- length(which(x.DEonly.T1T3$Model == "NoRes"))/length(x.DEonly.T1T3$Isoforms) * 100

##

KissDE.CT1.isoforms <- c(paste0(KissDE_diff$CT1$`1.ID`, ".1"),paste0(KissDE_diff$CT1$`1.ID`, ".2"))
DD.CT1 <- KissDE.CT1.isoforms[KissDE.CT1.isoforms %in% DESeq2_diff$CT1$KissID]
x.DDonly.CT1 <- x[which(x$Isoforms %in% Isoforms_vector$Intersection_doublediff[Isoforms_vector$Intersection_doublediff %in% DD.CT1]),]

DD_CDS.CT1 <- length(which(x.DDonly.CT1$Model == "CDS" | x.DDonly.CT1$Model == "StopCodon_CDS"))/length(x.DDonly.CT1$Isoforms) * 100
DD_3UTR.CT1 <- length(which(x.DDonly.CT1$Model == "3UTR"))/length(x.DDonly.CT1$Isoforms) * 100
DD_5UTR.CT1 <- length(which(x.DDonly.CT1$Model == "5UTR"))/length(x.DDonly.CT1$Isoforms) * 100
DD_NoRes.CT1 <- length(which(x.DDonly.CT1$Model == "NoRes"))/length(x.DDonly.CT1$Isoforms) * 100

KissDE.CT3.isoforms <- c(paste0(KissDE_diff$CT3$`1.ID`, ".1"),paste0(KissDE_diff$CT3$`1.ID`, ".2"))
DD.CT3 <- KissDE.CT3.isoforms[KissDE.CT3.isoforms %in% DESeq2_diff$CT3$KissID]
x.DDonly.CT3 <- x[which(x$Isoforms %in% Isoforms_vector$Intersection_doublediff[Isoforms_vector$Intersection_doublediff %in% DD.CT3]),]

DD_CDS.CT3 <- length(which(x.DDonly.CT3$Model == "CDS" | x.DDonly.CT3$Model == "StopCodon_CDS"))/length(x.DDonly.CT3$Isoforms) * 100
DD_3UTR.CT3 <- length(which(x.DDonly.CT3$Model == "3UTR"))/length(x.DDonly.CT3$Isoforms) * 100
DD_5UTR.CT3 <- length(which(x.DDonly.CT3$Model == "5UTR"))/length(x.DDonly.CT3$Isoforms) * 100
DD_NoRes.CT3 <- length(which(x.DDonly.CT3$Model == "NoRes"))/length(x.DDonly.CT3$Isoforms) * 100

KissDE.T1T3.isoforms <- c(paste0(KissDE_diff$T1T3$`1.ID`, ".1"),paste0(KissDE_diff$T1T3$`1.ID`, ".2"))
DD.T1T3 <- KissDE.T1T3.isoforms[KissDE.T1T3.isoforms %in% DESeq2_diff$T1T3$KissID]
x.DDonly.T1T3 <- x[which(x$Isoforms %in% Isoforms_vector$Intersection_doublediff[Isoforms_vector$Intersection_doublediff %in% DD.T1T3]),]

DD_CDS.T1T3 <- length(which(x.DDonly.T1T3$Model == "CDS" | x.DDonly.T1T3$Model == "StopCodon_CDS"))/length(x.DDonly.T1T3$Isoforms) * 100
DD_3UTR.T1T3 <- length(which(x.DDonly.T1T3$Model == "3UTR"))/length(x.DDonly.T1T3$Isoforms) * 100
DD_5UTR.T1T3 <- length(which(x.DDonly.T1T3$Model == "5UTR"))/length(x.DDonly.T1T3$Isoforms) * 100
DD_NoRes.T1T3 <- length(which(x.DDonly.T1T3$Model == "NoRes"))/length(x.DDonly.T1T3$Isoforms) * 100

##

All.Mod.desglose <- data.frame(level = c(rep("All",4), rep("DS_CT1", 4), rep("DS_CT3", 4), rep("DS_T1T3", 4),rep("DE_CT1", 4), rep("DE_CT3", 4), rep("DE_T1T3", 4),rep("DD_CT1", 4), rep("DD_CT3", 4), rep("DD_T1T3", 4)), variation = rep(c("CDS", "3UTR", "5UTR", "NoRes"),10), 
                      value = c(All_CDS, All_3UTR, All_5UTR, All_NoRes,
                                DV_CDS_CT1, DV_3UTR_CT1, DV_5UTR_CT1, DV_NoRes_CT1,
                                DV_CDS_CT3, DV_3UTR_CT3, DV_5UTR_CT3, DV_NoRes_CT3,
                                DV_CDS_T1T3, DV_3UTR_T1T3, DV_5UTR_T1T3, DV_NoRes_T1T3,
                                DE_CDS_CT1, DE_3UTR_CT1, DE_5UTR_CT1, DE_NoRes_CT1,
                                DE_CDS_CT3, DE_3UTR_CT3, DE_5UTR_CT3, DE_NoRes_CT3,
                                DE_CDS.T1T3, DE_3UTR.T1T3, DE_5UTR.T1T3, DE_NoRes.T1T3,
                                DD_CDS.CT1, DD_3UTR.CT1, DD_5UTR.CT1, DD_NoRes.CT1,
                                DD_CDS.CT3, DD_3UTR.CT3, DD_5UTR.CT3, DD_NoRes.CT3,
                                DD_CDS.T1T3, DD_3UTR.T1T3, DD_5UTR.T1T3, DD_NoRes.T1T3
                                ))

All.desglose<- ggplot(All.Mod.desglose, aes(fill=variation, y=value, x=level, color=variation, width=0.75), size=0.01) + 
  geom_bar(stat="identity", color="black", size=0.1)+theme_classic()+labs(x=NULL, y="Ratio")+
  scale_fill_npg()+scale_colour_manual(values=rep("#000000",11))+
  theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
        legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
        axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
        axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
        legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))



## PTC desglose


x.PTC <- x[which(x$Model == "CDS" | x$Model == "StopCodon_CDS"),]

All_PTC <- length(which(x.PTC$PTC == "PTC"))/length(x.PTC$Isoforms) * 100
All_CDS_change <-length(which(x.PTC$PTC == "CDS_change"))/length(x.PTC$Isoforms) * 100

x.PTC.DVonly.CT1 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$KissDE_diff_only[Isoforms_vector$KissDE_diff_only %in% c(paste0(KissDE_diff$CT1$`1.ID`,".1"), paste0(KissDE_diff$CT1$`1.ID`,".2"))]),]

DV_PTC_CT1.PTC <- length(which(x.PTC.DVonly.CT1$PTC == "PTC"))/length(x.PTC.DVonly.CT1$Isoforms) * 100
DV_CDSchange_CT1.PTC <- length(which(x.PTC.DVonly.CT1$PTC == "CDS_change"))/length(x.PTC.DVonly.CT1$Isoforms) * 100

x.PTC.DVonly.CT3 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$KissDE_diff_only[Isoforms_vector$KissDE_diff_only %in% c(paste0(KissDE_diff$CT3$`1.ID`,".1"), paste0(KissDE_diff$CT3$`1.ID`,".2"))]),]

DV_PTC_CT3.PTC <- length(which(x.PTC.DVonly.CT3$PTC == "PTC"))/length(x.PTC.DVonly.CT3$Isoforms) * 100
DV_CDSchange_CT3.PTC <- length(which(x.PTC.DVonly.CT3$PTC == "CDS_change"))/length(x.PTC.DVonly.CT3$Isoforms) * 100

x.PTC.DVonly.T1T3 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$KissDE_diff_only[Isoforms_vector$KissDE_diff_only %in% c(paste0(KissDE_diff$T1T3$`1.ID`,".1"), paste0(KissDE_diff$T1T3$`1.ID`,".2"))]),]

DV_PTC_T1T3.PTC <- length(which(x.PTC.DVonly.T1T3$PTC == "PTC"))/length(x.PTC.DVonly.T1T3$Isoforms) * 100
DV_CDSchange_T1T3.PTC <- length(which(x.PTC.DVonly.T1T3$PTC == "CDS_change"))/length(x.PTC.DVonly.T1T3$Isoforms) * 100

##

x.PTC.DEonly.CT1 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$DESeq2_diff_only[Isoforms_vector$DESeq2_diff_only %in% DESeq2_diff$CT1$KissID]),] 

DE_PTC_CT1.PTC <- length(which(x.PTC.DEonly.CT1$PTC == "PTC"))/length(x.PTC.DEonly.CT1$Isoforms) * 100
DE_CDSchange_CT1.PTC <- length(which(x.PTC.DEonly.CT1$PTC == "CDS_change"))/length(x.PTC.DEonly.CT1$Isoforms) * 100

x.PTC.DEonly.CT3 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$DESeq2_diff_only[Isoforms_vector$DESeq2_diff_only %in% DESeq2_diff$CT3$KissID]),] 

DE_PTC_CT3.PTC <- length(which(x.PTC.DEonly.CT3$PTC == "PTC"))/length(x.PTC.DEonly.CT3$Isoforms) * 100
DE_CDSchange_CT3.PTC <- length(which(x.PTC.DEonly.CT3$PTC == "CDS_change"))/length(x.PTC.DEonly.CT3$Isoforms) * 100

x.PTC.DEonly.T1T3 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$DESeq2_diff_only[Isoforms_vector$DESeq2_diff_only %in% DESeq2_diff$T1T3$KissID]),] 

DE_PTC_T1T3.PTC <- length(which(x.PTC.DEonly.T1T3$PTC == "PTC"))/length(x.PTC.DEonly.T1T3$Isoforms) * 100
DE_CDSchange_T1T3.PTC <- length(which(x.PTC.DEonly.T1T3$PTC == "CDS_change"))/length(x.PTC.DEonly.T1T3$Isoforms) * 100

##

KissDE.CT1.isoforms <- c(paste0(KissDE_diff$CT1$`1.ID`, ".1"),paste0(KissDE_diff$CT1$`1.ID`, ".2"))
DD.CT1 <- KissDE.CT1.isoforms[KissDE.CT1.isoforms %in% DESeq2_diff$CT1$KissID]
x.PTC.DDonly.CT1 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$Intersection_doublediff[Isoforms_vector$Intersection_doublediff %in% DD.CT1]),] 

DD_PTC.CT1 <- length(which(x.PTC.DDonly.CT1$PTC == "PTC"))/length(x.PTC.DDonly.CT1$Isoforms) * 100
DD_CDS_change.CT1 <- length(which(x.PTC.DDonly.CT1$PTC == "CDS_change"))/length(x.PTC.DDonly.CT1$Isoforms) * 100

KissDE.CT3.isoforms <- c(paste0(KissDE_diff$CT3$`1.ID`, ".1"),paste0(KissDE_diff$CT3$`1.ID`, ".2"))
DD.CT3 <- KissDE.CT3.isoforms[KissDE.CT3.isoforms %in% DESeq2_diff$CT3$KissID]
x.PTC.DDonly.CT3 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$Intersection_doublediff[Isoforms_vector$Intersection_doublediff %in% DD.CT3]),] 

DD_PTC.CT3 <- length(which(x.PTC.DDonly.CT3$PTC == "PTC"))/length(x.PTC.DDonly.CT3$Isoforms) * 100
DD_CDS_change.CT3 <- length(which(x.PTC.DDonly.CT3$PTC == "CDS_change"))/length(x.PTC.DDonly.CT3$Isoforms) * 100

KissDE.T1T3.isoforms <- c(paste0(KissDE_diff$T1T3$`1.ID`, ".1"),paste0(KissDE_diff$T1T3$`1.ID`, ".2"))
DD.T1T3 <- KissDE.T1T3.isoforms[KissDE.T1T3.isoforms %in% DESeq2_diff$T1T3$KissID]
x.PTC.DDonly.T1T3 <- x.PTC[which(x.PTC$Isoforms %in% Isoforms_vector$Intersection_doublediff[Isoforms_vector$Intersection_doublediff %in% DD.T1T3]),] 

DD_PTC.T1T3 <- length(which(x.PTC.DDonly.T1T3$PTC == "PTC"))/length(x.PTC.DDonly.T1T3$Isoforms) * 100
DD_CDS_change.T1T3 <- length(which(x.PTC.DDonly.T1T3$PTC == "CDS_change"))/length(x.PTC.DDonly.T1T3$Isoforms) * 100

All.Mod.desglose.PTC <- data.frame(level = c(rep("All",2), rep("DS_CT1", 2), rep("DS_CT3", 2), rep("DS_T1T3", 2),rep("DE_CT1", 2), rep("DE_CT3", 2), rep("DE_T1T3", 2),rep("DD_CT1", 2), rep("DD_CT3", 2), rep("DD_T1T3", 2)), variation = rep(c("Change", "PTC"),10), 
                               value = c(All_CDS_change, All_PTC,
                                         DV_CDSchange_CT1.PTC, DV_PTC_CT1.PTC,
                                         DV_CDSchange_CT3.PTC, DV_PTC_CT3.PTC,
                                         DV_CDSchange_T1T3.PTC, DV_PTC_T1T3.PTC,
                                         DE_CDSchange_CT1.PTC, DE_PTC_CT1.PTC,
                                         DE_CDSchange_CT3.PTC, DE_PTC_CT3.PTC,
                                         DE_CDSchange_T1T3.PTC, DE_PTC_T1T3.PTC,
                                         DD_CDS_change.CT1, DD_PTC.CT1,
                                         DD_CDS_change.CT3, DD_PTC.CT3,
                                         DD_CDS_change.T1T3, DD_PTC.T1T3
                               ))

All.desglose.PTC<- ggplot(All.Mod.desglose.PTC, aes(fill=variation, y=value, x=level, color=variation, width=0.75), size=0.01) + 
  geom_bar(stat="identity", color="black", size=0.1)+theme_classic()+labs(x=NULL, y="Ratio")+
  scale_fill_jama()+scale_colour_manual(values=rep("#000000",11))+
  theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
        legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
        axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
        axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
        legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))



pdf(file = "GlobalDescriptive_AllDVDEDD_desglose.pdf", paper = "special", height = 1.1, width = 2.2, onefile = T)

All.desglose
All.desglose.PTC

dev.off()

# Kmeans Global: Piechart for each Kmeans cluster (MODEL/w and wo NoRes ones and PTCs predict/CDS)

# Import def Kmeans -- norm counts 16 clusters

k16.kmeas <- read.delim("clipboard", header = T) # from 2.kmeans_all

donuts <- list()

for(i in levels(factor(k16.kmeas$groups))){
  cluster <- k16.kmeas[which(k16.kmeas$groups == i),]
  cluster.models <- x[x$Isoforms %in% cluster$ID,]
  cluster.ptc <- x.PTC[x.PTC$Isoforms %in% cluster$ID,]
  
  CDS <- length(which(cluster.models$Model == "CDS" | cluster.models$Model == "StopCodon_CDS"))/length(cluster.models$Isoforms) * 100
  threeUTR <- length(which(cluster.models$Model == "3UTR"))/length(cluster.models$Isoforms) * 100
  fiveUTR <- length(which(cluster.models$Model == "5UTR"))/length(cluster.models$Isoforms) * 100
  NoRes <- length(which(cluster.models$Model == "NoRes"))/length(cluster.models$Isoforms) * 100
  
  PTC <- length(which(cluster.ptc$PTC == "PTC"))/length(cluster.ptc$Isoforms) * 100
  Change <- length(which(cluster.ptc$PTC == "CDS_change"))/length(cluster.ptc$Isoforms) * 100
  
  mood <- data.frame(category = c("CDS", "3UTR", "5UTR", "NoRes"),
                     count = c(CDS, threeUTR, fiveUTR, NoRes))
  mood$ymax <- cumsum(mood$count)
  mood$ymin <- c(0, head(mood$ymax, n=-1))
  mood$labelPosition <- (mood$ymax + mood$ymin) / 2
  mood$label <- paste0(mood$category, "\n value: ", mood$count)
  mood$level <- rep(paste0("Cluster ", i),4)
  
  ptcs <- data.frame(category = c("Change", "PTC"),
                     count = c(Change, PTC))
  ptcs$ymax <- cumsum(ptcs$count)
  ptcs$ymin <- c(0, head(ptcs$ymax, n=-1))
  ptcs$labelPosition <- (ptcs$ymax + ptcs$ymin) / 2
  ptcs$label <- paste0(ptcs$category, "\n value: ", ptcs$count)
  ptcs$level <- rep(paste0("Cluster ", i),2)
  
  donut.mood <- ggplot(mood, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=3) +
    scale_fill_npg() +
    scale_color_npg() +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "none") + ggtitle(paste0("Cluster ",i))
  donuts[[i]][['Model']] <- donut.mood
  
  donut.ptcs <- ggplot(ptcs, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=3) +
    scale_fill_jama() +
    scale_color_jama() +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "none") + ggtitle(paste0("Cluster ",i))
  donuts[[i]][['PTCs']] <- donut.ptcs 
}


pdf("GlobalDescriptive_kmeans16.pdf")
grid.arrange(donuts$`1`$Model, donuts$`2`$Model, donuts$`3`$Model, donuts$`4`$Model,
             donuts$`5`$Model, donuts$`6`$Model, donuts$`7`$Model, donuts$`8`$Model,
             donuts$`9`$Model, donuts$`10`$Model, donuts$`11`$Model, donuts$`12`$Model,
             donuts$`13`$Model, donuts$`14`$Model, donuts$`15`$Model, donuts$`16`$Model, 
             ncol = 4, nrow = 4)

grid.arrange(donuts$`1`$PTCs, donuts$`2`$PTCs, donuts$`3`$PTCs, donuts$`4`$PTCs,
             donuts$`5`$PTCs, donuts$`6`$PTCs, donuts$`7`$PTCs, donuts$`8`$PTCs,
             donuts$`9`$PTCs, donuts$`10`$PTCs, donuts$`11`$PTCs, donuts$`12`$PTCs,
             donuts$`13`$PTCs, donuts$`14`$PTCs, donuts$`15`$PTCs, donuts$`16`$PTCs, 
             ncol = 4, nrow = 4)
dev.off()

## Transcripts models: Candidates (fourth) (after MM)

Candidates.KissIDs <- read.delim("clipboard", header = T)

Candidatesf <- c(paste0(Candidates.KissIDs$KissID, ".1"), paste0(Candidates.KissIDs$KissID, ".2"))[order(c(paste0(Candidates.KissIDs$KissID, ".1"), paste0(Candidates.KissIDs$KissID, ".2")), decreasing = T)]

CodAn.Candidates <- CodAn.gtf[which(CodAn.gtf$V1 %in% Candidatesf),]

# Conservative approach: take lower (.2) isoform insert variabe part ant grab coordinates

# CHMP1A

chmp1a.lower <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_9458|Cycle_1.1"),]
chmp1a.lower <- chmp1a.lower[,c(3,4,5)]
length(Transabyss[[which(names(Transabyss) == "bcc_9458|Cycle_1.1")]])
chmp1a.lower.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                              coordinates = c("1-404", "405-980", "980-1175"))
chmp1a.upper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                              coordinates = c("1-404", "405-980", "980-1423"))

par(mfrow=c(2,1))
genemodel.plot(chmp1a.upper.df, start=1, bpstop=1423, orientation="reverse", xaxis=T)
mutation.plot(1116, 1117, text="+152 nt", col="black", drop=-.15, haplotypes=c("red", "blue"))
genemodel.plot(chmp1a.lower.df, start=1, bpstop=1175, orientation="reverse", xaxis=F)

# HDT2

hdt2.lower <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_3998|Cycle_0.2"),]
hdt2.upper <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_3998|Cycle_0.1"),]

length(Transabyss[[which(names(Transabyss) == "bcc_3998|Cycle_0.2")]])
length(Transabyss[[which(names(Transabyss) == "bcc_3998|Cycle_0.1")]])

hdt2.lower.df <- data.frame(type = c("coding_region", "5' utr"),
                            coordinates = c("1-710", "711-830"))
hdt2.upper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                         coordinates = c("1-497", "498-1244", "1245-1787"))

hdt2.lowerupper.df <- data.frame(type = c("coding_region", "5' utr"),
                            coordinates = c("1-824", "825-944"))

length(VariableSeq[[which(names(VariableSeq) == "bcc_3998|Cycle_0")]])
bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == "bcc_3998|Cycle_0.2")]

par(mfrow=c(2,1))
genemodel.plot(hdt2.upper.df, start=1, bpstop=1787, orientation="reverse", xaxis=T)
genemodel.plot(hdt2.lower.df, start=1, bpstop=830, orientation="reverse", xaxis=F)

par(mfrow=c(2,1))
genemodel.plot(hdt2.lowerupper.df, start=1, bpstop=944, orientation="reverse", xaxis=T)
mutation.plot(407, 408, text="+114 nt", col="black", drop=-.15, haplotypes="red")
mutation.plot(466, 467, text="PTC", col="black", drop=-.15, haplotypes="blue")
genemodel.plot(hdt2.lower.df, start=1, bpstop=830, orientation="reverse", xaxis=F)

# MgProto

mgproto.lower <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_17045|Cycle_0.2"),]
mgproto.upper <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_17045|Cycle_0.1"),]

length(Transabyss[[which(names(Transabyss) == "bcc_17045|Cycle_0.2")]])
length(Transabyss[[which(names(Transabyss) == "bcc_17045|Cycle_0.1")]])

mgproto.lower.df <- data.frame(type = c("5' utr", "coding_region"),
                            coordinates = c("1-18", "19-678"))
mgproto.upper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                            coordinates = c("1-704", "705-1691", "1691-1816"))

length(VariableSeq[[which(names(VariableSeq) == "bcc_17045|Cycle_0")]])
bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == "bcc_17045|Cycle_0.2")]

mgproto.lowerupper.df <- data.frame(type = c("5' utr", "coding_region"),
                               coordinates = c("1-18", "19-786"))

x[x$Isoforms == "bcc_17045|Cycle_0.2",]

par(mfrow=c(2,1))
genemodel.plot(mgproto.upper.df, start=1, bpstop=1816, orientation="reverse", xaxis=T)
genemodel.plot(mgproto.lower.df, start=1, bpstop=678, orientation="reverse", xaxis=T)

par(mfrow=c(2,1))
genemodel.plot(mgproto.lowerupper.df, start=1, bpstop=786, orientation="reverse", xaxis=T)
mutation.plot(69, 70, text="+108 nt", col="black", drop=-.15, haplotypes="red")
mutation.plot(99, 100, text="PTC", col="black", drop=-.15, haplotypes="blue")
genemodel.plot(mgproto.lower.df, start=1, bpstop=678, orientation="reverse", xaxis=F)

# DUF4050: no lower prediction for this candidate 

duf.upper <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_29432|Cycle_2.1"),]

length(Transabyss[[which(names(Transabyss) == "bcc_29432|Cycle_2.1")]])

duf.upper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                               coordinates = c("1-1014", "1015-1500", "1501-1796"))

length(VariableSeq[[which(names(VariableSeq) == "bcc_29432|Cycle_2")]])
bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == "bcc_29432|Cycle_2.1")]


x[x$Isoforms == "bcc_29432|Cycle_2.1",]

par(mfrow=c(2,1))
genemodel.plot(duf.upper.df, start=1, bpstop=1796, orientation="reverse", xaxis=T)
mutation.plot(295+72, 295+72, text="+88 nt", col="black", drop=-.15, haplotypes="red")
mutation.plot(295+62+72, 295+62+72, text="PTC", col="black", drop=-.15, haplotypes="blue")

# RSZ22

rsz.lower <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_9268|Cycle_9.2"),]
rsz.upper <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_9268|Cycle_9.1"),]

length(Transabyss[[which(names(Transabyss) == "bcc_9268|Cycle_9.2")]])
length(Transabyss[[which(names(Transabyss) == "bcc_9268|Cycle_9.1")]])

rsz.lower.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                               coordinates = c("1-766", "767-1330", "1331-1627"))
rsz.upper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                               coordinates = c("1-929", "930-1598", "1599-1823"))

length(VariableSeq[[which(names(VariableSeq) == "bcc_9268|Cycle_9")]])
bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == "bcc_9268|Cycle_9.2")]

rsz.lowerupper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                                    coordinates = c("1-766", "767-1435", "1436-1732"))

x[x$Isoforms == "bcc_9268|Cycle_9.2",]

par(mfrow=c(2,1))
genemodel.plot(rsz.upper.df, start=1, bpstop=1823, orientation="reverse", xaxis=T)
genemodel.plot(rsz.lower.df, start=1, bpstop=1627, orientation="reverse", xaxis=T)

par(mfrow=c(2,1))
genemodel.plot(rsz.lowerupper.df, start=1, bpstop=1732, orientation="reverse", xaxis=T)
mutation.plot(808, 809, text="+105 nt", col="black", drop=-.15, haplotypes="red")
mutation.plot(808+50, 809+50, text="PTC", col="black", drop=-.15, haplotypes="blue")
genemodel.plot(rsz.lower.df, start=1, bpstop=678, orientation="reverse", xaxis=F)

# TFB5: no pred for this so all CDS; only for scheme purposes

tfb5.lower <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_12324|Cycle_0.2"),]
tfb5.upper <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_12324|Cycle_0.1"),]

length(Transabyss[[which(names(Transabyss) == "bcc_12324|Cycle_0.2")]])
length(Transabyss[[which(names(Transabyss) == "bcc_12324|Cycle_0.1")]])

tfb5.lower.df <- data.frame(type = c("3' utr", "coding_region"),
                           coordinates = c("0-0","1-759"))
tfb5.upper.df <- data.frame(type = c("3' utr","coding_region"),
                           coordinates = c("0-0","1-854"))

length(VariableSeq[[which(names(VariableSeq) == "bcc_12324|Cycle_0")]])
bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == "bcc_12324|Cycle_0.1")]



par(mfrow=c(2,1))
genemodel.plot(tfb5.upper.df, start=1, bpstop=854, orientation="reverse", xaxis=T)
genemodel.plot(tfb5.lower.df, start=1, bpstop=759, orientation="reverse", xaxis=T)

par(mfrow=c(2,1))
genemodel.plot(tfb5.upper.df, start=1, bpstop=854, orientation="reverse", xaxis=T)
mutation.plot(449, 450, text="+95 nt", col="black", drop=-.15, haplotypes="red")
mutation.plot(808+50, 809+50, text="PTC", col="black", drop=-.15, haplotypes="blue")
genemodel.plot(tfb5.lower.df, start=1, bpstop=678, orientation="reverse", xaxis=F)

# GlyRich

glyrich.lower <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_40022|Cycle_1.2"),]
glyrich.upper <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_40022|Cycle_1.1"),]

length(Transabyss[[which(names(Transabyss) == "bcc_40022|Cycle_1.2")]])
length(Transabyss[[which(names(Transabyss) == "bcc_40022|Cycle_1.1")]])

glyrich.lower.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                           coordinates = c("1-1482", "1483-2118", "2118-2370"))
glyrich.upper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                           coordinates = c("1-1727", "1728-2363", "2364-2575"))

length(VariableSeq[[which(names(VariableSeq) == "bcc_40022|Cycle_1")]])
bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == "bcc_40022|Cycle_1.2")]

glyrich.lowerupper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                                coordinates = c("1-1601", "1602-2237", "2238-2489"))

x[x$Isoforms == "bcc_40022|Cycle_1.2",]

par(mfrow=c(2,1))
genemodel.plot(glyrich.upper.df, start=1, bpstop=2575, orientation="reverse", xaxis=T)
genemodel.plot(glyrich.lower.df, start=1, bpstop=2370, orientation="reverse", xaxis=F)

par(mfrow=c(2,1))
genemodel.plot(glyrich.lowerupper.df, start=1, bpstop=2489, orientation="reverse", xaxis=T)
mutation.plot(1062, 1063, text="+119 nt", col="black", drop=-.15, haplotypes="red")
mutation.plot(808+50, 809+50, text="PTC", col="black", drop=-.15, haplotypes="blue")
genemodel.plot(rsz.lower.df, start=1, bpstop=678, orientation="reverse", xaxis=F)

# UbiqE2

e2.lower <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_15225|Cycle_0.2"),]
e2.upper <- CodAn.gtf[which(CodAn.gtf$V1 %in% "bcc_15225|Cycle_0.1"),]

length(Transabyss[[which(names(Transabyss) == "bcc_15225|Cycle_0.2")]])
length(Transabyss[[which(names(Transabyss) == "bcc_15225|Cycle_0.1")]])

e2.lower.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                               coordinates = c("1-330", "331-789", "790-1065"))
e2.upper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                               coordinates = c("1-657", "658-1119", "1120-1395"))

length(VariableSeq[[which(names(VariableSeq) == "bcc_15225|Cycle_0")]])
bam.coords[[1]]$pos[which(bam.coords[[1]]$qname == "bcc_15225|Cycle_0.2")]

e2.lowerupper.df <- data.frame(type = c("3' utr", "coding_region", "5' utr"),
                                    coordinates = c("1-330", "331-877", "878-1153"))

x[x$Isoforms == "bcc_15225|Cycle_0.2",]

par(mfrow=c(2,1))
genemodel.plot(e2.upper.df, start=1, bpstop=1395, orientation="reverse", xaxis=T)
genemodel.plot(e2.lower.df, start=1, bpstop=1065, orientation="reverse", xaxis=T)

par(mfrow=c(2,1))
genemodel.plot(e2.lowerupper.df, start=1, bpstop=1153, orientation="reverse", xaxis=T)
mutation.plot(695, 696, text="+81 nt", col="black", drop=-.15, haplotypes="red")
mutation.plot(695+40, 696+40, text="PTC", col="black", drop=-.15, haplotypes="blue")
genemodel.plot(e2.lower.df, start=1, bpstop=678, orientation="reverse", xaxis=F)


## Isoforms -- Transcript co-relation (fifth)

blastp.res <- read.table("D:/AA_PProject1_TFM_final/AA_GIT_WorkingDirectory/Transabyss_Full/blast_results/blast_result.tabular")
prot_CT1 <- read.table("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/5.Validation/Insilico/TPCompara/ProteinsLimma_CT1.txt", header = T)
prot_CT3 <- read.table("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/5.Validation/Insilico/TPCompara/ProteinsLimma_CT3.txt", header = T)
prot_T1T3 <- read.table("D:/AA_GIT_myproject/SplicingAnalysis_PinusRadiata/Data/5.Validation/Insilico/TPCompara/ProteinsLimma_T1T3.txt", header = T)
isoforms_CT1 <- DESeq2_all$CT1
isoforms_CT3 <- DESeq2_all$CT3
isoforms_T1T3 <- DESeq2_all$T1T3
All_prot <- read.delim("clipboard", header = T)# from mofa_inputs
All_isoforms <- read.delim("clipboard", header = T)# from mofa_inputs

# scatterplot: volcano

# CT1

targets <- blastp.res[which(blastp.res$V3 > 80 & blastp.res$V11 < 10^-5),c(1,2)]
targets <- targets[targets$V1 %in% isoforms_CT1$KissID,]
targets$Isoforms <- isoforms_CT1$log2FoldChange[which(isoforms_CT1$KissID %in% as.character(targets$V1))]
targets$Proteins <- NA
delete <- list()
for(i in 1:nrow(targets)){
  if(length(which(rownames(prot_CT1) == as.character(targets$V2[i]))) == 0){
    delete[[i]] <- i
    next
  }
  targets$Proteins[i] <- prot_CT1$logFC[which(rownames(prot_CT1) == as.character(targets$V2[i]))]
}
targets <- targets[-unlist(delete, use.names = F),]
targets$group <- "NoRelation"
targets$group[which(abs(targets$Isoforms) > 0.5 & abs(targets$Proteins) > 0.5)] <- "Medium"
targets$group[which(abs(targets$Isoforms) > 1 & abs(targets$Proteins) > 1)] <- "High"
targets$group[which(abs(targets$Isoforms) > 1.8 & abs(targets$Proteins) > 1.8)] <- "VeryHigh"
write.table(targets, file = "TPVolcano_CT1.txt",sep = "\t", col.names = T)
pal <- c( "lightcoral", "lightgoldenrod1", "lightsteelblue", "lightskyblue4")

p1 <- plot_ly(data = targets, x = targets$Isoforms, y = targets$Proteins, text = paste0(targets$V1, ";", targets$V2), mode = "markers", color = targets$group , colors = pal, xaxis = "Isoforms FoldChange[log2]",  yaxis = "Proteins FoldChange[log2]") %>% 
  layout(title = "Control vs T1")

orca(p1, file = "TPVolcano_CT1.pdf")

# CT3

targets <- blastp.res[which(blastp.res$V3 > 80 & blastp.res$V11 < 10^-5),c(1,2)]
targets <- targets[targets$V1 %in% isoforms_CT3$KissID,]
targets$Isoforms <- isoforms_CT3$log2FoldChange[which(isoforms_CT3$KissID %in% as.character(targets$V1))]
targets$Proteins <- NA
delete <- list()
for(i in 1:nrow(targets)){
  if(length(which(rownames(prot_CT3) == as.character(targets$V2[i]))) == 0){
    delete[[i]] <- i
    next
  }
  targets$Proteins[i] <- prot_CT3$logFC[which(rownames(prot_CT3) == as.character(targets$V2[i]))]
}
targets <- targets[-unlist(delete, use.names = F),]
targets$group <- "NoRelation"
targets$group[which(abs(targets$Isoforms) > 0.5 & abs(targets$Proteins) > 0.5)] <- "Medium"
targets$group[which(abs(targets$Isoforms) > 1 & abs(targets$Proteins) > 1)] <- "High"
targets$group[which(abs(targets$Isoforms) > 1.8 & abs(targets$Proteins) > 1.8)] <- "VeryHigh"
write.table(targets, file = "TPVolcano_CT3.txt",sep = "\t", col.names = T)
pal <- c( "lightcoral", "lightgoldenrod1", "lightsteelblue")

p2 <- plot_ly(data = targets, x = targets$Isoforms, y = targets$Proteins, text = paste0(targets$V1, ";", targets$V2), mode = "markers", color = targets$group , colors = pal, xaxis = "Isoforms FoldChange[log2]",  yaxis = "Proteins FoldChange[log2]") %>% 
  layout(title = "Control vs T3")

orca(p2, file = "TPVolcano_CT3.pdf")

# T1T3

targets <- blastp.res[which(blastp.res$V3 > 80 & blastp.res$V11 < 10^-5),c(1,2)]
targets <- targets[targets$V1 %in% isoforms_T1T3$KissID,]
targets$Isoforms <- isoforms_T1T3$log2FoldChange[which(isoforms_T1T3$KissID %in% as.character(targets$V1))]
targets$Proteins <- NA
delete <- list()
for(i in 1:nrow(targets)){
  if(length(which(rownames(prot_T1T3) == as.character(targets$V2[i]))) == 0){
    delete[[i]] <- i
    next
  }
  targets$Proteins[i] <- prot_T1T3$logFC[which(rownames(prot_T1T3) == as.character(targets$V2[i]))]
}
targets <- targets[-unlist(delete, use.names = F),]
targets$group <- "NoRelation"
targets$group[which(abs(targets$Isoforms) > 0.5 & abs(targets$Proteins) > 0.5)] <- "Medium"
targets$group[which(abs(targets$Isoforms) > 1 & abs(targets$Proteins) > 1)] <- "High"
targets$group[which(abs(targets$Isoforms) > 1.8 & abs(targets$Proteins) > 1.8)] <- "VeryHigh"
write.table(targets, file = "TPVolcano_T1T3.txt",sep = "\t", col.names = T)
pal <- c("lightsteelblue")

p3 <- plot_ly(data = targets, x = targets$Isoforms, y = targets$Proteins, text = paste0(targets$V1, ";", targets$V2), mode = "markers", color = targets$group , colors = pal, xaxis = "Isoforms FoldChange[log2]",  yaxis = "Proteins FoldChange[log2]") %>% 
  layout(title = "T1 vs T3")

orca(p3, file = "TPVolcano_T1T3.pdf")

# UMAP and project

# isoforms training predicts proteins

targets <- blastp.res[which(blastp.res$V3 > 80 & blastp.res$V11 < 10^-5),c(1,2)]
isoforms.matrix <- All_isoforms[which(All_isoforms$KissID %in% targets$V1),]
rownames(isoforms.matrix) <- isoforms.matrix$KissID
isoforms.matrix <- t(isoforms.matrix[,-1])
proteins.matrix <- as.data.frame(matrix(NA, nrow(isoforms.matrix), ncol(isoforms.matrix)))
colnames(proteins.matrix) <- colnames(isoforms.matrix)
rownames(proteins.matrix) <- rownames(isoforms.matrix)
for(i in 1:ncol(proteins.matrix)){
  prot.search <- as.character(targets$V2[which(targets$V1 == colnames(proteins.matrix)[i])])
  proteins.matrix[,i] <- as.numeric(All_prot[which(All_prot$ID == prot.search),-1])
}

isoforms.matrix <- isoforms.matrix[,-which(is.na(proteins.matrix[1,]))]
proteins.matrix <- proteins.matrix[,-which(is.na(proteins.matrix[1,]))]
set.seed(404)
iso_umap <- uwot::umap(scale(isoforms.matrix), n_components = 2, n_threads = 6, verbose = T, init = "spca", n_neighbors = 9, ret_model = T)
iso_umap

iso_prot_umap <- uwot::umap_transform(scale(proteins.matrix), iso_umap, verbose = TRUE, n_threads = 6)

fumap <- as.data.frame(rbind(iso_umap$embedding, iso_prot_umap))
fumap$Level <- c(rep("Isoforms", 9), rep("Proteins", 9))
fumap$Treatment <-  c(rep(c(rep("C",3), rep("T1",3), rep("T3",3)), 2))
pal_umap <- c("#00A08799", "#3C548899", "#E64B3599")
pal_umap <- setNames(pal_umap, c("C", "T1", "T3"))
colnames(fumap)[1:2] <- c("UMAP1", "UMAP2")

p4 <- plot_ly(data = fumap, x = ~UMAP1, y = ~UMAP2, color = ~Treatment, colors = pal_umap, symbol = ~Level, symbols = c('circle','x','o'), marker = list(size=10)) %>%
  layout(title = "UMAP Isoforms~Proteins")

orca(p4, file = "UMAP_isoforsmpredictsproteins.pdf")
