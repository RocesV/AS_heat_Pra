########################################################################################
############################## Gotubro Aux.Functions  ##################################
########################################################################################

gotubro.parser <- function(dammit.output = NULL, dammit.namemap = NULL, interpro.output = NULL, interpro = F){
  ###################################################################################################
  # gotubro.parser == select the best hit for each DB by sequence length, significance and redundancy
  ###################################################################################################
  # ARGS:
  # dammit.output   == dammit annotation path                   (character)
  # dammit.namemap  == dammit namemap output path               (character)
  # interpro.output == interpro annotation path                 (character)
  # interpro        == if TRUE, interpro parser is executed     (logical)
  ###################################################################################################
  # OUTPUT:
  # Data.frame: rows each sequence in the query of dammit/IP, columns best hit for each DB. Query or  
  # DB without hit are filled with NAs. If interpro is TRUE, the columns are DBs implemented in the 
  # functional gene ontology cross-over (hamap, panther, smart, prints, prosite, interpro, pirsf).
  ###################################################################################################
  # HELP:
  # interpro.output should be in chunks. All the files in the same directory should come from the 
  # same assembly. First log annotated refers to a query that has a hit for any of the DBs tested.
  ###################################################################################################
  
  #### 0 Load pkgs
  
  suppressPackageStartupMessages(library(cowsay)) 
  suppressPackageStartupMessages(library(multicolor)) 
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(progress))
  
  #### 1 Conditionals -- full stop
  
  if(is.null(dammit.namemap)){stop(say(what = "\n dammit.namemap is required \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))}
  if(isTRUE(interpro) & is.null(interpro.output)){stop(say(what = "\n interpro.output is required if interpro is TRUE \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))}
  if(isFALSE(interpro) & is.null(dammit.output)){stop(say(what = "\n dammit.output is required if interpro is FALSE \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))}
  
  ################################################################
  ######################## INTERPRO ##############################
  ################################################################
  
if(isTRUE(interpro)){
  
  #### 2 Introduction-Import-Data.frame
  
  ### Introduction
  
  say(what = paste0(Sys.time(), " || Welcome to Got u Bro! Victor´s auxiliar functions to process interpro outputs" ), by = "rabbit", what_color = "salmon", by_color = "salmon")
  
  ### Import: names and chunks
  
  cat("\n Importing Chunks and NameMap... \n")
  
  dammit.names <- read.csv(dammit.namemap, header = T, sep = ",")
  interpro.tables <- list()
  for(i in list.files(interpro.output)){
    interpro.tables[[i]] <- read.csv(paste0(interpro.output,i), header = F, sep = "\t")
    if(ncol(interpro.tables[[i]]) == 11){
      delete.rows <- list()
      add.columns <- list()
      for(j in 1:nrow(interpro.tables[[i]])){
        if(j %in% unlist(delete.rows,use.names = F)){next}
        condition <- grep("bcc_", interpro.tables[[i]][j+1,1])
        if(length(condition) == 0){
          delete.rows[[j]] <- j+1
          add.columns[[j]] <- interpro.tables[[i]][j+1,1:4]
        }else if(length(condition) == 1){next}
      }
      delete.rows <- unlist(delete.rows, use.names = F)
      add.columns <- do.call(rbind,add.columns)
      add.columns$V1 <- as.character(add.columns$V1)
      add.columns$V2 <- as.character(add.columns$V2)
      add.columns$V3 <- as.character(add.columns$V3)
      add.columns$V4 <- as.character(add.columns$V4)
      tocbind <- as.data.frame(matrix(data = NA, nrow = nrow(interpro.tables[[i]]), ncol = 4))
      tocbind[delete.rows -1,1] <- add.columns[,1]
      tocbind[delete.rows -1,2] <- add.columns[,2]
      tocbind[delete.rows -1,3] <- add.columns[,3]
      tocbind[delete.rows -1,4] <- add.columns[,4]
      interpro.tables[[i]] <- cbind(interpro.tables[[i]], tocbind)
      interpro.tables[[i]] <- interpro.tables[[i]][-delete.rows,]
      if(ncol(interpro.tables[[i]]) != 15){stop(cat("\n Something go wrong boy ... \n"))}
    }
    colnames(interpro.tables[[i]]) <- c("CustomID", "md5sum", "X", "DB", "DBid", "DBname", "start", "end", "q-value", "Y", "Z", "IPid", "IPname", "IPgo", "Pathway")
    interpro.tables[[i]][,1] <- do.call(rbind,strsplit(as.character(interpro.tables[[i]][,1]), split = "_orf", fixed = T))[,1]
  }
  interpro.tables <- lapply(interpro.tables, function(x){
    jj <- apply(x, 2, as.character)
    jj
  })
  interpro.tables <- do.call(rbind, interpro.tables)
  identified <- length(which(as.character(dammit.names[,1]) %in% interpro.tables[,1] == TRUE))
  percentage.identification <- round(identified/nrow(dammit.names)* 100, 2)
  cat(paste0("\n ", percentage.identification, " % sequences have been identified (", identified, ") \n"))
  
  #### 3 Best hit selection
  
  cat("\n Hit selection for each database ... (hamap, panther, smart, interpro, prints, prosite, pirsf) \n")
  
  total <- nrow(dammit.names)
  final.df <- as.data.frame(matrix(NA, nrow = total, ncol = 8))
  colnames(final.df) <- c("KissID", "Hamap", "PANTHER", "SMART", "PRINTS", "ProSite", "PIRSF", "InterPro")
  interpro.tables <- interpro.tables[,-c(2,3,10,11)]
  dammit.names[,1] <- as.character(dammit.names[,1])
  final.df[,1] <- dammit.names[,1]
  
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for(i in 1:total){ # for 1:total
    setTxtProgressBar(pb, i)
    condition1 <- which(interpro.tables[,1] == dammit.names[i,1])
    if(length(condition1) == 0){ # no hit for this transcript -- next
      next
    }else if(length(condition1) == 1){ # only one hit for this transcript
      for(j in colnames(final.df)[-1]){
        condition4 <- grep(j, interpro.tables[condition1,2])
        if(length(condition4) == 0){ # hit for a DB not taked into account
          next
        }else if(length(condition4) > 0){ # hit for a DB taked into account
          if(j == "InterPro"){
            final.df[i,j] <- interpro.tables[condition1, 8]
          }else if(j != "InterPro"){
             final.df[i,j] <- interpro.tables[condition1, 3]
           }
        }
      }
    }else if(length(condition1) > 1){
      for(j in colnames(final.df)[-1]){ #for colnames DB
        if(j == "InterPro"){
          IP.specific <- grep("IPR",interpro.tables[condition1,8])
          if(length(IP.specific) == 0){ # No hit for InterPro
            next
          }else if(length(IP.specific) == 1){ # One hit for InterPro
            final.df[i,8] <- interpro.tables[condition1,8][IP.specific]
          }else if(length(IP.specific) > 1){ # Multiple hits for InterPro: several domains are possible inside the same protein/transcript
            # no dups, collapse by "|" as dammit
            hits <- interpro.tables[condition1,8][IP.specific]
            final.df[i,8] <- paste0(hits[!duplicated(hits)], collapse = " | ")
          }
        }else if(j != "InterPro"){
          DB.specific <- grep(j, interpro.tables[condition1,2])
          if(length(DB.specific) == 0){ # No hit for j DB
            next
          }else if(length(DB.specific) == 1){ # Unique hit for j DB
            final.df[i,j] <- interpro.tables[condition1,3][DB.specific]
            if(j == "Hamap"){
                final.df[i,j] <- gsub("_B", "",final.df[i,j], fixed = T)  
            }
          }else if(length(DB.specific) > 1){ # Multiple hits for j DB
            hits <- interpro.tables[condition1,3][DB.specific]
            no.dups <- hits[!duplicated(hits)]
            if(j == "Hamap"){
              no.dups <- gsub("_B", "",no.dups, fixed = T)
            }
            if(j == "PANTHER"){
              condition2 <- grep(":SF", no.dups)
              if(length(condition2) == 0){
                no.dups <- no.dups
              }else if(length(condition2) > 0){
                no.dups <- no.dups[condition2]
              }
            } 
           final.df[i,j] <- paste0(no.dups, collapse = " | ") 
          } # Multiples hits for j DB
        } # InterPro | != InterPro
      } # for colnames DBs
      
    } # condition1 conditional
  } # for 1:total
  close(pb)
  cat(paste0("\n DISCLAIMER: annotated percentage take into account all DBs from InterPro and not only the ones considered for downstream applications \n"))
  say(what = paste0(Sys.time(), " || Thanks for using! <3" ), by = "rabbit", what_color = "salmon", by_color = "salmon")
  return(final.df)
  #### CHECKKKK SOME is with interprotables to check | SOMETHING go wrong check INTERPRO
  
  
  
} # INTERPRO-DAMMIT conditional  
  
  ################################################################
  ######################### DAMMIT ###############################
  ################################################################
  
if(isFALSE(interpro)){
  
  #### 2 Introduction-Import-Data.frame
  
  ### Introduction
  
  say(what = paste0(Sys.time(), " || Welcome to Got u Bro! Victor´s auxiliar functions to process dammit outputs" ), by = "rabbit", what_color = "salmon", by_color = "salmon")
  
  ### Import: GFF3 issues
  
  cat("\n Importing GFF3 and NameMap... \n")
  
  dammit.GFF3 <- read.delim(dammit.output, header = F, sep = "\t")
  dammit.GFF3 <- dammit.GFF3[-1,]
  dammit.names <- read.csv(dammit.namemap, header = T, sep = ",")
  dammit.GFF3.noTransdecoder <- dammit.GFF3[which(dammit.GFF3$V2 != "transdecoder"),]
  
  identified <- length(dammit.names$renamed[which(dammit.names$renamed %in% dammit.GFF3.noTransdecoder$V1)])
  total <- nrow(dammit.names)
  percentage.identification <- round(identified / total * 100,2)
  dammit.GFF3.noTransdecoder$V9 <- as.character(dammit.GFF3.noTransdecoder$V9)
  
  cat(paste0("\n ", percentage.identification, " % sequences have been identified (", identified, ") \n"))
  
  #### 3 Best hit selection
  
  cat("\n Best hit selection for each database ... \n")
  
  final.df <- as.data.frame(matrix(NA, nrow(dammit.names), ncol(dammit.names) + 6))
  colnames(final.df) <- c("KissID", "dammitID", "uniref", "=sprot","OrthoDB", "Pfam:", "Rfam:", "Epiomides")
  final.df[,c(1,2)] <- dammit.names
  rm(dammit.names)
  rm(dammit.GFF3)
  
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for(i in 1:total){ # for loop all query sequences
    setTxtProgressBar(pb, i)
    condition1 <- which(as.character(dammit.GFF3.noTransdecoder$V1) == as.character(final.df$dammitID[i]))
    
    if(length(condition1) == 0){ # No hit for this transcript
      next
    }else if(length(condition1) == 1){ # Only one hit for this transcript
      single.transcript <- dammit.GFF3.noTransdecoder[condition1,]
      description <- do.call(rbind,strsplit(as.character(single.transcript$V9), split = ";"))
      for(j in colnames(final.df)[-c(1,2)]){ # search DB that contains the only hit for this transcript
        DB.classification <- grep(j, description)
        if(length(DB.classification) == 0){ # j DB does not contain the hit
          next
        }else if(length(DB.classification) > 0){ # j DB contains the hit
          if(j != "Rfam:" & j != "Pfam:"){
            final.df[i,j] <- description[grep("Name=", description)]  
          }else if(j == "Pfam:" & j != "Rfam:"){
            final.df[i,j] <- description[grep("Dbxref=", description)]
          }else if(j != "Pfam:" & j == "Rfam:"){
            final.df[i,j] <- description[grep("Dbxref=", description)] 
          }}
        }
    }else if(length(condition1) > 1){ # Multiple hits for this transcript
      # for each DB search best hit
      for(j in colnames(final.df)[-c(1,2)]){ # search DB that contains the hits for this transcript
        DB.specific <- grep(j, dammit.GFF3.noTransdecoder$V9[condition1])
        if(length(DB.specific) == 0){ # j DB without any hit for this transcript
          next
        }else if(length(DB.specific) == 1){ # j DB with only one hit for this transcript
          description <- do.call(rbind, strsplit(dammit.GFF3.noTransdecoder$V9[condition1][DB.specific], split = ";"))   
          if(j != "Rfam:" & j != "Pfam:"){
            final.df[i,j] <- description[grep("Name=", description)]  
          }else if(j == "Pfam:" & j != "Rfam:"){
            final.df[i,j] <- description[grep("Dbxref=", description)]
          }else if(j != "Pfam:" & j == "Rfam:"){
            final.df[i,j] <- description[grep("Dbxref=", description)] 
          }
        }else if(length(DB.specific) > 1){ # j DB with multiple hits for this transcript
          # apply selection criteria: 1) sequence length, 2) significance 3) most common result || Except Pfam
          DB.specific.multiple <- dammit.GFF3.noTransdecoder[condition1[DB.specific],]
          if(j == "Pfam:" | j == "Rfam:"){ # avoid filtering because could be several protein domains in the same sequence
            if(j == "Pfam:"){
              check <- do.call(rbind, strsplit(as.character(DB.specific.multiple$V9), split = ";", fixed = T))[,8]  
            }else if(j == "Rfam:"){
              check <- do.call(rbind, strsplit(as.character(DB.specific.multiple$V9), split = ";", fixed = T))[,5]  
            }
            if(length(which(duplicated(check) == TRUE)) == 0){ # no dups domains
              final.df[i,j] <- paste0(check, collapse = " | ")
            }else if(length(which(duplicated(check) == TRUE)) > 0){ # dups domains, dont retain redundant information
              final.df[i,j] <- paste0(check[!duplicated(check)], collapse = " | ")
            }
          }else if(j != "Pfam:" & j != "Rfam:"){ # need filtering -- unique hit
            # 1) Sequence length & 2) significance
            filter12 <- DB.specific.multiple[which((DB.specific.multiple$V5 - DB.specific.multiple$V4) == max(DB.specific.multiple$V5 - DB.specific.multiple$V4) & as.numeric(as.character(DB.specific.multiple$V6)) == min(as.numeric(as.character(DB.specific.multiple$V6)))),9]
            if(length(filter12) == 0){ # not hit that  matches both conditions
              filter1 <- DB.specific.multiple[which((DB.specific.multiple$V5 - DB.specific.multiple$V4) == max(DB.specific.multiple$V5 - DB.specific.multiple$V4)),9]
              filter2 <- DB.specific.multiple[which(as.numeric(as.character(DB.specific.multiple$V6)) == min(as.numeric(as.character(DB.specific.multiple$V6)))),9]
              filter12.mod <- data.frame(filters = c(filter1, filter2))
              filter12.mod$filters <- as.character(filter12.mod$filters)
              description <- do.call(rbind, strsplit(filter12.mod$filters, split = ";", fixed = T))
              if(length(grep("|",filter12.mod[,1], fixed = T)) == 0 | length(grep("|",filter12.mod[,1], fixed = T)) == 1){ # no hit is annotated so take the first by default or only one hit is annotated
                  if(length(grep("|",filter12.mod[,1], fixed = T)) == 1){
                    select <- grep("|",filter12.mod[,1], fixed = T)
                    final.df[i,j] <- description[grep("Name=", description)][select] 
                  }else if(length(grep("|",filter12.mod[,1], fixed = T)) == 0){
                    final.df[i,j] <- description[grep("Name=", description)][1]
                  } 
              }else if(length(grep("|",filter12.mod[,1], fixed = T)) > 1){ # several hits: priority hits annotated and described
                  # conditional of description with only one or more hits || one just add || several take the first
                  if(nrow(description) == 1){
                   final.df[i,j] <- description[grep("Name=", description)] 
                  }else if(nrow(description) > 1){ ### Conditional search OrthoID for epiomides
                    if(j == "Epiomides"){
                      Ortho.1KP <- grep("OrthoID", description[,2])
                      if(length(Ortho.1KP) == 0){
                        final.df[i,j] <- description[grep("Name=", description)][1]
                      }else if(length(Ortho.1KP) == 1){
                        final.df[i,j] <- description[Ortho.1KP,2]
                      }else if(length(Ortho.1KP) > 1){ # several with OrthoIDs, select best common name
                        filter3 <- which(as.numeric(as.factor(do.call(rbind,strsplit(description[Ortho.1KP,2], "|", fixed = T))[,2])) == which.max(table(as.numeric(as.factor(do.call(rbind,strsplit(description[Ortho.1KP,2], "|", fixed = T))[,2]))))) 
                        final.df[i,j] <- description[filter3[1],2] 
                      }
                    }else if(j != "Epiomides"){
                      final.df[i,j] <- description[grep("Name=", description)][1]
                    }
                  }else if(nrow(description) == 0){stop("\n Something go wrong with hits in filter12: CHECK !! \n")}
              }
            }else if(length(filter12) == 1){ # hit  that matches both conditions
              description <- do.call(rbind, strsplit(filter12, split = ";"))   
              final.df[i,j] <- description[grep("Name=", description)]  
            }else if(length(filter12) > 1){ # several hits that matches both conditions || 3) most common result
              description <- do.call(rbind, strsplit(filter12, split = ";"))
              if(j == "Epiomides"){ # most common Annotation name that has OrthoID: Priority 1Kp Orthos >>> PLAZA
                Ortho.1KP <- grep("OrthoID", description[,2])
                if(length(Ortho.1KP) == 0){ # multiple hits but no one with OrthoID (PLAZA - not classified 1KP)
                  if(length(grep("|", description[,2], fixed = T)) == 0){ # need conditional with | maybe a hit from PLAZA without annotation
                    filter3 <- which(as.numeric(as.factor(description[,2])) == which.max(table(as.numeric(as.factor(description[,2]))))) 
                    final.df[i,j] <- description[filter3[1],2] 
                  }else if(length(grep("|", description[,2], fixed = T)) != 0){
                    filter3 <- which(as.numeric(as.factor(do.call(rbind,strsplit(description[,2], "|", fixed = T))[,2])) == which.max(table(as.numeric(as.factor(do.call(rbind,strsplit(description[,2], "|", fixed = T))[,2]))))) 
                    final.df[i,j] <- description[filter3[1],2]
                  }
                }else if(length(Ortho.1KP) > 1){ # multiple hits and OrthoID (1KP)
                  # dont need | conditional because if i have OrthoID it is supposed to be there
                  filter3 <- which(as.numeric(as.factor(do.call(rbind,strsplit(description[Ortho.1KP,2], "|", fixed = T))[,2])) == which.max(table(as.numeric(as.factor(do.call(rbind,strsplit(description[Ortho.1KP,2], "|", fixed = T))[,2]))))) 
                  final.df[i,j] <- description[filter3[1],2] 
                }else if(length(Ortho.1KP) == 1){ # only one hit with OrthoID (1KP)
                  # dont need | conditional because if i have OrthoID it is supposed to be there
                  final.df[i,j] <- description[Ortho.1KP,2]
                }
              }else if(j == "uniref"){
                stop("\n More than 1 Uniref90 hit for this transcript: CHECK !! \n")
              }else if(j == "=sprot"){
                stop("\n More than 1 sprot hit for this transcript: CHECK !! \n")
              }else if(j == "OrthoDB"){
                stop("\n More than 1 OrthoDB hit for this transcript: CHECK !! \n")
              } # most common Annotation name
            } # 3) most common result
           } # Pfam avoid filtering
          } # j DB with multiple hits for this transcript 
        } # search DB than contains the hits for this transcript
      } # Multiple hits for this transcript 
    } # for loop all query sequences
  close(pb)
  say(what = paste0(Sys.time(), " || Thanks for using! <3" ), by = "rabbit", what_color = "salmon", by_color = "salmon")
  return(final.df)
 } # INTERPRO-DAMMIT conditional
} # function end

gotubro.crossrefs <- function(gotubro.parser.dammit, gotubro.parser.ip , SQLite3db = "F:/SQLite3/gotubro", quickGO = TRUE, mode = c("reduce", "complete", "normal"), return.merged = T){
  ##################################################################################################################################################
  # gotubro.crossxrefs == for the best hit for each DB maps the gene ontology (GO) using a comprehensive approach taking advantage of SQL language
  ##################################################################################################################################################
  # ARGS:
  # goutbro.parser.dammit == output from goutbro.parser with dammit args                                                   (character)
  # goutbro.parser.ip     == output from goutbro.parser with interpro args                                                 (character)
  # SQLite3db             == SQLite3 database path                                                                         (character)
  # quickGO               == if TRUE, functional annotations are completed using quickGO                                   (logical)
  # mode                  == reduce, complete or normal                                                                    (character)
  # return.merged         == if TRUE, the gotubro.parser merged table is merged with the final GO table                    (logical)
  ##################################################################################################################################################
  # OUTPUT:
  # Data.frame: rows each sequence in the query of gotubro.parser; if mode is reduced, only one column with non-duplicated GOs; if mode is complete, 
  # multiple columns with GOs, IDs and desc for each DB; if mode is normal, multiple columns with GOs for each DB. Query or DB without hit are 
  # filled with NAs.
  ##################################################################################################################################################
  # HELP:
  # gotubro.parser.dammit or gotubro.parser.ip can be used as query alone or together. If you use both, tables are merged. quickGO is advisable always 
  # because nowdays (07/09/2020) it is the most complete DB for GOs. Epiomides and OrthoDB are not considered in the functional annotation
  # (bad resources), only useful as identification. Mode arg can be a vector with several modes
  ##################################################################################################################################################
  
  #### 0 Load pkgs
  
  suppressPackageStartupMessages(library(dplyr)) 
  suppressPackageStartupMessages(library(dbplyr)) 
  suppressPackageStartupMessages(library(RSQLite)) 
  suppressPackageStartupMessages(library(cowsay))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(progress))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(multicolor))
  
  #### 1 Conditionals -- full stop
  
  if(is.null(gotubro.parser.dammit) & is.null(gotubro.parser.ip)){stop(say(what = "\n At least one gotubro.parser output is needed \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))}
  if(is.null(mode) | length(which((mode %in% c("reduce", "complete", "normal")) == TRUE)) == 0){stop(say(what = "\n Select one of the possible mode options \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))}
  if(is.null(quickGO) | is.null(return.merged)){stop(say(what = "\n Set quickGO and return.merged args \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))}
  
  #### 2 Introduction-Import-Data.frames
  
  ### Introduction
  
  say(what = paste0(Sys.time(), " || Welcome to Got u Bro! Victor´s auxiliar functions to process dammit/interpro outputs" ), by = "rabbit", what_color = "salmon", by_color = "salmon")
  
  ### Import Data.frames
  
  cat(paste0(Sys.time(), " || Importing gotubro.parser outputs ... \n"))
  
  if(!is.null(gotubro.parser.dammit) & is.null(gotubro.parser.ip)){
   gotubro.parser.id <- read.delim(gotubro.parser.dammit, sep = "\t", header = T)
   gotubro.parser <- gotubro.parser.id[,-c(2,5,8)]
   gotubro.parser.id <- gotubro.parser.id[,-2]
  }else if(is.null(gotubro.parser.dammit) & !is.null(gotubro.parser.ip)){
    gotubro.parser <- read.delim(gotubro.parser.ip, sep = "\t", header = T) 
  }else if(!is.null(gotubro.parser.dammit) & !is.null(gotubro.parser.ip)){
    gotubro.parser.d <- read.delim(gotubro.parser.dammit, sep = "\t", header = T)
    gotubro.parser.i <- read.delim(gotubro.parser.ip, sep = "\t", header = T)
    if(length(which((gotubro.parser.d$KissID == gotubro.parser.i$KissID) == TRUE)) == nrow(gotubro.parser.d)){
      gotubro.parser <- cbind(gotubro.parser.d[,-c(2,5,8)], gotubro.parser.i[,-1])
      gotubro.parser.id <- cbind(gotubro.parser.d[,-2], gotubro.parser.i[,-1])
    }else if(length(which((gotubro.parser.d$KissID == gotubro.parser.i$KissID) == TRUE)) != nrow(gotubro.parser.d)){
      stop(say(what = "\n Something go wrong! :( Ask Victor: fernandezvictor@uniovi.es \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))
      }
  }
  
  gotubro.parser <- apply(gotubro.parser, 2, FUN = as.character)
  gotubro.parser.id <- apply(gotubro.parser.id, 2, FUN = as.character)
  identified <- length(which((rowSums(is.na(gotubro.parser.id[,-1])) != ncol(gotubro.parser.id[,-1])) == TRUE))
  percentage.identification <- round(identified/nrow(gotubro.parser.id) * 100, 2)
  cat(paste0("\n ", percentage.identification, " % sequences have been identified (", identified, ") \n"))
  
  #### 3 Connecting to gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Connecting to goutbro SQLite3 DB ... \n"))
  
  gotubro.sqlite <- DBI::dbConnect(RSQLite::SQLite(), SQLite3db)
  # check .schema -- src_dbi(gotubro.sqlite)
  SQLquery.results <- list()
  
  #### 4 Querying the SQLite3 DB (dplyr sintax)
  
  cat(paste0(Sys.time(), " || Querying goutbro SQLite3 DB ... \n"))
  
  final.df <- as.data.frame(matrix(NA, nrow(gotubro.parser), ncol(gotubro.parser)))
  colnames(final.df) <- colnames(gotubro.parser)
  final.df[,1] <- gotubro.parser[,1]
  
  ############################################################################################################################
  
  ### SMART Domains
  
  cat(paste0(Sys.time(), " || Tunning SMART DB subset ... \n"))
  
  ## gotubro.parser.subset.smart table
  
  gotubro.parser.subset.smart <- as.data.frame(gotubro.parser[,c(1,8)])
  colnames(gotubro.parser.subset.smart) <- c("KissID", "ID")
  complex.rows <- grep(" | ", gotubro.parser.subset.smart[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.smart$ID[grep(" | ", gotubro.parser.subset.smart$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.smart <- apply(gotubro.parser.subset.smart, 2, as.character)
  gotubro.parser.subset.smart[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.smart <- rbind(gotubro.parser.subset.smart, data.frame(KissID = gotubro.parser.subset.smart[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.smart <- rbind(gotubro.parser.subset.smart, data.frame(KissID = gotubro.parser.subset.smart[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }  
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.smart)
  gotubro.parser.subset.smart <- tbl(gotubro.sqlite, "gotubro.parser.subset.smart")
  
  ## smart.subset table
  
  smart <- tbl(gotubro.sqlite, "smart")
  
  ## Querying | Joining
  
  SQLquery.results[['smart']] <- gotubro.parser.subset.smart %>% left_join(smart)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['smart']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting SMART results query ... \n"))
    
  SQLquery.results[['smart']] %>% collect()
  SQLquery.results[['smart']] <- as_tibble(SQLquery.results[['smart']])
  SQLquery.results[['smart.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$smart), nrow = nrow(final.df)))
  colnames(SQLquery.results$smart.def) <- colnames(SQLquery.results$smart)
  SQLquery.results[['smart.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['smart.def']])){
    hits <- which(SQLquery.results$smart$KissID == SQLquery.results$smart.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['smart.def']])){
      if(colnames(SQLquery.results$smart.def)[j] == "KissID"){ next }
      SQLquery.results$smart.def[i,j] <- paste0(as.data.frame(SQLquery.results$smart)[hits,j][which(as.data.frame(SQLquery.results$smart)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$smart)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing SMART DB subset from goutubro SQLite3 DB \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.smart")
  rm(gotubro.parser.subset.smart)
  rm(smart)
  SQLquery.results$smart <- NULL
  
  ############################################################################################################################
  
  ### Hamap Domains
  
  cat(paste0(Sys.time(), " || Tunning Hamap DB subset ... \n"))
  
  ## gotubro.parser.subset.hamap table
  
  gotubro.parser.subset.hamap <- as.data.frame(gotubro.parser[,c(1,6)])
  colnames(gotubro.parser.subset.hamap) <- c("KissID", "ID")
  complex.rows <- grep(" | ", gotubro.parser.subset.hamap[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.hamap$ID[grep(" | ", gotubro.parser.subset.hamap$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.hamap <- apply(gotubro.parser.subset.hamap, 2, as.character)
  gotubro.parser.subset.hamap[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.hamap <- rbind(gotubro.parser.subset.hamap, data.frame(KissID = gotubro.parser.subset.hamap[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.hamap <- rbind(gotubro.parser.subset.hamap, data.frame(KissID = gotubro.parser.subset.hamap[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }  
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.hamap)
  gotubro.parser.subset.hamap <- tbl(gotubro.sqlite, "gotubro.parser.subset.hamap")
  
  ## hamap.subset table
  
  hamap <- tbl(gotubro.sqlite, "hamap")
  
  ## Querying | Joining
  
  SQLquery.results[['hamap']] <- gotubro.parser.subset.hamap %>% left_join(hamap)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['hamap']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting Hamap results query ... \n"))
  
  SQLquery.results[['hamap']] %>% collect()
  SQLquery.results[['hamap']] <- as_tibble(SQLquery.results[['hamap']])
  SQLquery.results[['hamap.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$hamap), nrow = nrow(final.df)))
  colnames(SQLquery.results$hamap.def) <- colnames(SQLquery.results$hamap)
  SQLquery.results[['hamap.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['hamap.def']])){
    hits <- which(SQLquery.results$hamap$KissID == SQLquery.results$hamap.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['hamap.def']])){
      if(colnames(SQLquery.results$hamap.def)[j] == "KissID"){ next }
      SQLquery.results$hamap.def[i,j] <- paste0(as.data.frame(SQLquery.results$hamap)[hits,j][which(as.data.frame(SQLquery.results$hamap)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$hamap)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing Hamap DB subset from goutubro SQLite3 DB \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.hamap")
  rm(gotubro.parser.subset.hamap)
  rm(hamap)
  SQLquery.results$hamap <- NULL
  
  ############################################################################################################################
  
  ### InterPro Domains
  
  cat(paste0(Sys.time(), " || Tunning InterPro DB subset ... \n"))
  
  ## gotubro.parser.subset.interpro table
  
  gotubro.parser.subset.interpro <- as.data.frame(gotubro.parser[,c(1,12)])
  colnames(gotubro.parser.subset.interpro) <- c("KissID", "ID")
  complex.rows <- grep(" | ", gotubro.parser.subset.interpro[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.interpro$ID[grep(" | ", gotubro.parser.subset.interpro$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.interpro <- apply(gotubro.parser.subset.interpro, 2, as.character)
  gotubro.parser.subset.interpro[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.interpro <- rbind(gotubro.parser.subset.interpro, data.frame(KissID = gotubro.parser.subset.interpro[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.interpro <- rbind(gotubro.parser.subset.interpro, data.frame(KissID = gotubro.parser.subset.interpro[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }  
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.interpro)
  gotubro.parser.subset.interpro <- tbl(gotubro.sqlite, "gotubro.parser.subset.interpro")
  
  ## interpro.subset table
  
  interpro <- tbl(gotubro.sqlite, "interpro")
  
  ## Querying | Joining
  
  SQLquery.results[['interpro']] <- gotubro.parser.subset.interpro %>% left_join(interpro)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['interpro']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting InterPro results query ... \n"))
  
  SQLquery.results[['interpro']] %>% collect()
  SQLquery.results[['interpro']] <- as_tibble(SQLquery.results[['interpro']])
  SQLquery.results[['interpro.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$interpro), nrow = nrow(final.df)))
  colnames(SQLquery.results$interpro.def) <- colnames(SQLquery.results$interpro)
  SQLquery.results[['interpro.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['interpro.def']])){
    hits <- which(SQLquery.results$interpro$KissID == SQLquery.results$interpro.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['interpro.def']])){
      if(colnames(SQLquery.results$interpro.def)[j] == "KissID"){ next }
      SQLquery.results$interpro.def[i,j] <- paste0(as.data.frame(SQLquery.results$interpro)[hits,j][which(as.data.frame(SQLquery.results$interpro)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$interpro)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing InterPro DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.interpro")
  rm(gotubro.parser.subset.interpro)
  rm(interpro)
  SQLquery.results$interpro <- NULL
  
  ############################################################################################################################
  
  ### PANTHER Domains
  
  cat(paste0(Sys.time(), " || Tunning PANTHER DB subset ... \n"))
  
  ## gotubro.parser.subset.panther table
  
  gotubro.parser.subset.panther <- as.data.frame(gotubro.parser[,c(1,7)])
  colnames(gotubro.parser.subset.panther) <- c("KissID", "ID")
  complex.rows <- grep(" | ", gotubro.parser.subset.panther[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.panther$ID[grep(" | ", gotubro.parser.subset.panther$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.panther <- apply(gotubro.parser.subset.panther, 2, as.character)
  gotubro.parser.subset.panther[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.panther <- rbind(gotubro.parser.subset.panther, data.frame(KissID = gotubro.parser.subset.panther[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.panther <- rbind(gotubro.parser.subset.panther, data.frame(KissID = gotubro.parser.subset.panther[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }  
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.panther)
  gotubro.parser.subset.panther <- tbl(gotubro.sqlite, "gotubro.parser.subset.panther")
  
  ## panther.subset table
  
  panther <- tbl(gotubro.sqlite, "panther")
  
  ## Querying | Joining
  
  SQLquery.results[['panther']] <- gotubro.parser.subset.panther %>% left_join(panther)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['panther']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting PANTHER results query ... \n"))
  
  SQLquery.results[['panther']] %>% collect()
  SQLquery.results[['panther']] <- as_tibble(SQLquery.results[['panther']])
  SQLquery.results[['panther']]$GO <- paste0(SQLquery.results$panther$GO_MF, SQLquery.results$panther$GO_BP, SQLquery.results$panther$GO_CC, collapse = " | ")
  SQLquery.results[['panther.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$panther), nrow = nrow(final.df)))
  colnames(SQLquery.results$panther.def) <- colnames(SQLquery.results$panther)
  SQLquery.results[['panther.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['panther.def']])){
    hits <- which(SQLquery.results$panther$KissID == SQLquery.results$panther.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['panther.def']])){
      if(colnames(SQLquery.results$panther.def)[j] == "KissID"){ next }
      SQLquery.results$panther.def[i,j] <- paste0(as.data.frame(SQLquery.results$panther)[hits,j][which(as.data.frame(SQLquery.results$panther)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$panther)[hits,j]))], collapse = " | ")
    }
  }
  SQLquery.results$panther.def$GO <- paste(SQLquery.results$panther.def$GO_MF, SQLquery.results$panther.def$GO_BP, SQLquery.results$panther.def$GO_CC)
  SQLquery.results$panther.def$GO <-  gsub(x = SQLquery.results$panther.def$GO, pattern = " ", replacement = "", fixed = T) 
  SQLquery.results$panther.def$GO <- gsub(x = SQLquery.results$panther.def$GO, pattern =  ";", replacement = " | ", fixed = T)
  SQLquery.results$panther.def$GO <- gsub(x = SQLquery.results$panther.def$GO, pattern =  "|G", replacement = "G", fixed = T)
  SQLquery.results$panther.def <- SQLquery.results$panther.def[,-c(4,5,6)]
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing PANTHER DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.panther")
  rm(gotubro.parser.subset.panther)
  rm(panther)
  SQLquery.results$panther <-NULL
  
  ############################################################################################################################
  
  ### PIRSF Domains
  
  cat(paste0(Sys.time(), " || Tunning PIRSF DB subset ... \n"))
  
  ## gotubro.parser.subset.pirsf table
  
  gotubro.parser.subset.pirsf <- as.data.frame(gotubro.parser[,c(1,11)])
  colnames(gotubro.parser.subset.pirsf) <- c("KissID", "ID")
  complex.rows <- grep(" | ", gotubro.parser.subset.pirsf[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.pirsf$ID[grep(" | ", gotubro.parser.subset.pirsf$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.pirsf <- apply(gotubro.parser.subset.pirsf, 2, as.character)
  gotubro.parser.subset.pirsf[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.pirsf <- rbind(gotubro.parser.subset.pirsf, data.frame(KissID = gotubro.parser.subset.pirsf[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.pirsf <- rbind(gotubro.parser.subset.pirsf, data.frame(KissID = gotubro.parser.subset.pirsf[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.pirsf)
  gotubro.parser.subset.pirsf <- tbl(gotubro.sqlite, "gotubro.parser.subset.pirsf")
  
  ## pirsf.subset table
  
  pirsf <- tbl(gotubro.sqlite, "pirsf")
  
  ## Querying | Joining
  
  SQLquery.results[['pirsf']] <- gotubro.parser.subset.pirsf %>% left_join(pirsf)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['pirsf']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting PIRSF results query ... \n"))
  
  SQLquery.results[['pirsf']] %>% collect()
  SQLquery.results[['pirsf']] <- as_tibble(SQLquery.results[['pirsf']])
  SQLquery.results[['pirsf.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$pirsf), nrow = nrow(final.df)))
  colnames(SQLquery.results$pirsf.def) <- colnames(SQLquery.results$pirsf)
  SQLquery.results[['pirsf.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['pirsf.def']])){
    hits <- which(SQLquery.results$pirsf$KissID == SQLquery.results$pirsf.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['pirsf.def']])){
      if(colnames(SQLquery.results$pirsf.def)[j] == "KissID"){ next }
      SQLquery.results$pirsf.def[i,j] <- paste0(as.data.frame(SQLquery.results$pirsf)[hits,j][which(as.data.frame(SQLquery.results$pirsf)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$pirsf)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing PIRSF DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.pirsf")
  rm(gotubro.parser.subset.pirsf)
  rm(pirsf)
  SQLquery.results$pirsf <- NULL
  
  ############################################################################################################################
  
  ### PRINTS Domains
  
  cat(paste0(Sys.time(), " || Tunning PRINTS DB subset ... \n"))
  
  ## gotubro.parser.subset.prints table
  
  gotubro.parser.subset.prints <- as.data.frame(gotubro.parser[,c(1,9)])
  colnames(gotubro.parser.subset.prints) <- c("KissID", "ID")
  complex.rows <- grep(" | ", gotubro.parser.subset.prints[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.prints$ID[grep(" | ", gotubro.parser.subset.prints$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.prints <- apply(gotubro.parser.subset.prints, 2, as.character)
  gotubro.parser.subset.prints[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.prints <- rbind(gotubro.parser.subset.prints, data.frame(KissID = gotubro.parser.subset.prints[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.prints <- rbind(gotubro.parser.subset.prints, data.frame(KissID = gotubro.parser.subset.prints[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.prints)
  gotubro.parser.subset.prints <- tbl(gotubro.sqlite, "gotubro.parser.subset.prints")
  
  ## prints.subset table
  
  prints <- tbl(gotubro.sqlite, "prints")
  
  ## Querying | Joining
  
  SQLquery.results[['prints']] <- gotubro.parser.subset.prints %>% left_join(prints)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['prints']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting PRINTS results query ... \n"))
  
  SQLquery.results[['prints']] %>% collect()
  SQLquery.results[['prints']] <- as_tibble(SQLquery.results[['prints']])
  SQLquery.results[['prints.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$prints), nrow = nrow(final.df)))
  colnames(SQLquery.results$prints.def) <- colnames(SQLquery.results$prints)
  SQLquery.results[['prints.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['prints.def']])){
    hits <- which(SQLquery.results$prints$KissID == SQLquery.results$prints.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['prints.def']])){
      if(colnames(SQLquery.results$prints.def)[j] == "KissID"){ next }
      SQLquery.results$prints.def[i,j] <- paste0(as.data.frame(SQLquery.results$prints)[hits,j][which(as.data.frame(SQLquery.results$prints)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$prints)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing PRINTS DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.prints")
  rm(gotubro.parser.subset.prints)
  rm(prints)
  SQLquery.results$prints <- NULL
  
  ############################################################################################################################
  
  ### PROSITE Domains
  
  cat(paste0(Sys.time(), " || Tunning ProSite DB subset ... \n"))
  
  ## gotubro.parser.subset.prosite table
  
  gotubro.parser.subset.prosite <- as.data.frame(gotubro.parser[,c(1,10)])
  colnames(gotubro.parser.subset.prosite) <- c("KissID", "ID")
  complex.rows <- grep(" | ", gotubro.parser.subset.prosite[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.prosite$ID[grep(" | ", gotubro.parser.subset.prosite$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.prosite <- apply(gotubro.parser.subset.prosite, 2, as.character)
  gotubro.parser.subset.prosite[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.prosite <- rbind(gotubro.parser.subset.prosite, data.frame(KissID = gotubro.parser.subset.prosite[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.prosite <- rbind(gotubro.parser.subset.prosite, data.frame(KissID = gotubro.parser.subset.prosite[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.prosite)
  gotubro.parser.subset.prosite <- tbl(gotubro.sqlite, "gotubro.parser.subset.prosite")
  
  ## prosite.subset table
  
  prosite <- tbl(gotubro.sqlite, "prosite")
  
  ## Querying | Joining
  
  SQLquery.results[['prosite']] <- gotubro.parser.subset.prosite %>% left_join(prosite)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['prosite']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting ProSite results query ... \n"))
  
  SQLquery.results[['prosite']] %>% collect()
  SQLquery.results[['prosite']] <- as_tibble(SQLquery.results[['prosite']])
  SQLquery.results[['prosite.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$prosite), nrow = nrow(final.df)))
  colnames(SQLquery.results$prosite.def) <- colnames(SQLquery.results$prosite)
  SQLquery.results[['prosite.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['prosite.def']])){
    hits <- which(SQLquery.results$prosite$KissID == SQLquery.results$prosite.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['prosite.def']])){
      if(colnames(SQLquery.results$prosite.def)[j] == "KissID"){ next }
      SQLquery.results$prosite.def[i,j] <- paste0(as.data.frame(SQLquery.results$prosite)[hits,j][which(as.data.frame(SQLquery.results$prosite)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$prosite)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing ProSite DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.prosite")
  rm(gotubro.parser.subset.prosite)
  rm(prosite)
  SQLquery.results$prosite <- NULL
  
  ############################################################################################################################
  
  ### Pfam Domains
  
  cat(paste0(Sys.time(), " || Tunning Pfam DB subset ... \n"))
  
  ## gotubro.parser.subset.pfam table
  
  gotubro.parser.subset.pfam <- as.data.frame(gotubro.parser[,c(1,4)])
  colnames(gotubro.parser.subset.pfam) <- c("KissID", "ID")
  gotubro.parser.subset.pfam$ID <- gsub(x = gotubro.parser.subset.pfam$ID, pattern = "Dbxref=Pfam:", replacement = "")
  complex.rows <- grep(" | ", gotubro.parser.subset.pfam[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.pfam$ID[grep(" | ", gotubro.parser.subset.pfam$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.pfam <- apply(gotubro.parser.subset.pfam, 2, as.character)
  gotubro.parser.subset.pfam[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.pfam <- rbind(gotubro.parser.subset.pfam, data.frame(KissID = gotubro.parser.subset.pfam[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.pfam <- rbind(gotubro.parser.subset.pfam, data.frame(KissID = gotubro.parser.subset.pfam[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }
  }
  gotubro.parser.subset.pfam$ID <- do.call(rbind, strsplit(as.character(gotubro.parser.subset.pfam$ID), ".", fixed = T))[,1]
  copy_to(gotubro.sqlite, gotubro.parser.subset.pfam)
  gotubro.parser.subset.pfam <- tbl(gotubro.sqlite, "gotubro.parser.subset.pfam")
  
  ## pfam.subset table
  
  pfam <- tbl(gotubro.sqlite, "pfam")
  
  ## Querying | Joining
  
  SQLquery.results[['pfam']] <- gotubro.parser.subset.pfam %>% left_join(pfam)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['pfam']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting Pfam results query ... \n"))
  
  SQLquery.results[['pfam']] %>% collect()
  SQLquery.results[['pfam']] <- as_tibble(SQLquery.results[['pfam']])
  SQLquery.results[['pfam.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$pfam), nrow = nrow(final.df)))
  colnames(SQLquery.results$pfam.def) <- colnames(SQLquery.results$pfam)
  SQLquery.results[['pfam.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['pfam.def']])){
    hits <- which(SQLquery.results$pfam$KissID == SQLquery.results$pfam.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['pfam.def']])){
      if(colnames(SQLquery.results$pfam.def)[j] == "KissID"){ next }
      SQLquery.results$pfam.def[i,j] <- paste0(as.data.frame(SQLquery.results$pfam)[hits,j][which(as.data.frame(SQLquery.results$pfam)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$pfam)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing Pfam DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.pfam")
  rm(gotubro.parser.subset.pfam)
  rm(pfam)
  SQLquery.results$pfam <- NULL
  
  ############################################################################################################################
  
  ### Rfam Domains
  
  cat(paste0(Sys.time(), " || Tunning Rfam DB subset ... \n"))
  
  ## gotubro.parser.subset.rfam table
  
  gotubro.parser.subset.rfam <- as.data.frame(gotubro.parser[,c(1,5)])
  colnames(gotubro.parser.subset.rfam) <- c("KissID", "ID")
  gotubro.parser.subset.rfam$ID <- gsub(x = gotubro.parser.subset.rfam$ID, pattern = "Dbxref=Rfam:", replacement = "")
  complex.rows <- grep(" | ", gotubro.parser.subset.rfam[,2], fixed = T)
  complex.names <- str_split_fixed(string = gotubro.parser.subset.rfam$ID[grep(" | ", gotubro.parser.subset.rfam$ID, fixed = T)], " | ", n = Inf)
  complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
  gotubro.parser.subset.rfam <- apply(gotubro.parser.subset.rfam, 2, as.character)
  gotubro.parser.subset.rfam[complex.rows,2] <- complex.names[,1]
  gotubro.parser.subset.rfam <- rbind(gotubro.parser.subset.rfam, data.frame(KissID = gotubro.parser.subset.rfam[complex.rows,1], ID = complex.names[,2]))
  if(ncol(complex.names) > 2){
    for(Zz in 3:ncol(complex.names)){
      gotubro.parser.subset.rfam <- rbind(gotubro.parser.subset.rfam, data.frame(KissID = gotubro.parser.subset.rfam[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
    }
  }
  copy_to(gotubro.sqlite, gotubro.parser.subset.rfam)
  gotubro.parser.subset.rfam <- tbl(gotubro.sqlite, "gotubro.parser.subset.rfam")
  
  ## rfam.subset table
  
  rfam <- tbl(gotubro.sqlite, "rfam")
  
  ## Querying | Joining
  
  SQLquery.results[['rfam']] <- gotubro.parser.subset.rfam %>% left_join(rfam)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['rfam']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting Rfam results query ... \n"))
  
  SQLquery.results[['rfam']] %>% collect()
  SQLquery.results[['rfam']] <- as_tibble(SQLquery.results[['rfam']])
  SQLquery.results[['rfam.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$rfam), nrow = nrow(final.df)))
  colnames(SQLquery.results$rfam.def) <- colnames(SQLquery.results$rfam)
  SQLquery.results[['rfam.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['rfam.def']])){
    hits <- which(SQLquery.results$rfam$KissID == SQLquery.results$rfam.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['rfam.def']])){
      if(colnames(SQLquery.results$rfam.def)[j] == "KissID"){ next }
      SQLquery.results$rfam.def[i,j] <- paste0(as.data.frame(SQLquery.results$rfam)[hits,j][which(as.data.frame(SQLquery.results$rfam)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$rfam)[hits,j]))], collapse = " | ")
    }
  }
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing Rfam DB subset from goutubro SQLite3 DB \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.rfam")
  rm(gotubro.parser.subset.rfam)
  rm(rfam)
  SQLquery.results$rfam <- NULL
  
  ############################################################################################################################
  
  ### Sprot DB
  
  cat(paste0(Sys.time(), " || Tunning Sprot DB subset ... \n"))
  
  ## gotubro.parser.subset.sprot table
  
  gotubro.parser.subset.sprot <- as.data.frame(gotubro.parser[,c(1,3)])
  gotubro.parser.subset.sprot[,2] <- gsub("Name=","",gotubro.parser.subset.sprot[,2], fixed = T)
  colnames(gotubro.parser.subset.sprot) <- c("KissID", "UniprotKB-AC")
  gotubro.parser.subset.sprot[,2] <- do.call(rbind, strsplit(x = as.character(gotubro.parser.subset.sprot$`UniprotKB-AC`), split = "|", fixed = T))[,2]
  copy_to(gotubro.sqlite, gotubro.parser.subset.sprot)
  gotubro.parser.subset.sprot <- tbl(gotubro.sqlite, "gotubro.parser.subset.sprot")
  
  ## sprot.subset table
  
  sprot <- tbl(gotubro.sqlite, "uniprot")
  sprot <- sprot %>%
    select(`UniprotKB-AC`, GO)
  
  ## Querying | Joining
  
  SQLquery.results[['sprot']] <- gotubro.parser.subset.sprot %>% left_join(sprot)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['sprot']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting Sprot results query ... \n"))
  crawl(glue::glue("This process could take a bit (3-4 h). While you are waiting here´s a chicken.
                   {things$chicken %>% stringr::str_replace('%s', 'Fvck Mario Casas y que viva el pirri')}"),
        colors = RColorBrewer::brewer.pal(9, "YlOrRd")[3:7],
        direction = "vertical",
        pause = 0.04) 
  
  SQLquery.results[['sprot']] %>% collect()
  SQLquery.results[['sprot']] <- as_tibble(SQLquery.results[['sprot']])
  SQLquery.results[['sprot.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$sprot), nrow = nrow(final.df)))
  colnames(SQLquery.results$sprot.def) <- colnames(SQLquery.results$sprot)
  SQLquery.results[['sprot.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['sprot.def']])){
    hits <- which(SQLquery.results$sprot$KissID == SQLquery.results$sprot.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['sprot.def']])){
      if(colnames(SQLquery.results$sprot.def)[j] == "KissID"){ next }
      SQLquery.results$sprot.def[i,j] <- paste0(as.data.frame(SQLquery.results$sprot)[hits,j][which(as.data.frame(SQLquery.results$sprot)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$sprot)[hits,j]))], collapse = " | ")
    }
  }
  SQLquery.results$sprot.def$GO <- gsub("; ", " | ", SQLquery.results$sprot.def$GO, fixed = T)
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing Sprot DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.sprot")
  rm(gotubro.parser.subset.sprot)
  rm(sprot)
  SQLquery.results$sprot <- NULL
  
  ############################################################################################################################
  
  ### Uniref90 Cluster
  
  cat(paste0(Sys.time(), " || Tunning Uniref90 DB subset ... \n"))
  
  ## gotubro.parser.subset.uniref table
  
  gotubro.parser.subset.uniref <- as.data.frame(gotubro.parser[,c(1,2)])
  gotubro.parser.subset.uniref[,2] <- gsub("Name=","",gotubro.parser[,2], fixed = T)
  colnames(gotubro.parser.subset.uniref) <- c("KissID", "Uniref90")
  copy_to(gotubro.sqlite, gotubro.parser.subset.uniref)
  gotubro.parser.subset.uniref <- tbl(gotubro.sqlite, "gotubro.parser.subset.uniref")
  
  ## uniref.subset table
  
  uniref <- tbl(gotubro.sqlite, "uniprot")
  uniref <- uniref %>%
    select(`UniprotKB-AC`, GO, Uniref90)
  
  ## Querying | Joining
  
  SQLquery.results[['uniref']] <- gotubro.parser.subset.uniref %>% left_join(uniref)
  
  cat(paste0(Sys.time(), " || SQL code executed; \n"))
  SQLquery.results[['uniref']] %>% show_query()
  cat(paste0(Sys.time(), " || Collecting Uniref90 results query ... \n"))
  crawl(glue::glue("This process could take a bit (9-10 h). While you are waiting here´s a high af toad.
                   {things$hypnotoad %>% stringr::str_replace('%s', 'Everybody go to the discoteke')}"),
        colors = RColorBrewer::brewer.pal(9, "YlOrRd")[3:7],
        direction = "vertical",
        pause = 0.04) 
  
  SQLquery.results[['uniref']] %>% collect()
  SQLquery.results[['uniref']] <- as_tibble(SQLquery.results[['uniref']])
  SQLquery.results[['uniref.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$uniref), nrow = nrow(final.df)))
  colnames(SQLquery.results$uniref.def) <- colnames(SQLquery.results$uniref)
  SQLquery.results[['uniref.def']][,1] <- final.df[,1]
  for(i in 1:nrow(SQLquery.results[['uniref.def']])){
    hits <- which(SQLquery.results$uniref$KissID == SQLquery.results$uniref.def$KissID[i])
    for(j in 1:ncol(SQLquery.results[['uniref.def']])){
      if(colnames(SQLquery.results$uniref.def)[j] == "KissID"){ next }
      SQLquery.results$uniref.def[i,j] <- paste0(as.data.frame(SQLquery.results$uniref)[hits,j][which(as.data.frame(SQLquery.results$uniref)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$uniref)[hits,j]))], collapse = " | ")
    }
  }
  SQLquery.results$uniref.def$GO <- gsub("; ", " | ", SQLquery.results$uniref.def$GO, fixed = T)
  gc()
  
  ## Remove remote-table from gotubro SQLite3 DB
  
  cat(paste0(Sys.time(), " || Removing Uniref90 DB subset from goutubro SQLite3 DB  \n"))
  DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.uniref")
  rm(gotubro.parser.subset.uniref)
  rm(uniref)
  SQLquery.results$uniref <- NULL
  
  ############################################################################################################################
  
  ### QuickGO complete
  
  if(quickGO){
    
    cat(paste0(Sys.time(), " || QuickGO DB subset ... \n"))
    
    ## gotubro.parser.subset.quick table
    
    gotubro.parser.subset.quick <- SQLquery.results$sprot.def[,c(1,2)]
    if(length(which((gotubro.parser.subset.quick[,1] == SQLquery.results$uniref.def[,1]) == TRUE)) == nrow(gotubro.parser.subset.quick)){
      gotubro.parser.subset.quick[which(gotubro.parser.subset.quick[,2] == "" & SQLquery.results$uniref.def$`UniprotKB-AC` != ""),2] <- SQLquery.results$uniref.def$`UniprotKB-AC`[which(gotubro.parser.subset.quick[,2] == "" & SQLquery.results$uniref.def$`UniprotKB-AC` != "")]  
    }else if(length(which((gotubro.parser.subset.quick[,1] == SQLquery.results$uniref.def[,1]) == TRUE)) != nrow(gotubro.parser.subset.quick)){
      stop(say(what = "\n Something went wrong! :( Ask Victor: fernandezvictor@uniovi.es \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))
    }
    colnames(gotubro.parser.subset.quick) <- c("KissID", "ID")
    complex.rows <- grep(" | ", gotubro.parser.subset.quick[,2], fixed = T)
    complex.names <- str_split_fixed(string = gotubro.parser.subset.quick$ID[grep(" | ", gotubro.parser.subset.quick$ID, fixed = T)], " | ", n = Inf)
    complex.names <- complex.names[,which(1:ncol(complex.names) %% 2 == 1)]
    gotubro.parser.subset.quick <- apply(gotubro.parser.subset.quick, 2, as.character)
    gotubro.parser.subset.quick[complex.rows,2] <- complex.names[,1]
    gotubro.parser.subset.quick <- rbind(gotubro.parser.subset.quick, data.frame(KissID = gotubro.parser.subset.quick[complex.rows,1], ID = complex.names[,2]))
    if(ncol(complex.names) > 2){
      for(Zz in 3:ncol(complex.names)){
        gotubro.parser.subset.quick <- rbind(gotubro.parser.subset.quick, data.frame(KissID = gotubro.parser.subset.quick[complex.rows[-which(complex.names[,Zz] == "")],1], ID = complex.names[-which(complex.names[,Zz] == ""), Zz]))
      }  
    }
    copy_to(gotubro.sqlite, gotubro.parser.subset.quick)
    gotubro.parser.subset.quick <- tbl(gotubro.sqlite, "gotubro.parser.subset.quick")
    
    ## quick.subset table
    
    quick <- tbl(gotubro.sqlite, "quick")
    quick <- quick %>%
      select(ID, Symbol, GO, DBXrefs)
    
    ## Querying | Joining
    
    SQLquery.results[['quick']] <- gotubro.parser.subset.quick %>% left_join(quick)
    
    cat(paste0(Sys.time(), " || SQL code executed; \n"))
    SQLquery.results[['quick']] %>% show_query()
    cat(paste0(Sys.time(), " || Collecting QuickGO results query ... \n"))
    crawl(glue::glue("This process could take a bit (3-4 h). While you are waiting here´s a useful advice.
                   {things$yoda %>% stringr::str_replace('%s', 'Wise choice! A good analyst if you are, quickGO you always chose')}"),
          colors = RColorBrewer::brewer.pal(9, "YlOrRd")[3:7],
          direction = "vertical",
          pause = 0.04) 
    
    SQLquery.results[['quick']] %>% collect()
    SQLquery.results[['quick']] <- as_tibble(SQLquery.results[['quick']])
    SQLquery.results[['quick.def']] <- as.data.frame(matrix(NA, ncol = ncol(SQLquery.results$quick), nrow = nrow(final.df)))
    colnames(SQLquery.results$quick.def) <- colnames(SQLquery.results$quick)
    SQLquery.results[['quick.def']][,1] <- final.df[,1]
    for(i in 1:nrow(SQLquery.results[['quick.def']])){
      hits <- which(SQLquery.results$quick$KissID == SQLquery.results$quick.def$KissID[i])
      for(j in 1:ncol(SQLquery.results[['quick.def']])){
        if(colnames(SQLquery.results$quick.def)[j] == "KissID"){ next }
        SQLquery.results$quick.def[i,j] <- paste0(as.data.frame(SQLquery.results$quick)[hits,j][which(as.data.frame(SQLquery.results$quick)[hits,j] != "NA" & !duplicated(as.data.frame(SQLquery.results$quick)[hits,j]))], collapse = " | ")
      }
    }
    gc()
    
    ## Remove remote-table from gotubro SQLite3 DB
    
    cat(paste0(Sys.time(), " || Removing QuickGO DB subset from goutubro SQLite3 DB  \n"))
    DBI::dbRemoveTable(gotubro.sqlite, "gotubro.parser.subset.quick")
    rm(gotubro.parser.subset.quick)
    rm(quick)
    SQLquery.results$quick <- NULL
  }
  
  ### Output formatting # I NEED REVIEW WITH FINAL RESULTS
  
  cat(paste0(Sys.time(), " || Processing final outputs ...  \n"))
  final.list <- list()
  
  ## Quality control
  
  for(Zzz in 1:length(SQLquery.results)){
    test <- length(which((final.df[,1] == SQLquery.results[[Zzz]][,1]) == TRUE))
    if(test != nrow(final.df)){
      stop(say(what = "\n Something went wrong! :( Ask Victor: fernandezvictor@uniovi.es \n", by = "poop", type = "warning", what_color = "gold", by_color = "gold"))
    }
  }
  
  ## Formatting
  
  if("reduce" %in% mode){ 
   for(Zzzz in 1:length(SQLquery.results)){
    colnames(SQLquery.results[[Zzzz]]) <- paste0(colnames(SQLquery.results[[Zzzz]]), "_", names(SQLquery.results)[Zzzz]) 
   }
    final.list[['reduce']] <- do.call(cbind, SQLquery.results)
    final.list[['reduce']] <- final.list[['reduce']][,c(1,grep("GO_",colnames(final.list[['reduce']]), fixed = T))]
    final.list$reduce$GO_gotubro <- rep(NA, nrow(final.list$reduce))
    for(i in 1:nrow(final.list$reduce)){
      complete.row <- paste0(as.character(unlist(final.list$reduce[i,which(final.list$reduce[i,] != "")[-1]], use.names = F)), collapse = " | ")
      no.dups <- paste0(unlist(strsplit(complete.row, split = " | ", fixed = T))[!duplicated(unlist(strsplit(complete.row, split = " | ", fixed = T)))], collapse = " | ")
      final.list$reduce$GO_gotubro[i] <- no.dups
    }
    final.list[['reduce']] <- final.list[['reduce']][,c(1, ncol(final.list[['reduce']]))]
  }
      
  if("complete" %in% mode){
    for(Zzzz in 1:length(SQLquery.results)){
      colnames(SQLquery.results[[Zzzz]]) <- paste0(colnames(SQLquery.results[[Zzzz]]), "_", names(SQLquery.results)[Zzzz]) 
    }
    final.list[['complete']] <- do.call(cbind, SQLquery.results)
    final.list[['complete']] <- final.list[['complete']][,-grep("ID", colnames(final.list$complete))[-1]]
    final.list[['complete']] <- final.list[['complete']][,-grep("UniprotKB", colnames(final.list$complete))]
    final.list[['complete']] <- final.list[['complete']][,-grep("Uniref90", colnames(final.list$complete))]
  }
  
  if("normal" %in% mode){
    for(Zzzz in 1:length(SQLquery.results)){
      colnames(SQLquery.results[[Zzzz]]) <- paste0(colnames(SQLquery.results[[Zzzz]]), "_", names(SQLquery.results)[Zzzz]) 
    }
    final.list[['normal']] <- do.call(cbind, SQLquery.results)
    final.list[['normal']] <- final.list[['normal']][,c(1,grep("GO_",colnames(final.list[['normal']]), fixed = T))]
  }
  
  final.list <- lapply(final.list, function(x){
    x <- apply(x, 2, function(y){
       y[which(y == "")] <- NA
       y
    })
    x
  })
  
  if(return.merged){
    final.list <- lapply(final.list, function(x){
      sql.df <- cbind(gotubro.parser.id, x[,-1])
      colnames(sql.df) <- c(colnames(gotubro.parser.id), colnames(x)[-1])
      sql.df
    })
  }
  
  ### Annotated Percentage || Same as reduce with final list results
  
  annotated <- length(which((rowSums(is.na(final.list[[1]][,grep("GO_", colnames(final.list[[1]])), drop = FALSE])) != ncol(final.list[[1]][,grep("GO_", colnames(final.list[[1]])), drop = FALSE])) == TRUE))
  percentage.annotation <- round(annotated/nrow(gotubro.parser) * 100, 2)
  cat(paste0("\n ", percentage.annotation, " % sequences have been annotated (", annotated, ") \n"))
  
  ### Say Bye
  
  say(what = paste0(Sys.time(), " || Thanks for using it! <3" ), by = "rabbit", what_color = "salmon", by_color = "salmon")

  ### Disconnect from gotubro SQLite3 DB
  
  DBI::dbDisconnect(gotubro.sqlite)
  
  ### Return final list
  
  return(final.list)
}