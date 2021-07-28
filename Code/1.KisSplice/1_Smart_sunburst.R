########################################################################################
############################ 1. KISSPLICE ##############################################
########################################################################################

# Auxiliar function needed to pass from each kmer raw counts to each kmer presence/absence matrix

#SMARTSUNBURST
sunburstSMART <- function(..., treatmentcolumn, initialcolumn, initialrow, nsamples){
  smart_list <- list(...) 
  smart_mat <- do.call(bind_rows, smart_list)
  smart_mat[is.na(smart_mat)] <- 0
  nmat<-length(smart_list)
  list_mat <- list()
  for (i in 1:nmat) {
    r <- i-1
    k <- colSums(smart_mat[(r*nsamples +initialrow): (i*nsamples),initialcolumn:ncol(smart_mat)])
    list_mat[[i]] <- k
    }
  FINAL <- do.call(bind_rows, list_mat)
  FINAL [FINAL >0] <- 1
  FINAL <- cbind("id" = c(paste0("matrix", 1:nmat)), FINAL)
  return(FINAL)
}
 