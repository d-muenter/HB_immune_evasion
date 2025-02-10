#.libPaths('~/R/library/')
library(Seurat)
library(tidyverse)
library(future)
library(writexl)

setwd("")
source('seurat_functions_public.R')

sample_list <- c("286", "287", "288", "289", "386", "387", "388", "389", "391", "392")

Genes_nmf_w_basis <- list()

for (i in sample_list){
  
  load(paste0("Data/H-", i, "/res.RData"))  # load the "final" result based on the Barkley analysis
  rank <- summary(res)["rank"]  # find out the rank of that result
  
  res_list_file <- list.files(paste0("Data/H-", i, "/"))[grepl("res.list", list.files(paste0("Data/H-", i, "/")))]
  load(paste0("Data/H-", i, "/", res_list_file))  # load the list with all calculated results of that sample
  
  res.list <- res.list[1:(rank-1)]  # subset the list on all results with the rank of the final results and below
  sample_allNMF <- data.frame(matrix(nrow = 2000))  # create dataframe that will contain all nmf scores from all ranks and modules and will be one element of the final genes_nmf_w_basis list

  # extract nmf scores
  for(j in 1:length(res.list)){
    scores <- basis(res.list[[j]]) # extract scores for each rank
    cols_pre <- ncol(sample_allNMF) # get number of columns of dataframe in order to later be able to change the correct column names
    sample_allNMF <- cbind(sample_allNMF, scores) # add scores to the score dataframe for the respective sample
    cols <- ncol(sample_allNMF)
    nmodules <- cols - cols_pre
    colnames(sample_allNMF)[(cols_pre+1):cols] <- paste0("H-", i, "_rank2_", rank, "_nruns1.RDS.", j+1, ".", 1:nmodules)
    
    
  }
  
  sample_allNMF <- sample_allNMF[,-1]
  
  if (i == sample_list[1]){
    Genes_nmf_w_basis <- list(sample_allNMF)
  } else {
    temp <- list(sample_allNMF)
    Genes_nmf_w_basis <- c(Genes_nmf_w_basis, temp)
  }

  names(Genes_nmf_w_basis)[length(Genes_nmf_w_basis)] <- paste0("H-", i, "_rank2_", rank, "_nruns1.RDS")
}

save(Genes_nmf_w_basis, file = "Output/01_Genes_nmf_w_basis_hb10x.RData")





