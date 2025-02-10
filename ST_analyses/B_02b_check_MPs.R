library(Seurat)
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(writexl)

setwd("")

load("Data/merged_integration_ST_postSCIntegration.RData")
DefaultAssay(dataIntegrated) <- "SCT"


method = "major_8-3_no-res"

mp_genes <- as.data.frame(read_xlsx(paste0("Output/", method, "/02a_NMF+MP_hb10x_", method, ".xlsx"), sheet = 3))

writeTo <- paste0("Output/", method, "/check_MPs/")
dir.create(writeTo)

# calculate MP scores
rawData <- GetAssayData(dataIntegrated, slot = "data")
for(i in 1:ncol(mp_genes)){
  genes <- mp_genes[, i]
  gene.set <- genes[genes %in% row.names(rawData)]
  
  mean.exp <- colMeans(rawData[gene.set, ], na.rm = TRUE)
  if (all(names(x = mean.exp) == rownames(x = dataIntegrated@meta.data))) {
    cat(paste0("\n", colnames(mp_genes)[i]))
    dataIntegrated@meta.data[colnames(mp_genes)[i]] <- mean.exp
  }
}

mp_genes_long <- data.frame(genes = as.vector(as.matrix(mp_genes)), MP = rep(colnames(mp_genes), each=50))
MP_geneExpr <- t(as.data.frame(rawData[rownames(rawData) %in% mp_genes_long$genes,])) # get gene expression for all MP genes
MP_geneExpr <- MP_geneExpr[,unique(mp_genes_long$genes[mp_genes_long$genes %in% colnames(MP_geneExpr)])] # order genes by MP gene list

MP_scores <- dataIntegrated@meta.data[,colnames(mp_genes)] # get MP scores

sum(rownames(MP_geneExpr) != rownames(MP_scores)) == 0 # make sure rownames are identical

MP_cor <- as.data.frame(cor(MP_geneExpr, MP_scores)) # calculate correlations


gene_anno <- mp_genes_long[!duplicated(mp_genes_long$genes),]
rownames(gene_anno) <- gene_anno$genes
gene_anno <- gene_anno[gene_anno$genes %in% colnames(MP_geneExpr),]

pheatmap(MP_cor, cluster_rows = F, cluster_cols = F, show_rownames = F, 
         annotation_row = gene_anno[c("MP")], 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         filename = paste0(writeTo, "MP_gene_correlation.tiff"), 
         width = 10, height = 10, dpi = 300)

# filter genes to only include genes showing a Pearson correlation coefficient >= 0.2
# with their respective MP and MPs with >= 10 filtered genes
mp_genes_fil <- data.frame(genes=NA, MP=NA)
for (i in colnames(MP_cor)){
  genes <- gene_anno$genes[gene_anno$MP == i]
  MP_cor_fil <- MP_cor[genes,]
  genes_fil <- genes[MP_cor_fil[i] >= 0.2]
  mp_genes_temp <- gene_anno[genes_fil,]
  if (nrow(mp_genes_temp) >= 10){
    mp_genes_fil <- rbind(mp_genes_fil, mp_genes_temp)
  }
}
mp_genes_fil <- mp_genes_fil[-1,]

MP_cor_fil <- MP_cor[mp_genes_fil$genes,]

pheatmap(MP_cor_fil, cluster_rows = F, cluster_cols = F, show_rownames = F, 
         annotation_row = gene_anno[c("MP")], 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         filename = paste0(writeTo, "MP_gene_correlation_filtered.tiff"), 
         width = 10, height = 10, dpi = 300)

mp_genes_fil_df <- data.frame(MP_1 = c(mp_genes_fil$genes[mp_genes_fil$MP == "MP_1"], rep(NA, (50 - length(mp_genes_fil$genes[mp_genes_fil$MP == "MP_1"])))), 
                              MP_2 = c(mp_genes_fil$genes[mp_genes_fil$MP == "MP_2"], rep(NA, (50 - length(mp_genes_fil$genes[mp_genes_fil$MP == "MP_2"])))), 
                              MP_3 = c(mp_genes_fil$genes[mp_genes_fil$MP == "MP_3"], rep(NA, (50 - length(mp_genes_fil$genes[mp_genes_fil$MP == "MP_3"])))), 
                              MP_4 = c(mp_genes_fil$genes[mp_genes_fil$MP == "MP_4"], rep(NA, (50 - length(mp_genes_fil$genes[mp_genes_fil$MP == "MP_4"])))), 
                              MP_5 = c(mp_genes_fil$genes[mp_genes_fil$MP == "MP_5"], rep(NA, (50 - length(mp_genes_fil$genes[mp_genes_fil$MP == "MP_5"])))), 
                              MP_6 = c(mp_genes_fil$genes[mp_genes_fil$MP == "MP_6"], rep(NA, (50 - length(mp_genes_fil$genes[mp_genes_fil$MP == "MP_6"])))), 
                              MP_7 = c(mp_genes_fil$genes[mp_genes_fil$MP == "MP_7"], rep(NA, (50 - length(mp_genes_fil$genes[mp_genes_fil$MP == "MP_7"])))))



write_xlsx(mp_genes_fil_df, path = paste0("Output/", method, "/02b_filtered_MPs_", method, ".xlsx"))
