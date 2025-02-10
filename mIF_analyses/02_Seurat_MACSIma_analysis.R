library(Seurat)
library(tidyverse)


for (i in 1:7){ # for each sample
  
  sample <- readRDS(paste0("sample_", i, "/Sample-", i, "_Seurat_raw.rds"))
  DefaultAssay(sample) <- "cell_exp"
  
  writeTo <- paste0("sample_", i, "/SeuratOutput/")
  dir.create(paste0(writeTo, "TumorDefinition/"))
  
  
  # add tumor classification
  coord <- as.data.frame(sample@reductions$cell_coord@cell.embeddings)
  
  tumor_r1 <- read.csv(paste0("sample_", i, "/sample_", i, "_Crop_1/sample_", i, "_Crop_1_expr_tumor_gate.csv"))
  tumor_r2 <- read.csv(paste0("sample_", i, "/sample_", i, "_Crop_2/sample_", i, "_Crop_2_expr_tumor_gate.csv"))
  tumor_r3 <- read.csv(paste0("sample_", i, "/sample_", i, "_Crop_3/sample_", i, "_Crop_3_expr_tumor_gate.csv"))
  tumor_r1 <- paste(tumor_r1[,73], tumor_r1[,75], sep="_")
  tumor_r2 <- paste(tumor_r2[,73], tumor_r2[,75], sep="_")
  tumor_r3 <- paste(tumor_r3[,73], tumor_r3[,75], sep="_")
  
  
  coord$pos <- paste(coord[,1], coord[,2], sep="_")
  coord$classification <- "non-tumor"
  coord$classification[coord$pos %in% c(tumor_r1, tumor_r2, tumor_r3)] <- "tumor"
  
  
  sample_meta <- sample@meta.data
  sample_meta$cellID <- rownames(sample_meta)
  sample_meta$order <- 1:nrow(sample_meta)
  coord$cellID <- rownames(coord)
  coord_merge <- coord[,4:5]
  sum(sample_meta$cellID != coord_merge$cellID)
  sample_meta_merge <- merge(sample_meta, coord_merge, by="cellID")
  sample_meta_merge <- sample_meta_merge[order(sample_meta_merge$order),]
  rownames(sample_meta_merge) <- sample_meta_merge$cellID
  sample_meta_merge$order <- NULL
  sample@meta.data <- sample_meta_merge
  
  
  d1 <- DimPlot(sample, reduction = "cell_coord", pt.size=0.5, group.by="classification", raster=F)
  
  tiff(paste0(writeTo, "TumorDefinition/tumor_non-tumor_DimPlot.tiff"), width=700, height=800)
  print(d1)
  dev.off()
  

  saveRDS(sample, file = paste0("sample_", i, "/Sample-", i, "_Seurat_raw.rds"))
  
  
  
  
  # coordinates for density plots
  coord <- as.data.frame(sample@reductions$cell_coord@cell.embeddings)
  x_min <- min(coord$coord_1)
  x_max <- max(coord$coord_1)
  y_min <- min(coord$coord_2)
  y_max <- max(coord$coord_2)
  
  
  # quantification based on lymphoid cell signature
  dir.create(paste0(writeTo, "immune_quantification/lymphoid_signature/"), r=T)
  
  lymphoid_sign <- c("CD3", "CD4", "CD8a", "CD20", "CD45", "CD56")
  
  # select proteins of interest
  lymphoid <- rownames(sample@assays$SCT[rownames(sample@assays$SCT) %in% lymphoid_sign])
  # get mean expression of proteins of interest per cell
  mean.exp <- colMeans(sample@assays$SCT@data[lymphoid, ], na.rm = TRUE)
  if (all(names(x = mean.exp) == rownames(x = sample@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'sample@meta.data':\n", 
        "adding gene set mean expression values in 'sample@meta.data$lymphoid.score'")
    sample@meta.data$lymphoid.score <- mean.exp
  }
  
  
  tiff(paste0(writeTo, "immune_quantification/lymphoid_signature/01_full_signature_SCT_featurePlot.tiff"), width=700, height=800)
  print(FeaturePlot(sample, features = "lymphoid.score", raster=F, reduction = "cell_coord"))
  dev.off()
  
  tiff(paste0(writeTo, "immune_quantification/lymphoid_signature/02_full_signature_SCT_vlnPlot.tiff"), width=500, height=500)
  print(VlnPlot(sample, "lymphoid.score", group.by="classification", pt.size=0, raster=F))
  dev.off()
  
  tiff(paste0(writeTo, "immune_quantification/lymphoid_signature/03_full_signature_SCT_scatter_EpCAM.tiff"), width=600, height=500)
  print(FeatureScatter(sample, "lymphoid.score", "EpCAM", plot.cor = F))
  dev.off()
  
  sample$lymphoid <- NA
  sample$lymphoid[sample$lymphoid.score > 6.5] <- "lymphoid"
  
  tiff(paste0(writeTo, "immune_quantification/lymphoid_signature/04_full_signature_SCT_lymphoid_defined.tiff"), width=700, height=800)
  print(DimPlot(sample, group.by="lymphoid", reduction = "cell_coord", 
                pt.size=.5, raster=F, shuffle=T, cols="blue", na.value="lightgrey"))
  dev.off()
  
  # lymphoid density
  coord_lymphoid <- coord[!is.na(sample$lymphoid),]
  tiff(paste0(writeTo, "immune_quantification/lymphoid_signature/05_full_signature_SCT_lymphoid_density.tiff"), width=600, height=800)
  print(ggplot(coord_lymphoid, aes(coord_1, coord_2)) + 
          geom_density_2d_filled() + 
          theme_minimal() + 
          theme(panel.grid = element_blank(), 
                axis.line = element_line(), 
                axis.text = element_blank()) + 
          xlim(x_min, x_max) + ylim(y_min, y_max))
  dev.off()
  
  lymphoid_table <- as.data.frame.matrix(table(sample$lymphoid, sample$classification, useNA="always"))
  
  # lymphoid in tumor
  lymphoid_table[1,2] / (lymphoid_table[1,2] + lymphoid_table[2,2])
  # lymphoid in non-tumor
  lymphoid_table[1,1] / (lymphoid_table[1,1] + lymphoid_table[2,1])
  
  
  saveRDS(sample, file = paste0("sample_", i, "/Sample-", i, "_Seurat_raw.rds"))
  
  
  
  
  
  
  # quantification based on myeloid cell signature
  dir.create(paste0(writeTo, "immune_quantification/myeloid_signature/"), r=T)
  
  myeloid_sign <- c("CD1c", "CD11b", "CD13", "CD14", "CD66b", "CD68", "CD163", "CD204")
  
  # select proteins of interest
  myeloid <- rownames(sample@assays$SCT[rownames(sample@assays$SCT) %in% myeloid_sign])
  # get mean expression of proteins of interest per cell
  mean.exp <- colMeans(sample@assays$SCT@data[myeloid, ], na.rm = TRUE)
  if (all(names(x = mean.exp) == rownames(x = sample@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'sample@meta.data':\n", 
        "adding gene set mean expression values in 'sample@meta.data$myeloid.score'")
    sample@meta.data$myeloid.score <- mean.exp
  }
  
  
  tiff(paste0(writeTo, "immune_quantification/myeloid_signature/01_full_signature_SCT_featurePlot.tiff"), width=700, height=800)
  print(FeaturePlot(sample, features = "myeloid.score", raster=F, reduction = "cell_coord"))
  dev.off()
  
  tiff(paste0(writeTo, "immune_quantification/myeloid_signature/02_full_signature_SCT_vlnPlot.tiff"), width=500, height=500)
  print(VlnPlot(sample, "myeloid.score", group.by="classification", pt.size=0, raster=F))
  dev.off()
  
  tiff(paste0(writeTo, "immune_quantification/myeloid_signature/03_full_signature_SCT_scatter_EpCAM.tiff"), width=600, height=500)
  print(FeatureScatter(sample, "myeloid.score", "EpCAM", plot.cor = F))
  dev.off()
  
  sample$myeloid <- NA
  sample$myeloid[sample$myeloid.score > 6.1] <- "myeloid"
  
  tiff(paste0(writeTo, "immune_quantification/myeloid_signature/04_full_signature_SCT_myeloid_defined.tiff"), width=700, height=800)
  print(DimPlot(sample, group.by="myeloid", reduction = "cell_coord", 
                pt.size=.5, raster=F, shuffle=T, cols="blue", na.value="lightgrey"))
  dev.off()
  
  # myeloid density
  coord_myeloid <- coord[!is.na(sample$myeloid),]
  tiff(paste0(writeTo, "immune_quantification/myeloid_signature/05_full_signature_SCT_myeloid_density.tiff"), width=600, height=800)
  print(ggplot(coord_myeloid, aes(coord_1, coord_2)) + 
          geom_density_2d_filled() + 
          #geom_density_2d() + 
          theme_minimal() + 
          theme(panel.grid = element_blank(), 
                axis.line = element_line(), 
                axis.text = element_blank(), 
                legend.position = "none") + # add legend manually with a continuous rectangle of the respective color scheme (low -> high density)
          xlim(x_min, x_max) + ylim(y_min, y_max))
  dev.off()
  
  myeloid_table <- as.data.frame.matrix(table(sample$myeloid, sample$classification, useNA="always"))
  
  # myeloid in tumor
  myeloid_table[1,2] / (myeloid_table[1,2] + myeloid_table[2,2])
  # myeloid in non-tumor
  myeloid_table[1,1] / (myeloid_table[1,1] + myeloid_table[2,1])
  
  
  saveRDS(sample, file = paste0("sample_", i, "/Sample-", i, "_Seurat_raw.rds"))
  
  
  
  
  
  
  
  # quantification of CD204+ cells
  
  source("marker_gating_function.R")

  sample <- marker_gating(sample, "CD204", gate_against = "EpCAM", save_plots=T)
  
  tiff(paste0(writeTo, "marker_gating/CD204_DimPlot_defined.tiff"), width=700, height=800)
  print(DimPlot(sample, group.by="CD204_pos", reduction = "cell_coord", 
                pt.size=.5, raster=F, shuffle=T, cols="blue", na.value="lightgrey"))
  dev.off()
  
  # CD204+ density
  coord_CD204 <- coord[!is.na(sample$CD204_pos),]
  tiff(paste0(writeTo, "marker_gating/CD204_DimPlot_density.tiff"), width=600, height=800)
  print(ggplot(coord_CD204, aes(coord_1, coord_2)) + 
          geom_density_2d_filled() + 
          theme_minimal() + 
          theme(panel.grid = element_blank(), 
                axis.line = element_line(), 
                axis.text = element_blank()) + 
          xlim(x_min, x_max) + ylim(y_min, y_max))
  dev.off()
  
  CD204_table <- as.data.frame.matrix(table(sample$CD204_pos, sample$classification, useNA="always"))
  
  # CD204 in tumor
  CD204_table[1,2] / (CD204_table[1,2] + CD204_table[2,2])
  # CD204 in non-tumor
  CD204_table[1,1] / (CD204_table[1,1] + CD204_table[2,1])
  
  saveRDS(sample, file = paste0("sample_", i, "/Sample-", i, "_Seurat_raw.rds"))

  
}


