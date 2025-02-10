library(tidyverse)
library(Seurat)
library(Matrix)


for (i in 1:7){ # for each sample
  
  writeTo <- paste0("sample_", i, "/SeuratOutput/")
  dir.create(writeTo)

  
  sample_crops <- c("Crop_1", "Crop_2", "Crop_3") # for computational efficiency samples were divided into multiple regions ("crops") for cell segmentation which now need to be recombined
  
  sample_expr <- list()
  features <- list()
  
  
  for (j in sample_crops) {
    
    # read MACSiQ View output data matrix
    sample_expr[[j]] <- read.csv(paste0("sample_", i, "/sample_", i, "_", j, "/sample_", i, "_", j, "_expr.csv"))
    
    # extract feature names
    features[[j]] <- colnames(sample_expr[[j]])[grepl("Cell.Exp", colnames(sample_expr[[j]]))]
    features[[j]] <- gsub(".Cell.Exp", "", features[[j]])
    features[[j]] <- gsub("\\.", "_", features[[j]])
    features[[j]] <- gsub("Mouse_IgG1_", "", features[[j]])
    features[[j]] <- gsub("Goat_IgG_", "", features[[j]])
    
    
  }
  
  if (sum(features[[1]] != features[[2]]) == 0 & sum(features[[1]] != features[[3]]) == 0) { # check feature names from all crops are identical
    features <- features[[1]]
  }
  
  
  barcodes <- list()
  
  for (j in sample_crops) {
    
    # generate cell IDs
    IDs <- as.character(1:nrow(sample_expr[[j]]))
    long_IDs <- as.vector(sapply(IDs, function(x) {
      n <- 5 - nchar(x)
      zeros <- paste0(rep("0", n), collapse="")
      return(paste0(zeros, x))
    }))
    
    barcodes[[j]] <- paste0("Sample-", i, "-", j, "_", long_IDs)
    barcodes[[j]] <- gsub("Crop_", "Crop-", barcodes[[j]])
    
  }
  
  barcodes <- c(barcodes[[1]], barcodes[[2]], barcodes[[3]])
  
  
  # combine x and y coordinates from all crops
  coord <- rbind(sample_expr[[1]][,73:75], sample_expr[[2]][,73:75], sample_expr[[3]][,73:75])
  rownames(coord) <- barcodes
  coord$Crop <- gsub("_.*", "", rownames(coord))
  
  
  min_1 <- min(coord$Nuc.Y.Inv[coord$Crop == paste0("Sample-", i, "-Crop-1")])
  max_1 <- max(coord$Nuc.Y.Inv[coord$Crop == paste0("Sample-", i, "-Crop-1")])
  min_2 <- min(coord$Nuc.Y.Inv[coord$Crop == paste0("Sample-", i, "-Crop-2")])
  max_2 <- max(coord$Nuc.Y.Inv[coord$Crop == paste0("Sample-", i, "-Crop-2")])
  min_3 <- min(coord$Nuc.Y.Inv[coord$Crop == paste0("Sample-", i, "-Crop-3")])
  max_3 <- max(coord$Nuc.Y.Inv[coord$Crop == paste0("Sample-", i, "-Crop-3")])
  
  print(c(min_1, max_1, min_2, max_2, min_3, max_3))
  
  g1_coord <- ggplot(coord, aes(Nuc.X, Nuc.Y.Inv)) + geom_point(aes(color=Crop), alpha=0.2)
  
  
  tiff(paste0(writeTo, "coord_crop_paste.tiff"), width=600, height=800)
  print(g1_coord)
  dev.off()
  
  
  # generate data matrix combining expression, feature names and cell IDs
  data_raw <- rbind(sample_expr[["Crop_1"]][,37:72], # raw expression data was in columns 37-72 of the exported data matrix
                    sample_expr[["Crop_2"]][,37:72], 
                    sample_expr[["Crop_3"]][,37:72])
  colnames(data_raw) <- features
  rownames(data_raw) <- barcodes
  data_raw <- t(data_raw)
  
  matrix_sample <- Matrix(
    data = data_raw, 
    ncol = ncol(data_raw), 
    sparse = TRUE
  )
  
  
  # generate a Seurat object from the data matrix
  sample <- CreateSeuratObject(counts = matrix_sample, 
                               project = paste0("Sample-", i), 
                               assay="cell_exp")
  
  tiff(paste0(writeTo, "nFeature_nCount.tiff"), width=800, height=500)
  print(VlnPlot(sample, features = c("nFeature_cell_exp", "nCount_cell_exp"), ncol = 2, pt.size = 0, group.by="orig.ident"))
  dev.off()
  
  sample <- SCTransform(sample, assay = "cell_exp")
  sample <- RunPCA(sample, features = rownames(sample), npcs = 30)
  sample <- RunUMAP(sample, dims = 1:30)
  sample <- FindNeighbors(sample, reduction = "pca", dims = 1:30)
  sample <- FindClusters(sample, resolution = 0.3)
  
  tiff(paste0(writeTo, "umap_clustering.tiff"), width=700, height=500)
  print(DimPlot(sample, raster=F))
  dev.off()
  
  DefaultAssay(sample) <- "cell_exp"

  
  # add x and y coordinates
  matrix_coord <- coord_sub[,c("Nuc.X", "Nuc.Y.Inv")]
  colnames(matrix_coord) <- c("Nuc_1", "Nuc_2")
  matrix_coord <- as.matrix(matrix_coord)
  
  sample[["cell_coord"]] <- CreateDimReducObject(embeddings = matrix_coord, key = "coord_", assay = "cell_exp")
  
  tiff(paste0(writeTo, "cell_coord_clusters.tiff"), width=700, height=800)
  print(DimPlot(sample, reduction = "cell_coord", pt.size = .5, raster=F))
  dev.off()
  

  saveRDS(sample, file = paste0("sample_", i, "/Sample-", i, "_Seurat_raw.rds"))
  
  
}



