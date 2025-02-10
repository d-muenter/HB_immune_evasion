library(Seurat)
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)

setwd("")

load("Data/merged_integration_ST_postSCIntegration.RData")
DefaultAssay(dataIntegrated) <- "SCT"

writeTo <- "Output/merged_integration/neighborhood_proximity/"
dir.create(writeTo)


### calculate distances & neighbors
all_neighbors <- c()
all_distances <- list()
images <- names(dataIntegrated@images)
for (i in 1:length(images)){
  
  imageData <- dataIntegrated@images[[i]]
  coord = imageData@coordinates[,c("imagerow", "imagecol")] # get spot pixel coordinates
  
  distances = as.matrix(dist(coord)) # calculate interspot distances based on pixel coordinates
  distances <- distances/sort(distances[distances > 0], decreasing = FALSE)[1] # divide by lowest non-zero interspot distance (= distance between two neighboring spots) to scale such that the unit is interspot distance
  
  all_distances[[i]] <- distances
  names(all_distances)[i] <- images[i]
  
  neighbors = apply(distances, 1, function(d){ # find neighbors: all spots within each row that have a distance < 1.5 are a neighbor to the spot of that row
    names(d)[d >= 0 & d <= 1.5]
  })
  
  print(paste(images[i], sum(names(neighbors) != rownames(coord)) == 0))
  
  all_neighbors <- c(all_neighbors, neighbors)

}


### calculate neighborhood scores
predictions <- t(as.matrix(GetAssayData(dataIntegrated, slot = "data", assay = "predictions_overview")))
predictions <- predictions[,1:5]  # get the celltype estimates; rows are spots, columns are cell types

# setup a matrix for the neighborhood scores; rows will be spots, columns will be cell types -> every spot will have a neighborhood score for each cell type
coef_nei <- matrix(nrow = nrow(predictions), ncol = ncol(predictions))
rownames(coef_nei) <- rownames(predictions)
colnames(coef_nei) <- colnames(predictions)

for (i in colnames(predictions)){  # for each different cell type
  
  for (j in names(all_neighbors)){  # and for each individual spot
    
    coef_nei[j, i] <- mean(predictions[all_neighbors[[j]], i])   # calculate the mean abundance estimate between all neighbor spots for that respective spot and cell type
                                                                          
  }
}

saveRDS(all_neighbors, file = paste0(writeTo, "neighbors.rds"))
saveRDS(coef_nei, file = paste0(writeTo, "neighbor_scores.rds"))
saveRDS(all_distances, file = paste0(writeTo, "distances.rds"))



### calculate proximity scores (based on previously calculated distances)

pred_sample <- list()  # get the celltype estimates - separately for each sample; rows are spots, columns are cell types
for (i in unique(dataIntegrated$orig.ident)){
  dataTemp <- subset(dataIntegrated, subset = orig.ident == i)
  pred_sample[[i]] <- t(as.matrix(GetAssayData(dataTemp, slot = "data", assay = "predictions_overview")))
  pred_sample[[i]] <- pred_sample[[i]][,1:5]
}


source("findProximityThresholds.R")
thresholds <- findProximityThresholds(dataIntegrated, assay="predictions_overview")

all_coef_prox <- list() # create a list. Then for each image/sample create a matrix within that list which has the same number of columns as there are cell types and the number of rows as there are spots in that sample
                        # then in the end rbind all of these matrices -> one large matrix with has the proximity score for each cell type and each spot within each sample

for (i in 1:length(images)){  # for each different sample
  
  coef_prox <- matrix(ncol = ncol(pred_sample[[i]]), 
                      nrow = nrow(pred_sample[[i]]))
  rownames(coef_prox) <- rownames(pred_sample[[i]])
  colnames(coef_prox) <- colnames(pred_sample[[i]])
  
  for (j in colnames(pred_sample[[i]])){  # and for each different cell type
    
    for (k in rownames(all_distances[[i]])){  # and for each individual spot
      
      # get names of spots with that celltype over the threshold
      spots <- rownames(pred_sample[[i]][pred_sample[[i]][,j] > thresholds[j],])

      # of all those potential spots, find the spot with the smallest distance to the current spot of interest
      if (length(spots) >= 1) { # in case there are no spots with that cell type over the threshold in that sample
        # take min spot
        mi <- min(all_distances[[i]][k, spots])
        
        prox <- 1/(1+mi)
        coef_prox[k, j] <- prox
        
      }
      
    }
    
  }
  
  all_coef_prox[[i]] <- coef_prox
  names(all_coef_prox)[i] <- images[i]
  
  print(paste0(images[i], ": Done"))
  
}

all_coef_prox_df <- do.call(rbind, all_coef_prox)

saveRDS(all_coef_prox_df, file = paste0(writeTo, "proximity_scores_autoTH_min1.rds"))


