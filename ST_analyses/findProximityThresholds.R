findProximityThresholds <- function(seurat_object, assay, slot = "data") {
  
  predictions <- t(as.matrix(GetAssayData(seurat_object, slot = slot, assay = assay)))
  predictions <- predictions[,-which(colnames(predictions) == "max")]  # get the celltype estimates for the respective annotation; rows are spots, columns are cell types
  
  thresholds <- sqrt(apply(log(predictions+1), 2, mean))
  
  return(thresholds)
}
