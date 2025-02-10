try(suppressMessages(library(tidyverse)))
try(suppressMessages(library(Seurat)))
try(suppressMessages(library(gridExtra)))


marker_gating <- function(seurat_obj, 
                          marker_to_gate,
                          gate_against = "b-Catenin", 
                          gate = NULL, 
                          save_plots = FALSE, 
                          save_dir = paste0(writeTo, "marker_gating/")) {
  
  if(save_plots==T) dir.create(save_dir, r=T)
  if(sum(marker_to_gate %in% rownames(seurat_obj)) != length(marker_to_gate)){
    stop(paste0(marker_to_gate[!(marker_to_gate %in% rownames(seurat_obj))], " => Not a valid marker name."))
  }
  
  for (marker in marker_to_gate) {
    
    if (paste0(marker, "_pos") %in% colnames(seurat_obj@meta.data)) {
      overwrite <- readline(prompt = paste0("Gate for the chosen marker (", marker, ") seems to already exist. Overwrite? (Y/N) "))
      if (overwrite %in% c("n", "N")) {
        stop() } else {
          seurat_obj@meta.data[[marker]] <- NULL
        }
    }
    
    marker_exp <- as.vector(seurat_obj@assays$cell_exp[marker,])
    
    scatter_plot <- FeatureScatter(seurat_obj, marker, gate_against, group.by="classification", shuffle=T, plot.cor=F)
    print(scatter_plot)
    if(save_plots==T){
      tiff(paste0(save_dir, marker, "_", gate_against, "_scatter.tiff"), width=600, height=500)
      print(scatter_plot)
      dev.off()
    }
    
    x = readline(prompt = "Enter threshold intensity value or s to skip marker: ")
    if(x == "s") next
    gate = as.numeric(x)
    
    marker_pos <- colnames(seurat_obj)[marker_exp > gate]
    
    seurat_obj[[paste0(marker, "_pos")]] <- NA
    seurat_obj[[paste0(marker, "_pos")]][rownames(seurat_obj@meta.data) %in% marker_pos,] <- marker
    
    umap <- DimPlot(seurat_obj, group.by=paste0(marker, "_pos"), raster=F, shuffle=T, cols="blue", na.value="lightgrey")
    
    coord_plot <- DimPlot(seurat_obj, reduction = "cell_coord", group.by = paste0(marker, "_pos"), pt.size=.5, 
                          raster=F, shuffle=T, cols="blue", na.value="lightgrey") + 
      annotate(geom="text", label=paste0("Threshold: ", gate), x=-Inf, y=Inf, hjust=-0.1, vjust=1)
    
    print(grid.arrange(umap, coord_plot, ncol=2))
    if(save_plots==T){
      tiff(paste0(save_dir, marker, "_DimPlot.tiff"), width=1400, height=700)
      print(grid.arrange(umap, coord_plot, ncol=2))
      dev.off()
    }
    
    if(marker != marker_to_gate[length(marker_to_gate)]){
      readline(prompt = "Press [enter] to continue to next marker.")
      
    }
    
  }
  
  return(seurat_obj)
  
}
