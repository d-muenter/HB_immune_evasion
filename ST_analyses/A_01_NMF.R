#.libPaths('~/R/library/')
library(Seurat)
library(tidyverse)
library(future)
library(writexl)

setwd("")
source('seurat_functions_public.R')

sample_list <- c("286", "287", "288", "289", "386", "387", "388", "389", "391", "392")

for (i in sample_list){
  
  writeTo <- (paste0('Output/H-', i, '/'))
  dir.create(writeTo)
  
  # load spatial data
  load(paste0('Data/H-', i, '/H-', i, '_ST.RData'))
  
  # NMF
  # run NMF
  data = as.matrix(GetAssayData(st, assay = 'SCT', slot = 'scale.data'))
  data = data[SpatiallyVariableFeatures(st),]
  data[data < 0] = 0  # set negative values to zero
  data = data[apply(data, 1, var) > 0,]   # subsetting on genes with a variance above zero
  print(dim(data))
  range = 2:25    # apply nsNMF for ranks between 2 and 25
  #plan(multicore)
  plan(multisession)
  res.list = lapply(range, function(r){
    nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
  })
  names(res.list) = range
  dir.create(paste0('Data/H-', i))
  save(res.list, file = paste0('Data/H-', i, '/res.list_2-25.RData'))
  plan(sequential)
  # select rank
  modules.list = lapply(res.list, NMFToModules, gmin = 5) 
  print(sapply(modules.list, length))
  comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
  comp  
  mi = min(comp)
  r = names(which(comp == mi))
  r = r[length(r)]
  print(r)
  res = res.list[[r]]
  save(res, file = paste0('Data/H-', i, '/res.RData'))
  # process output
  tryCatch(expr = {
    modules = NMFToModules(res, gmin = 5)
    scores = basis(res)
    colnames(scores) = names(modules)
    coefs = coefficients(res)
    rownames(coefs) = names(modules)
    # order modules
    h = Heatmap(coefs, clustering_distance_columns = "euclidean")
    o = row_order(h)
    scores = scores[, o]
    coefs = coefs[o, ]
    modules = modules[o]
    print(modules)
    st = AddMetaData(st, t(coefs), col.name = rownames(coefs))
    # cluster NMF output
    h = Heatmap(coefs, clustering_distance_columns = "euclidean")
    hcl = as.hclust(column_dend(h))
    sig = cutree(hcl, k = length(modules))
    nmf = c(by(t(coefs), INDICES = sig, FUN = function(x){names(modules)[which.max(colMeans(x))]}))[sig]
    st$nmf = factor(nmf, levels = names(modules))
    col_nmf = c(brewer.pal(12, "Set3"), brewer.pal(8, "Set1"))[1:nlevels(st$nmf)]
    names(col_nmf) = levels(st$nmf)
  }, error = function(e){c()})
  save(st, file = paste0('Data/H-', i, '/H-', i, '_ST_postNMF.RData')) # save new seurat object with nmf modules
  
  # saving the modules to an excel file
  max_length <- max(sapply(modules, length))
  modulesSave <- lapply(modules, function(vec){
    if(length(vec) < max_length){
      c(vec, rep(NA, max_length - length(vec)))
    } else {
      vec
    }
  })
  data_frame <- data.frame(modulesSave)
  write_xlsx(data_frame, path=paste0(writeTo, "modules_", i, ".xlsx"))
  
  
  # plotting
  
  # get gene enrichment database
  gsc <- getMsigdb()
  db <- subsetCollection(gsc, "c5", c("GO:BP", "GO:CC", "GO:MF"))
  
  col_cluster = hue_pal()(nlevels(st$cluster)) # saving the cluster colors right away
  names(col_cluster) = levels(st$cluster)
  
  pdf(paste0(writeTo, 'dimplots_', i, '.pdf'), height = 15, width = 15)
  for (reduction in c('pca','umap')){
    h = DimPlot(st, pt.size = 2, reduction = reduction, group.by = 'cluster', cols = col_cluster, label = TRUE, label.size=6)
    print(h)
    h = DimPlot(st, pt.size = 2, reduction = reduction, group.by = 'nmf', cols = col_nmf)
    print(h)
    h = FeaturePlot(st, reduction = reduction, features = names(modules))
    print(h)
  }
  h = SpatialDimPlot(st, pt.size = 1.5, group.by = 'cluster', cols = col_cluster)
  print(h)
  h = SpatialDimPlot(st, pt.size = 1.5, group.by = 'nmf', cols = col_nmf)
  print(h)
  h = SpatialFeaturePlot(st, features = names(modules))
  print(h)
  dev.off()
  
  pdf(paste0(writeTo, 'heatmaps_', i, '.pdf'), height = 20, width = 10)
  sink(paste0(writeTo, 'heatmaps_', i, '.txt'))
  # Clusters
  markers_clusters = FindAllMarkers2(st, do.print = TRUE, do.plot = TRUE, enrichment.type = 'GO', group.by = 'cluster', cols = col_cluster, print.bar = FALSE, do.enrichment=F)
  # NMF
  markers_nmf = FindAllMarkers2(st, do.print = TRUE, do.plot = TRUE, enrichment.type = 'GO', group.by = 'nmf', cols = col_nmf, print.bar = FALSE, do.enrichment=F)
  # Modules
  top_ann = HeatmapAnnotation(df = data.frame('cluster' = st$cluster,
                                              'nmf' = st$nmf), 
                              col = list('cluster' = col_cluster,
                                         'nmf' = col_nmf), which = 'column')
  h = Heatmap(coefs, name = 'coefs', top_ann = top_ann,
              show_row_names = TRUE, 
              cluster_rows = FALSE,
              cluster_columns = FALSE, column_order = order(st$nmf, st$cluster),
              breaks = c(0,max(coefs)/2), colors = c('white','red'))
  print(h)
  h = Heatmap(scores[unlist(modules),], name = 'scores', 
              show_row_names = TRUE, row_names_gp = gpar(cex = 0.5),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              split = factor(unlist(mapply(rep, names(modules), sapply(modules, length))), levels = names(modules)))
  print(h)
  print(modules)
  enrichment.matrix = BuildEnrichmentMatrix(rownames(st), type = 'GO', db=geneIds(db))
  go_10 = sapply(modules, function(tbl){
    e = Enrichment(geneset = tbl, enrichment.matrix = enrichment.matrix)
    names(sort(e))[1:10]
  })
  colnames(go_10) = names(modules)
  print(go_10)
  sink()
  dev.off()
  
  
}
