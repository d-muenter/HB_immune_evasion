# Source Code for Münter et al.  

 
This repository contains code that was used for our manuscript "Multiomic analysis uncovers a continuous spectrum of differentiation and Wnt-MDK-driven immune evasion in hepatoblastoma" and deviated from previously published workflows.  

### mIF Analysis  
Code used for the analysis of MACSima multiplex immunofluorescence cell segmentation data.  

- **01_Seurat_from_MACSima.R:** Importing cell segmentation data into R and generating respective Seurat objects. Each row in the input data corresponds to one cell, columns contain individual protein intensities and cell coordinates.  
- **02_Seurat_MACSima_analysis.R:** Classifying each cell into tumor/non-tumor based on whether it is located inside or outside a tumor cell area (and not whether it is a tumor cell itself). The input data format is identical to that of '01_Seurat_from_MACSima.R' but only contains cells in tumor cell areas based on previous gating in the MACSiQ View software. Lymphoid, myeloid and CD204+ cells get quantified according to marker expression.  
- **03_immune_quantification.R:** Comparison of the proportions of lymphoid, myeloid and CD204+ cells between tumor and non-tumor cell areas.  
- **marker_gating_function.R:** R function to gate cells based on the expression of two markers.  

### ST Analyis  
Code used for the analysis of Visium spatial transcriptomics data. NMF, neighborhood and proximity analyses heavily relied on the approach described by Barkley et al. (https://doi.org/10.1038/s41588-022-01141-9) and code was adapted only slightly for our data. Similarly, generation of metaprograms  heavily relied on the approach introduced by Gavish et al. (https://doi.org/10.1038/s41586-023-06130-4).  

- **A_01_NMF.R:** Running non-negative matrix factorization (NMF) and generating corresponding gene modules. Inputs are Seurat objects of individual spatial transcriptomics samples. 'seurat_functions_public.R' refers to a file of that name published by Barkley et al..  
- **B_01_make_genes_nmf_list.R:** Generating a data matrix containing NMF module genes for the generation of metaprograms. Inputs are the final results from running NMF.  
- **B_02a_generate_MPs.R:** Generating metaprograms based on the previously calculated NMF modules and the code published by Gavish et al.. 'robust_nmf_programs.R' refers to a file of that name published by Gavish et al..  
- **B_02b_check_MPs.R:** Filtering metaprogram genes based on the correlation of their expression to the overall metaprogram expression score in each spot of the spatial transcriptomics samples.  
- **C_02a_neighborhood_proximity.R:** Calculating cell type neighborhood and proximity scores. Input is a Seurat object containing all samples and respective cell type abundance estimates.  
- **findProximityThresholds.R:** R function to calculate threshold scores based on the cell type abundance estimates.  
