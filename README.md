# SA_bone_marrow
scRNA-seq analysis for Radhouani et al. paper, Science Immunology, 2025.

The running order of the scripts is as follows:

1. Reading the data in R using:
   * R/1_read_in_pilot.Rmd
   * R/2_read_in_final.Rmd

2. Performing single-cell integration using scVI, initial annotation using scNym, and leiden clustering of mature and hematopoietic stem cells (HSC). 
   * python/1_scvi_integration.ipynb
   * python/2_scnym_annotation.ipynb
   * python/3_clustering_for_manual_annotation.ipynb
  
3. Marker detection of the clusters.
   * R/3_hvg_DEGs_HSC_scnym_annotation.Rmd
   * R/4_hvg_DEGs_mature_scnym_annotation.Rmd
  
These markers have been used for manual annotation of the clusters.

4. Transfer and plotting of the data with manually annotated clusters. This script relies on python/4_dynamo_cell_cycle.ipynb for cell cycle annotation.
   * python/5_transfering_manual_annotation.ipynb
   
5. Differential expression between SA- and PBS-treated mice across cell types.
   * R/5_DEGs_HSC_scnym_annotation_manual_annotation.Rmd
   * R/6_DEGs_mature_scnym_annotation_manual_annotation.Rmd

[![DOI](https://zenodo.org/badge/667063786.svg)](https://doi.org/10.5281/zenodo.14634569)
