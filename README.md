# Spatially-resolved-transcriptomics-of-macrophage-heterogeneity-and-prognostic-significance-in-DLBCL


1. DSP analysis 


Differential expression analyses (DEA_DLBCL.R, DEA_Tonsil.R) were carried out using the limma package's moderated t-test.1 Upregulated/downregulated genes were selected by applying the Benjamini-Hochberg correction on the p-values (adjusted p-values < 0.05). In case of the absence of statistically significant genes, upregulated/downregulated genes were selected by applying a threshold on the absolute log2-fold-change (|log2FC| > 0.58).

The gene set enrichment analysis (function_hallmark_enrichment.R) has been manually implemented by calculating the Fisher exact test p-values on the GSEA Molecular Signature Database (MSigDB). The BH correction for multiple comparisons has been applied to evaluate the statistical significance of the pathways. An additional gene set enrichment analysis has been carried out using the GSEA software (version 4.2.3)

Kaplan-Meier method has been used to estimate the survival functions among groups, and the log-rank test has been used to test the differences in the overall survival between the selected groups (survival.R)

2. Single cell analysis 

Seurat (v2.3.0)11 was used for the analysis of the single cell datasets. All functions were run with default parameters, unless specified otherwise. Low quality cells (<200 genes/cell, <3 cells/gene and >10% mitochondrial genes) were excluded. Harmony12 was used to integrate the single cell datasets and remove any batch effects. SingleR13 was used for annotation of clusters. Refinement and validation of annotation was conducted by projecting and evaluating a curated B cell, T cell and macrophage marker list (scRNASeq_code.Rmd)

For each macrophage signature (GC- and IF like, LZ- and DZ-like, GC- and DLBCL-like, Non-relapse- and relapse-like), the top 50 or all genes were taken to create a Module Score to analyze feature expression. This score was then projected onto the uniform manifold approximation and projection (UMAP) space named “Mo-MacVERSE”. The region that showed enrichment was further subsetted and analyzed. (Signatures_Plot.R, Signatures_Plot_List.R)

### Authors
Shruti Sridhar: [shruti.sridhar@u.nus.edu](mailto::shruti.sridhar@u.nus.edu). 
Giorgio Bertolazzi: [giorgio.bertolaz@gmail.com]((mailto::giorgio.bertolaz@gmail.com)
