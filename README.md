# Upregulated IL-1 Gene Expression in Epithelialized Tunnels of Hidradenitis Suppurativa

The repository contains codes used to process the data used in the publication and generate all plots. These codes only use existing data packages and no new method/package has been developed. The codes are well commented. A brief disciption of each file is provided below:


1) Visium_Spatial_cluster_analysis.rmd - This pipeline processes and analyzes visiusm spatial transcriptomics data from multiple samples using Seurat in R. It covers data import, metadata integration, quality control, clustering, marker analysis, and generates visualizations and marker tables for downstream analysis.

2) Single cell workflow for RCTD reference data used.rmd - This pipeline performs integration of single-cell RNA-seq data from multiple studies and samples using Seurat. It includes steps like data loading, quality control, doublet removal, normalization, and batch correction with Harmony, followed by clustering and identification of marker genes. The resulting annotated clusters are then used as reference data for deconvolution of Visium Spatial transcriptomics datasets.

3) Run_RCTD_full_mode.Rmd - This pipeline runs cell type deconvolution for Visium spatial transcriptomics samples using RCTD with an integrated single-cell reference. It loads sample metadata, processes each sample through quality control and annotation, runs RCTD for cell-type assignment, and merges the results into a combined Seurat object for downstream analysis and visualization.

4) Pseudobulk_gene_expression_Hs_Visium.Rmd - This pipeline performs pseudobulk spatial transcriptomics analysis on HS Visium data: it aggregates gene expression by region and treatment, filters low-quality samples, prepares metadata, carries out normalization and differential expression analysis (using edgeR and limma), and exports lists of significant up- and down-regulated genes for publication-ready downstream analyses.

5) Over_representation_analysis.Rmd  - This pipeline performs GO-based pathway enrichment (over-representation) analysis on significant DEGs, using the clusterProfiler and related R packages. It identifies biological processes enriched in your gene set, saves results, and visualizes the top pathways with a dotplot for downstream interpretation.

6) Geneset_enrichment_analysis.Rmd - This pipeline performs gene set enrichment analysis (GSEA) on differential expression results using GO Biological Process terms, identifying pathways enriched in your data. It saves the full GSEA output and visualizes the top pathways in dotplots, supporting interpretative reporting for spatial transcriptomics studies.

7) Jaccard_index_function.R - Calculates Jaccard similarity between two users (columns) in a binary matrix

8) Gene_specific_plots_for_all_sections.Rmd - This pipeline visualizes gene expression across HS skin Visium sections using Seurat in R. For each sample, it displays the H&E image, annotated tissue regions, and a spatial feature plot of gene expression, enabling comparative visualization of gene localization within annotated skin regions.

9) Correaltion_gene_expression_and_RCTD_cell_types.R - This workflow calculates per-sample Pearson correlations between the abundance of each cell type and cytokine expression (TNF, IL1A, IL1B, IL17A, IL17F, IL23A) in spatial transcriptomics HS Visium data. It outputs correlation matrices for all cytokines across all samples and cell types to CSV, enabling downstream analysis of cellular-cytokine relationships.

10) Bulk_RNAseq_Analysis.Rmd - This bulk RNA-seq pipeline loads and organizes count data with metadata, conducts quality checks, and performs normalization and filtering. Differential expression analysis via limma/voom yields gene lists for key contrasts, with PCA plots to visualize sample clustering across conditions.






