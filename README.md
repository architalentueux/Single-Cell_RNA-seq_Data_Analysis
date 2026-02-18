# Single-Cell_RNA-seq_Data_Analysis
I performed a single-cell RNA sequencing (scRNA-seq) analysis using the publicly available dataset from GEO (Using R and Python)
The scientific article entitled “A single cell RNAseq benchmark experiment embedding ‘controlled’ cancer heterogeneity” describes the creation and publication of a single-cell RNA-seq dataset designed as a benchmark resource to evaluate and develop bioinformatics methods for analyzing cancer cell heterogeneity.

The authors present a single-cell RNA sequencing experiment performed using 10x Genomics technology on a panel of seven lung cancer cell lines (including A549, PC9, HCC78, etc.), each characterized by distinct driver mutations (EGFR, KRAS, BRAF, ALK, MET, ERBB2, ROS1).

The objective of the study was to create a controlled experimental framework in which cell composition is known, while still mimicking real-world tumor heterogeneity. This enables benchmarking and comparison of scRNA-seq analysis tools. By generating mixtures of cells from different cell lines, the dataset allows researchers to evaluate how well computational methods can:

Detect subpopulations

Classify cells according to their origin

Identify transcriptional profiles associated with specific mutations

Overall, this dataset serves as a standardized reference resource for the scientific community, facilitating objective assessment of bioinformatics methods in oncology.

# Describing My Single-Cell Analysis

In my project, I performed a single-cell RNA sequencing (scRNA-seq) analysis using the publicly available dataset from GEO accession GSE243665, derived from the experiment described in the above study.

The dataset includes the standard 10x Genomics output files:

GSE243665_A549_barcodes.tsv.gz — cell barcode identifiers

GSE243665_A549_features.tsv.gz — gene annotation (features)

GSE243665_A549_matrix.mtx.gz — sparse gene expression count matrix

These files were imported into a single-cell analysis pipeline (e.g., using Seurat or Scanpy), where I performed quality control filtering, normalization, dimensionality reduction, clustering, and cell population characterization.

The objective of this analysis was to explore transcriptional heterogeneity within A549 cells and identify potential transcriptionally distinct subpopulations or cellular states.
