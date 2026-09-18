# Single-cell Transcriptomic Analysis of OMD-associated Fibroblasts in Oral Submucous Fibrosis

Analysis code for the single-cell transcriptomic analysis presented in **Figure 2** of the manuscript:

**Danshensu attenuates Osteomodulin-driven fibrogenesis in oral submucous fibrosis**

## Description

This repository contains the R code used to analyze single-cell RNA sequencing (scRNA-seq) data from oral submucous fibrosis (OSF) and healthy oral tissues.

The analysis focuses on the cellular distribution of **osteomodulin (OMD)** and the transcriptional characteristics of **OMD-positive fibroblasts**, with particular emphasis on extracellular matrix remodeling and TGF-β-related signaling programs.

The scRNA-seq dataset used in this analysis is available from the **Genome Sequence Archive for Human (GSA-Human)** under accession number:

**HRA004340**

## Key analyses

The `Fig2.R` script includes:

* UMAP visualization of OMD-positive and OMD-negative fibroblasts
* Differential expression analysis between OMD-positive and OMD-negative fibroblasts
* Volcano plot visualization of differentially expressed genes
* Gene Ontology (GO) enrichment analysis
* Gene Set Enrichment Analysis (GSEA)
* Analysis of TGF-β-related signaling programs
* Intersection of genes from TGF-β-related GO terms
* Heatmap visualization of representative differentially expressed genes
* Correlation analysis between OMD and collagen-associated genes

These analyses were used to characterize the transcriptional state associated with OMD expression in fibroblasts and to investigate its relationship with extracellular matrix remodeling and TGF-β signaling.

## Files

```text
Fig2.R
```

`Fig2.R` contains the analysis and visualization code used for the single-cell transcriptomic results shown in Figure 2.



## Requirements

The analyses were performed in **R**.

Major R packages include:

```r
Seurat
tidyverse
dplyr
ggplot2
ggpubr
ggstatsplot
ggprism
gground
ggvenn
ggnewscale
patchwork
clusterProfiler
enrichplot
org.Hs.eg.db
EnhancedVolcano
ComplexHeatmap
circlize
```

## Data availability

The single-cell RNA sequencing data analyzed in this study are available from the **Genome Sequence Archive for Human (GSA-Human)**.

**Accession: HRA004340**

Access to the human sequencing data is subject to the applicable data-access policies of GSA-Human.
