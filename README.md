# Single-cell and Spatial Transcriptomics Analysis Workflow

This repository contains the analysis code used for single-cell, single-nucleus, spatial transcriptomics, and spatial multi-omics data processing in our study of Alzheimer’s disease-associated spinal cord neuroimmune remodeling.

The workflow supports multiple data modalities, including human single-nucleus RNA sequencing, 10x Genomics Visium HD spatial transcriptomics, mouse single-cell RNA sequencing, mouse single-nucleus RNA sequencing, reference dataset integration, cell–cell communication analysis, and MiP-seq-based spatial proximity analysis.

## Overview

The repository is organized as a modular and reproducible workflow covering the following analytical steps:

1. **Single-cell and single-nucleus RNA-seq preprocessing**
   - Cell Ranger-based alignment and UMI counting
   - Quality control of cells and nuclei
   - Removal of ambient RNA, empty droplets, and doublets
   - Normalization, dimensionality reduction, clustering, and cell-type annotation

2. **Dataset integration and cell-type annotation**
   - Integration of scRNA-seq and snRNA-seq datasets
   - Reciprocal PCA-based integration
   - scVI-based integration for microglial subclustering
   - Annotation using canonical marker genes and public reference atlases

3. **Microglial and neuronal subclustering**
   - Identification of microglial transcriptional states, including homeostatic, IFN-I-responsive, cycling, DAM I, and DAM II populations
   - Neuronal subtype annotation using published spinal cord reference datasets
   - Cross-dataset comparison of brain and spinal cord microglial states

4. **Differential expression and pathway analysis**
   - Pseudobulk differential expression analysis using DESeq2
   - Cell-type-level pairwise comparisons using Seurat
   - Functional enrichment analysis using clusterProfiler and Enrichr
   - Gene-set activity scoring using AUCell and UCell

5. **Spatial transcriptomics analysis**
   - 10x Genomics Visium HD spatial data preprocessing
   - Spatial scoring of amyloid-associated gene signatures
   - Definition of Aβ-associated spatial niches
   - RCTD-based spatial deconvolution using snRNA-seq reference annotations
   - Spatial differential expression analysis using Scanpy
   - Visualization of spatial gene expression and pathway activity

6. **Cell–cell communication analysis**
   - Ligand–receptor inference using CellChat v2.2.0
   - Orthogonal validation using LIANA+
   - Comparative analysis of signaling networks between experimental groups

7. **MiP-seq spatial analysis**
   - DAPI-based nuclear segmentation using StarDist
   - Cell-level spatial expression matrix construction
   - Cell-type annotation based on marker-gene expression
   - Spatial proximity analysis between microglia and neurons

## Main software and packages

The workflow uses commonly adopted tools in single-cell and spatial omics analysis, including:

- Cell Ranger
- Seurat
- Scanpy
- scDblFinder
- CellBender
- DropletQC
- scVI
- SeuratWrappers
- DESeq2
- clusterProfiler
- Enrichr
- AUCell
- UCell
- RCTD
- CellChat
- LIANA+
- StarDist
- Fiji/ImageJ

## Scope

This repository is designed for researchers working with single-cell, single-nucleus, spatial transcriptomics, and spatial multi-omics datasets. Although the workflow was developed for spinal cord and Alzheimer’s disease-related analyses, most modules can be adapted to other tissues, disease models, and spatial transcriptomic studies.

## Notes

Raw sequencing data, protected human sample information, and large intermediate files are not included in this repository unless otherwise specified. Users should modify file paths, metadata tables, and sample information according to their own dataset structure.

---

All repository code and analyses were completed by Deng Yusen.
