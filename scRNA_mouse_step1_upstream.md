# scRNA / snRNA Mouse — Step 1: Upstream Analysis

Quality control, doublet removal, integration, clustering, cell-type
annotation, and microglial subclustering for the WT vs 5XFAD mouse spinal
cord (and additional brain comparison) datasets.

---

## 0. Method Summary

- **Alignment**: Cell Ranger v8.0.1; mouse GRCm39 + Ensembl release 113.
- **scRNA-seq QC**: ≥ 300 detected genes, ≥ 1,000 UMIs, ≤ 25,000 UMIs, mt% ≤ 10%.
- **snRNA-seq QC**: ≥ 300 detected genes, ≥ 1,000 UMIs, ≤ 20,000 UMIs, mt% ≤ 2%.
- **Doublet removal**: scDblFinder.
- **Integration**: Seurat v5; counts split per `orig.ident`; LogNormalize;
  top 3,000 highly variable features; reciprocal-PCA (rPCA) integration;
  first 20 integrated dimensions used for SNN graph and UMAP.
- **Annotation**: canonical markers cross-referenced against CellMarker2.0.
- **Microglial subclustering**: scVI (via SeuratWrappers); 30-D latent space;
  clustering resolution 0.3.

---

## 1. Per-Sample Loading + Doublet Removal

A shared helper reads the Cell Ranger `filtered_feature_bc_matrix`, runs
scDblFinder with the 10x heuristic doublet rate, drops doublets, and
attaches sample-level metadata. The same helper is used for scRNA-seq and
snRNA-seq libraries.

```r
library(Seurat)
library(DropletUtils)
library(scDblFinder)
library(SeuratWrappers)

no_double <- function(samplename) {
  countpath <- file.path("xx/Cellranger_output", samplename,
                         "outs/filtered_feature_bc_matrix")
  sce <- DropletUtils::read10xCounts(countpath)
  sce$Sample <- samplename
  colnames(sce) <- sce$Barcode

  doubrate <- (ncol(sce) * 0.008) / 1000  # 10x empirical doublet rate
  sce <- scDblFinder(sce, dbr = doubrate)

  counts <- Read10X(data.dir = countpath)
  obj <- CreateSeuratObject(counts    = counts,
                            meta.data = as.data.frame(colData(sce)))
  obj <- subset(obj, scDblFinder.class == "singlet")
  obj <- RenameCells(obj, add.cell.id = samplename)
  obj$orig.ident <- samplename
  obj
}

sample_info <- function(obj, info1, info2, info3) {
  obj$Batch1 <- info1   # AD / WT
  obj$Batch2 <- info2   # detailed group
  obj$tissue <- info3   # brain / spinal cord
  obj
}
```

### 1.1 scRNA-seq (spinal cord; 5XFAD vs WT)

```r
setwd("/data/dengys/yhc_pain_sp/DZOE2024110757-b1/out")
sampleinfo <- read.csv("sample_info_pain_new.csv", header = TRUE)

scRNAlist <- list()
for (i in seq_len(nrow(sampleinfo))) {
  obj <- no_double(sampleinfo[i, 1])
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- sample_info(obj, sampleinfo[i, 2], sampleinfo[i, 3], sampleinfo[i, 4])
  scRNAlist[[i]] <- obj
}
scRNA_AD_SC <- merge(scRNAlist[[1]], y = scRNAlist[-1])
scRNA_AD_SC <- JoinLayers(scRNA_AD_SC)
save(scRNA_AD_SC, file = "RAW_AD_SC.rdata")
```

### 1.2 snRNA-seq

snRNA-seq libraries are processed with the same `no_double` /
`sample_info` helpers; only the sample sheet and Cell Ranger output paths
differ. The mitochondrial threshold is tightened to 2% at the joint-QC
step (Section 2) to reflect the lower cytoplasmic contamination expected
from nuclei.

---

## 2. Merge and Joint Quality Control

scRNA-seq and snRNA-seq Seurat objects are merged into one object before
QC filtering, so that downstream filtering, normalization, and
integration are all performed on the combined dataset.

```r
objs <- list(scRNA_AD_brain = scRNA_AD_brain,
             scRNA_AD_SC    = scRNA_AD_SC,
             scRNA_brain    = scRNA_brain)
total <- merge(objs[[1]], objs[-1])
total <- JoinLayers(total)
Idents(total) <- "orig.ident"

# Pre-QC visualization
pdf("sampleinfo1.pdf", width = 16, height = 9)
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0, ncol = 2)
dev.off()

# scRNA-seq cohort: nCount ≤ 25,000 and mt < 10%
total <- subset(total, subset =
                  nFeature_RNA > 300 &
                  nCount_RNA   > 1000 &
                  nCount_RNA   < 25000 &
                  percent.mt   < 10)

# snRNA-seq cohort: identical code, thresholds nCount ≤ 20,000 and mt < 2%
# total <- subset(total, subset =
#                   nFeature_RNA > 300 &
#                   nCount_RNA   > 1000 &
#                   nCount_RNA   < 20000 &
#                   percent.mt   < 2)

pdf("sampleinfo_after_flit.pdf", width = 16, height = 9)
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0, ncol = 2)
dev.off()
save(total, file = "RAW_total.rdata")
```

The QC violin plots include the four metrics
`nCount_RNA / nFeature_RNA / percent.mt / scDblFinder.score`, all of
which are already attached to `total` after Sections 1–2.

---

## 3. rPCA Integration, Clustering, and 13-Cell-Type Annotation

```r
# Seurat v5 layered representation
total[["RNA"]] <- split(total[["RNA"]], f = total$orig.ident)
total <- NormalizeData(total, normalization.method = "LogNormalize",
                       scale.factor = 10000)
total <- FindVariableFeatures(total, selection.method = "vst", nfeatures = 3000)
total <- ScaleData(total)
total <- RunPCA(total)

# Reciprocal-PCA integration
total <- IntegrateLayers(
  object         = total,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  verbose        = FALSE
)

# Clustering and UMAP
total <- FindNeighbors(total, reduction = "integrated.rpca", dims = 1:30)
total <- FindClusters (total, resolution = 0.1, cluster.name = "rpca_clusters")
total <- RunUMAP     (total, reduction = "integrated.rpca", dims = 1:30,
                       reduction.name = "umap.rpca", min.dist = 0.2)
total <- JoinLayers(total)
```

### 3.1 Cluster-to-Cell-Type Mapping (13 populations)

```r
celltype <- data.frame(ClusterID = 0:33, celltype = 0:33)

celltype[celltype$ClusterID %in% c(0,1,3,17,29,22,33),            2] <- "Microglia"
celltype[celltype$ClusterID %in% c(15),                            2] <- "Peri"
celltype[celltype$ClusterID %in% c(5),                             2] <- "Endo"
celltype[celltype$ClusterID %in% c(8,11,14,18,19,20,24,26,28,31), 2] <- "Neuron"
celltype[celltype$ClusterID %in% c(21),                            2] <- "Ependymal"
celltype[celltype$ClusterID %in% c(23),                            2] <- "BAM"
celltype[celltype$ClusterID %in% c(2,4,6,7,9,10,12),               2] <- "Oligo"
celltype[celltype$ClusterID %in% c(13),                            2] <- "Astro"
celltype[celltype$ClusterID %in% c(16),                            2] <- "OPC"
celltype[celltype$ClusterID %in% c(27),                            2] <- "Fibro"
celltype[celltype$ClusterID %in% c(30),                            2] <- "Cycling"
celltype[celltype$ClusterID %in% c(25),                            2] <- "T/NK"
celltype[celltype$ClusterID %in% c(32),                            2] <- "Mono"

total$celltype <- "NA"
for (i in seq_len(nrow(celltype))) {
  total@meta.data[total$seurat_clusters == celltype$ClusterID[i], "celltype"] <-
    celltype$celltype[i]
}

ctypes_order <- c("Microglia","Neuron","Oligo","Astro","OPC","Endo",
                  "Ependymal","BAM","Mono","Cycling","T/NK","Peri","Fibro")
total$celltype <- factor(total$celltype, levels = ctypes_order)

sc_color <- c("#9F987F","#D5B26C","#A499CC","#E07370","#4C8588","#56642a",
              "#A56BA7","#5E4D9A","#B6766C","#8A8CA7","#6D7A57","#c6a5c1","#941456")

DimPlot(total, group.by = "celltype",
        cols = sc_color, pt.size = 0.5) + NoAxes()
```

### 3.2 Joint sc/sn Sample Label

```r
clean_ident <- gsub("_SC$", "", total$orig.ident)
total$Batch3 <- dplyr::case_when(
  grepl("^AD_10M_MG", clean_ident) ~ paste0("scRNA_AD", gsub("^AD_10M_MG", "", clean_ident)),
  grepl("^WT_10M_MG", clean_ident) ~ paste0("scRNA_WT", gsub("^WT_10M_MG", "", clean_ident)),
  grepl("^AD[0-9]$",  clean_ident) ~ paste0("snRNA_AD", gsub("^AD",        "", clean_ident)),
  grepl("^WT[0-9]$",  clean_ident) ~ paste0("snRNA_WT", gsub("^WT",        "", clean_ident)),
  TRUE ~ "Unknown"
)
total$Batch3 <- factor(total$Batch3, levels = c(
  "scRNA_WT1","scRNA_WT2","scRNA_WT3","scRNA_AD1","scRNA_AD2","scRNA_AD3",
  "snRNA_WT1","snRNA_WT2","snRNA_AD1","snRNA_AD2"))
```

### 3.3 Marker-Gene Heatmap and AD vs WT DEGs

```r
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

Idents(total) <- "celltype"

# Per-cluster positive markers
markers  <- FindAllMarkers(total, only.pos = TRUE, logfc.threshold = 0.1)
markers2 <- subset(markers,
                   avg_log2FC >= 0.5 & pct.1 >= 0.5 & pct.2 <= 0.2 &
                   p_val_adj  <= 0.05)

target_genes         <- unique(markers2$gene)
av                   <- AggregateExpression(total, group.by = "celltype")$RNA
expression_data      <- as.data.frame(av)[target_genes, ]
expression_data_clean <- na.omit(expression_data)

# AD vs WT DEGs per cell type
celltypes <- unique(total$celltype)
deg_stats <- data.frame(celltype = character(), direction = character(),
                        DEG_count = integer())
for (ct in celltypes) {
  sub <- subset(total, subset = celltype == ct)
  Idents(sub) <- "Batch1"
  degs <- tryCatch(
    FindMarkers(sub, ident.1 = "AD", ident.2 = "WT",
                min.pct = 0.1, logfc.threshold = 0),
    error = function(e) data.frame())
  if (nrow(degs) > 0) {
    degs <- dplyr::filter(degs, p_val_adj < 0.05)
    up   <- sum(degs$avg_log2FC >=  1)
    down <- sum(degs$avg_log2FC <= -1)
    deg_stats <- rbind(
      deg_stats,
      data.frame(celltype = ct, direction = "Upregulated",   DEG_count =  up),
      data.frame(celltype = ct, direction = "Downregulated", DEG_count = -down))
  }
}
```

---

## 4. Microglial Subclustering (DAM I / DAM II / Homeostatic / IFN-I / Cycling)

### 4.1 Subset and Re-Integrate with scVI

```r
# Load the 13-cell-type annotated object
load(".../total_nodouble.Rdata")

MG <- subset(total, subset = celltype == "Microglia")

# Seurat v5: split layers per sample before re-integration
MG[["RNA"]] <- split(MG[["RNA"]], f = MG$orig.ident)

source(".../seurat_integration_pipeline.R")
norm_methods      <- c("LogNormalize")
integrate_methods <- c("scVI")  # scVI removes batch effects while preserving microglial state heterogeneity

process_seurat_integrate(
  MG,
  norm_methods      = norm_methods,
  integrate_methods = integrate_methods,
  output_dir        = "./Results/integrated_data/",
  verbose           = TRUE
)
# Output: clustered_obj with reduction "integrated.scvi"
```

### 4.2 Five-Subtype Annotation

```r
MG <- subset(MG, scvi_clusters %in% 0:5)

celltype <- data.frame(ClusterID = 0:5, celltype = 0:5)
celltype[celltype$ClusterID %in% c(0),   2] <- "Homeostatic"
celltype[celltype$ClusterID %in% c(1,3), 2] <- "DAM I"
celltype[celltype$ClusterID %in% c(2),   2] <- "DAM II"
celltype[celltype$ClusterID %in% c(4),   2] <- "Cycling"
celltype[celltype$ClusterID %in% c(5),   2] <- "IFN-I responsive"

MG@meta.data$celltype <- "NA"
for (i in seq_len(nrow(celltype))) {
  MG@meta.data[MG$seurat_clusters == celltype$ClusterID[i], "celltype"] <-
    celltype$celltype[i]
}
MG$celltype <- factor(MG$celltype,
  levels = c("Homeostatic","DAM I","DAM II","IFN-I responsive","Cycling"))

DimPlot(MG, group.by = "celltype",
        cols = c("#BFBE91","#E5CC96","#bca4cc","#74ADD1","#D39394","#84a494"),
        pt.size = 0.5) + NoAxes()
```

### 4.3 Microglial Marker Dot Plot

```r
genes <- c(
  "Tmem119","Sall1","Siglech","P2ry12","Hexb","Aif1","Cx3cr1",
  "Lyz2","Spp1","Apoe","Cstb","Fabp5","Prdx1","Fth1",
  "Cst7","Igf1","Etl4","Cox6a2","Wfdc17","Clec7a","Trem2",
  "Ifit3","Isg15","Slfn5","Ifi209","Ifit2","Irf7","Ifit3b",
  "Top2a","Birc5","Mki67","Cenpe","Cdca8","Lockd","Cdca3"
)
DotPlot_scCustom(MG, features = genes,
  colors_use = c("#547298","#8E9DBA","#DADEE7","#F1DDD6","#DA9087","#D24846")) +
  theme(axis.text.x = element_text(face = "italic", angle = 90,
                                   vjust = 0.5, hjust = 1))
```
