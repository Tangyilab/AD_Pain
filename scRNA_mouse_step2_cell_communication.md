# scRNA / snRNA Mouse — Step 2: Cell-Cell Communication

Cell-cell communication analysis between DAM II microglia and the
remaining annotated populations of the mouse spinal cord, using CellChat
v2.2.0 with the mouse `Secreted Signaling` database, plus an orthogonal
LIANA+ analysis on the same annotated object.

This step starts from the integrated, annotated Seurat object produced
by Step 1 (13-cell-type whole-tissue object `total` and the five-subtype
microglial object `MG`).

---

## 1. CellChat — DAM II vs All Other Cell Types

### 1.1 Build the Comparison Object (DAM II + 12 other major populations)

```r
total_subset <- subset(total, subset = celltype %in%
  c("Neuron","Oligo","Astro","OPC","Endo","Ependymal",
    "BAM","Mono","Cycling","T/NK","Peri","Fibro"))

MG_subset <- subset(MG, subset = celltype == "DAM II")

# 13 populations: DAM II replaces the original Microglia label
total <- merge(total_subset, MG_subset)
```

### 1.2 CellChat Wrapper

The wrapper builds and saves an independent CellChat object per group
and, in multi-group mode, also returns the merged object. For each
group it runs:

`createCellChat → subsetData → identifyOverExpressedGenes /
identifyOverExpressedInteractions → computeCommunProb(triMean) →
filterCommunication(min.cells = 10) → computeCommunProbPathway →
aggregateNet`.

```r
scRAN_Cellchat <- function(
  seurat_obj,
  assay          = "RNA",
  celltype.col   = "celltype",
  compare.by     = NULL,
  compare.groups = NULL,
  db.type        = c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor"),
  ncores         = 10,
  min.cells      = 10,
  out.prefix     = "cellchat"
) {
  library(Seurat); library(CellChat); library(future); library(rlang)

  DefaultAssay(seurat_obj) <- assay
  if (assay == "RNA") seurat_obj <- JoinLayers(seurat_obj)
  Idents(seurat_obj) <- celltype.col

  CellChatDB     <- CellChatDB.mouse
  CellChatDB.use <- subsetDB(CellChatDB, search = db.type)

  options(future.globals.maxSize = 300 * 1024^3)
  future::plan("multisession", workers = ncores)

  run_single_cellchat <- function(obj, name = NULL) {
    data_input <- obj[[assay]]$data
    labels     <- Idents(obj)
    meta       <- data.frame(labels = labels, row.names = names(labels))

    cellchat <- createCellChat(object = data_input, meta = meta,
                               group.by = "labels")
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = min.cells)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    if (!is.null(name))
      save(cellchat, file = paste0(out.prefix, "_", name, ".rdata"))
    cellchat
  }

  # Single-sample mode
  if (is.null(compare.by)) {
    cellchat <- run_single_cellchat(seurat_obj)
    save(cellchat, file = paste0(out.prefix, ".rdata"))
    return(cellchat)
  }

  # Multi-sample mode
  if (!compare.by %in% colnames(seurat_obj@meta.data))
    stop("compare.by not found in Seurat meta.data")
  all.groups <- unique(seurat_obj@meta.data[[compare.by]])
  if (is.null(compare.groups)) compare.groups <- all.groups

  cellchat.list <- list()
  for (s in compare.groups) {
    obj.sub <- subset(seurat_obj, subset = !!sym(compare.by) == s)
    cellchat.list[[s]] <- run_single_cellchat(obj.sub, name = s)
  }
  cellchat.merge <- mergeCellChat(cellchat.list,
                                  add.names = names(cellchat.list))
  save(cellchat.merge, file = paste0(out.prefix, ".rdata"))
  cellchat.merge
}
```

### 1.3 Run CellChat Separately for AD and WT, then Merge

```r
setwd("D:/2512_YHC_IGF1/cellchatDAMtoALL")

scRAN_Cellchat(
  total,
  assay          = "RNA",
  celltype.col   = "celltype",
  compare.by     = "Batch1",            # split by AD / WT
  compare.groups = c("AD", "WT"),
  db.type        = c("Secreted Signaling"),
  ncores         = 10,
  min.cells      = 10,
  out.prefix     = "cellchat"
)
```

### 1.4 Merge, Compute Centrality, and Visualize

```r
load("cellchat_AD.rdata"); cellchat_AD <- cellchat
load("cellchat_WT.rdata"); cellchat_WT <- cellchat

object.list <- list(WT = cellchat_WT, AD = cellchat_AD)
for (i in seq_along(object.list))
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# DAM II: signaling-changes scatter (incoming vs outgoing strength)
netAnalysis_signalingChanges_scatter(
  cellchat,
  idents.use = "DAM II",
  color.use  = c("grey10","#547298","#DA9087"),
  dot.size   = 3, label.size = 3,
  xlims = c(0, 1), ylims = c(0, 0.25)
)

# 13-cell-type palette (Microglia recoloured red because the slot is now DAM II)
sc_color <- c("#D24846","#D5B26C","#A499CC","#E07370","#4C8588","#56642a",
              "#A56BA7","#5E4D9A","#B6766C","#8A8CA7","#6D7A57","#c6a5c1","#941456")
ctypes_order <- c("DAM II","Neuron","Oligo","Astro","OPC","Endo",
                  "Ependymal","BAM","Mono","Cycling","T/NK","Peri","Fibro")
sc_color_named <- setNames(sc_color, ctypes_order)

# IGF and SPP1 pathway chord plots, AD vs WT side-by-side
for (pw in c("SPP1", "IGF")) {
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in seq_along(object.list)) {
    netVisual_aggregate(
      object.list[[i]],
      signaling      = pw,
      layout         = "chord",
      color.use      = sc_color_named,
      pt.title       = 20,
      signaling.name = paste(pw, names(object.list)[i]))
  }
}

# Bubble plot of DAM II → target signaling, AD vs WT
netVisual_bubble(
  cellchat,
  sources.use    = 13,
  targets.use    = c(10, 13),
  comparison     = c(1, 2),
  signaling      = "IGF",
  max.dataset    = 2,
  angle.x        = 90,
  remove.isolate = TRUE,
  color.text.use = FALSE
)
```

### 1.5 Per-Cell-Type Igf1r AD-vs-WT Volcano

```r
library(Seurat); library(dplyr); library(ggplot2); library(ggrepel)

seurat_obj       <- total
gene_of_interest <- "Igf1r"
celltype_col     <- "celltype"
group_col        <- "Batch1"
ident_1          <- "AD"
ident_2          <- "WT"

celltypes      <- unique(seurat_obj@meta.data[[celltype_col]])
igf1r_results  <- list()

for (ct in celltypes) {
  cells_subset <- subset(seurat_obj,
                         subset = !!sym(celltype_col) == ct)
  n_ad <- sum(cells_subset@meta.data[[group_col]] == ident_1)
  n_wt <- sum(cells_subset@meta.data[[group_col]] == ident_2)
  if (n_ad < 3 || n_wt < 3) next
  Idents(cells_subset) <- group_col

  markers <- tryCatch(
    FindMarkers(cells_subset,
                ident.1         = ident_1,
                ident.2         = ident_2,
                features        = gene_of_interest,
                logfc.threshold = 0,
                min.pct         = 0,
                test.use        = "wilcox"),
    error = function(e) NULL)

  if (!is.null(markers)) {
    markers$celltype <- ct
    markers$gene     <- rownames(markers)
    markers$n_AD     <- n_ad
    markers$n_WT     <- n_wt
    igf1r_results[[ct]] <- markers
  }
}

igf1r_df <- bind_rows(igf1r_results) %>%
  dplyr::select(celltype, gene, avg_log2FC, pct.1, pct.2,
                p_val, p_val_adj, n_AD, n_WT) %>%
  arrange(p_val_adj)

cell_type_colors <- c(
  "Microglia"="#9F987F","Neuron"="#D5B26C","Oligo"="#A499CC","Astro"="#E07370",
  "OPC"="#4C8588","Endo"="#56642a","Ependymal"="#A56BA7","BAM"="#5E4D9A",
  "Mono"="#B6766C","Cycling"="#8A8CA7","T/NK"="#6D7A57","Peri"="#c6a5c1",
  "Fibro"="#941456")

igf1r_df$neg_log10_padj <- -log10(igf1r_df$p_val_adj)
igf1r_df$neg_log10_padj[is.infinite(igf1r_df$neg_log10_padj)] <- 0

ggplot(igf1r_df, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray40", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray40", linewidth = 0.6) +
  geom_point(aes(color = celltype, size = pct.1 * 100), alpha = 0.8) +
  scale_color_manual(values = cell_type_colors,
                     breaks = names(cell_type_colors)) +
  scale_size_continuous(range = c(2, 10), name = "% Cells\nExpressing\n(AD)",
                        breaks = c(20, 40, 60, 80),
                        labels = c("20%", "40%", "60%", "80%")) +
  geom_text_repel(
    data = dplyr::filter(igf1r_df, p_val_adj < 0.05),
    aes(label = celltype), size = 4, fontface = "bold",
    box.padding = 0.6, point.padding = 0.4,
    segment.color = "gray50", segment.size = 0.4,
    max.overlaps = 20) +
  theme_classic(base_size = 14) +
  labs(title = "Igf1r Differential Expression: AD vs WT",
       x = "log2 Fold Change (AD / WT)",
       y = "-log10(Adjusted P-value)")
```

---

## 2. LIANA+ — Orthogonal Ligand–Receptor Analysis

LIANA+ is run on the same annotated object after exporting it to 10x mtx
format and reloading it with Scanpy. Human ligand–receptor pairs from
the LIANA+ `consensus` resource are translated to mouse orthologs via
HCOP (`min_evidence = 3`).

### 2.1 Export Seurat → 10x mtx and Reload as AnnData

```r
SaveSeuratMatrix <- function(seurat_obj, outdir = ".",
                             filetype = c("10x", "csv")) {
  library(Seurat)
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- JoinLayers(seurat_obj)
  filetype <- match.arg(filetype)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  mat <- Seurat::GetAssayData(seurat_obj, layer = "counts")

  if (filetype == "10x") {
    mtx_dir <- file.path(outdir, "10x_mtx")
    if (!dir.exists(mtx_dir)) dir.create(mtx_dir)

    tmp_mtx <- tempfile(fileext = ".mtx")
    Matrix::writeMM(mat, tmp_mtx)
    R.utils::gzip(tmp_mtx, file.path(mtx_dir, "matrix.mtx.gz"),
                  overwrite = TRUE)

    tmp_feat <- tempfile()
    write.table(
      data.frame(gene_id = rownames(mat),
                 gene_name = rownames(mat),
                 feature_type = "Gene Expression"),
      tmp_feat, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE)
    R.utils::gzip(tmp_feat, file.path(mtx_dir, "features.tsv.gz"),
                  overwrite = TRUE)

    tmp_bar <- tempfile()
    write.table(
      data.frame(barcode = colnames(mat)),
      tmp_bar, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE)
    R.utils::gzip(tmp_bar, file.path(mtx_dir, "barcodes.tsv.gz"),
                  overwrite = TRUE)
  }

  if (filetype == "csv")
    write.csv(as.matrix(mat), file.path(outdir, "matrix.csv"), quote = FALSE)
  invisible(TRUE)
}

# Export the annotated Seurat object plus its meta.data
SaveSeuratMatrix(total, outdir = ".../LIANA", filetype = "10x")
write.csv(total@meta.data, ".../LIANA/10xmtx/meta.csv")
```

### 2.2 Preprocess in Scanpy and Translate Resource to Mouse

```python
import scanpy as sc
import liana as li
import pandas as pd

adata = sc.read_10x_mtx(".../LIANA/10xmtx")
meta  = pd.read_csv(".../LIANA/10xmtx/meta.csv", index_col=0)
adata.obs = meta.copy()

adata.var["mt"]   = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.startswith(("Rpl", "Rps"))
adata.var["hb"]   = adata.var_names.str.contains("^Hb[^(p)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"],
                           inplace=True, log1p=True)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="orig.ident")
adata.write_h5ad("all_celltype.h5ad")

# Human consensus LR resource → mouse orthologs via HCOP (min 3 evidences)
resource = li.rs.select_resource("consensus")
map_df = li.rs.get_hcop_orthologs(
    url="https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/"
        "human_mouse_hcop_fifteen_column.txt.gz",
    columns=["human_symbol", "mouse_symbol"],
    min_evidence=3,
)
map_df = map_df.rename(columns={"human_symbol": "source",
                                "mouse_symbol": "target"})
mouse = li.rs.translate_resource(
    resource, map_df=map_df,
    columns=["ligand", "receptor"], replace=True, one_to_many=3)
```

### 2.3 Per-Sample Rank Aggregate

```python
li.mt.rank_aggregate.by_sample(
    adata         = adata,
    groupby       = "celltype",
    resource_name = "consensus",
    resource      = mouse,
    min_cells     = 5,
    sample_key    = "Batch1",
    expr_prop     = 0.1,
    key_added     = "liana_res",
    use_raw       = False,
    layer         = "counts",
    inplace       = True,
)
adata.write_h5ad("all_celltype_Liana.h5ad")
```

### 2.4 Igf1 → Igf1r Circle Plot, AD vs WT

```python
import scanpy as sc
import liana as li
import pandas as pd

sc_color = [
    "#D24846","#D5B26C","#A499CC","#E07370","#4C8588","#56642a",
    "#A56BA7","#5E4D9A","#B6766C","#8A8CA7","#6D7A57","#c6a5c1","#941456",
]
ctypes_order = [
    "DAM II","Neuron","Oligo","Astro","OPC","Endo",
    "Ependymal","BAM","Mono","Cycling","T/NK","Peri","Fibro",
]

adata    = sc.read_h5ad(".../LIANA/all_celltype_Liana.h5ad")
adata_ad = adata[adata.obs["Batch1"] == "AD"].copy()
adata_wt = adata[adata.obs["Batch1"] == "WT"].copy()

liana_res = adata.uns["liana_res"]
igf1_ad = liana_res[
    (liana_res["Batch1"]          == "AD")    &
    (liana_res["ligand_complex"]  == "Igf1")  &
    (liana_res["receptor_complex"] == "Igf1r")
]
igf1_wt = liana_res[
    (liana_res["Batch1"]          == "WT")    &
    (liana_res["ligand_complex"]  == "Igf1")  &
    (liana_res["receptor_complex"] == "Igf1r")
]

for sub in (adata_ad, adata_wt):
    sub.uns["celltype_colors"] = sc_color
    sub.obs["celltype"] = pd.Categorical(
        sub.obs["celltype"], categories=ctypes_order, ordered=True)

ax_ad = li.pl.circle_plot(
    adata=adata_ad,
    source_labels=["DAM II", "Fibro"],
    liana_res=igf1_ad,
    groupby="celltype",
    score_key="magnitude_rank",
    inverse_score=True,
    pivot_mode="mean",
    figure_size=(6, 6),
)
ax_ad.figure.savefig("IGF1_Igf1r_circle_AD.pdf", bbox_inches="tight")

ax_wt = li.pl.circle_plot(
    adata=adata_wt,
    source_labels=["Fibro"],
    liana_res=igf1_wt,
    groupby="celltype",
    score_key="magnitude_rank",
    inverse_score=True,
    pivot_mode="mean",
    figure_size=(6, 6),
)
ax_wt.figure.savefig("IGF1_Igf1r_circle_WT.pdf", bbox_inches="tight")
```

### 2.5 DAM II → Neuronal Subtypes (Subset Run)

The same LIANA+ pipeline is rerun on a subset that combines DAM II
microglia with neuronal subtypes carrying elevated Igf1r, to resolve the
Igf1 → Igf1r interactions at the subtype level.

```python
sc_color = [
    "#D24846","#9F987F","#A499CC","#E07370","#D5B26C",
    "#2D5C33","#E069A6","#A56BA7","#5E4D9A",
]
ctypes_order = [
    "DAM II",
    "DE_Reln_Gabra2","DE_Rerb1_Zim1","DE_Sox5_Col5a2",
    "DE_Sox5_Col5a2_Enpp1","DE_Sox5_Tac1",
    "MEL","MEP","DI_Pdyn_Gal_Mlxipl",
]

adata = sc.read_h5ad(".../LIANA/DAMtosubset/all_celltype_Liana.h5ad")
adata_ad = adata[adata.obs["Batch1"] == "AD"].copy()
adata_wt = adata[adata.obs["Batch1"] == "WT"].copy()

liana_res = adata.uns["liana_res"]
igf1_ad = liana_res[
    (liana_res["Batch1"]           == "AD")   &
    (liana_res["ligand_complex"]   == "Igf1") &
    (liana_res["receptor_complex"] == "Igf1r")
]
igf1_wt = liana_res[
    (liana_res["Batch1"]           == "WT")   &
    (liana_res["ligand_complex"]   == "Igf1") &
    (liana_res["receptor_complex"] == "Igf1r")
]

for sub in (adata_ad, adata_wt):
    sub.uns["celltype_colors"] = sc_color
    sub.obs["celltype"] = pd.Categorical(
        sub.obs["celltype"], categories=ctypes_order, ordered=True)

ax_ad = li.pl.circle_plot(
    adata=adata_ad,
    source_labels=["DAM II"],
    liana_res=igf1_ad,
    groupby="celltype",
    score_key="magnitude_rank",
    inverse_score=True,
    pivot_mode="mean",
    figure_size=(6, 6),
)
ax_ad.figure.savefig("IGF1_Igf1r_circle_AD.pdf", bbox_inches="tight")

ax_wt = li.pl.circle_plot(
    adata=adata_wt,
    source_labels=["DI_Pdyn_Gal_Mlxipl"],
    liana_res=igf1_wt,
    groupby="celltype",
    score_key="magnitude_rank",
    inverse_score=True,
    pivot_mode="mean",
    figure_size=(6, 6),
)
ax_wt.figure.savefig("IGF1_Igf1r_circle_WT.pdf", bbox_inches="tight")
```
