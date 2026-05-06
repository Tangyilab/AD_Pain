# MiP-seq — Spinal Cord Spatial Transcriptomics Pipeline

Image preprocessing, DAPI-based nuclear segmentation, transcript-to-cell
assignment, expression-matrix construction, and downstream clustering
for the MiP-seq (multiplexed in situ pairwise sequencing) data
(WT_1–3 and 5xFAD_1–3 mouse spinal cord sections).

---

## 0. Per-Sample Folder Layout

Each biological sample lives under
`analysis_result/<group>/<sample>/`. A representative example
(`5xFAD/5xFAD_1/`) is shown below; the same structure applies to
`5xFAD_2`, `5xFAD_3`, `WT_1`, `WT_2`, and `WT_3`.

```
analysis_result/
├── 5xFAD/
│   ├── 5xFAD_1/
│   │   ├── 01_raw_data/                 # multichannel TIFF input (all_tif.tif)
│   │   ├── 02_gene_location/            # transcript-spot table (all_gene_location.csv)
│   │   ├── 03_gene_plot/                # per-gene 2-D scatter (pdf/, png/)
│   │   ├── 04_cell_segment/             # StarDist segmentation outputs
│   │   │   └── STARdist/
│   │   │       ├── figures/
│   │   │       └── out_stardist_qc/
│   │   ├── 05_expression/               # gene × cell counts
│   │   ├── 06_dapi_gene_plot/           # DAPI + gene overlay (pdf/, png/)
│   │   ├── 07_top_20_gene_plot/         # top-20 genes per cluster
│   │   ├── 08_spatial_heatmap_plot/     # spatial expression heatmaps (pdf/, png/)
│   │   ├── stardist/                    # per-sample StarDist run (canonical)
│   │   ├── stardist_0/                  # 0-px dilation variant
│   │   └── stardist_5/                  # 5-px dilation variant
│   ├── 5xFAD_2/        # (same layout)
│   └── 5xFAD_3/        # (same layout)
└── WT/
    ├── WT_1/           # (same layout)
    ├── WT_2/           # (same layout)
    └── WT_3/           # (same layout)
```

`stardist/`, `stardist_0/`, and `stardist_5/` correspond to the same
StarDist segmentation evaluated at three different transcript-assignment
radii (default = 10 px, plus 0 px = strict nuclear assignment, and 5 px
= half radius). Each subfolder holds its own `gene_by_cell_long.csv` and
`cell_centroids.csv`, so the Seurat / Scanpy merge step can switch
between dilation strategies without re-running segmentation.

---

## 1. Inputs

| File | Folder | Description |
|------|--------|-------------|
| `all_tif.tif` | `01_raw_data/` | Multichannel TIFF; channel 0 = DAPI |
| `all_gene_location.csv` | `02_gene_location/` | Long-format transcript table with columns `gene`, `axis-1` (= y / row), `axis-2` (= x / col) |

---

## 2. Per-Sample StarDist Pipeline

The script below is run once per sample folder and produces all
artefacts inside `<sample>/stardist/` (or the dilation-variant
subfolders). It can be executed as:

```bash
python run_stardist_pipeline.py \
  --base_dir /data_result/dengys/spRNA/analysis_result/5xFAD/5xFAD_1 \
  --dapi_channel 0 \
  --prob_thresh 0.4 \
  --nms_thresh 0.2 \
  --radius_px 10
```

```python
# -*- coding: utf-8 -*-
"""
Pipeline:
  base_dir/
    01_raw_data/all_tif.tif
    02_gene_location/all_gene_location.csv

Outputs (under base_dir/stardist/):
  - DAPI.tif
  - nuclei_labels.tif
  - nuclei_boundaries.png
  - stardist_segmentation.png
  - density_smooth.png
  - RNA_density_with_nuclei_boundaries.png
  - QC_RNA_points_on_DAPI.png
  - QC_dilated_range.png
  - gene_by_cell_long.csv          (FINAL 1)
  - cell_centroids.csv             (FINAL 2)
"""

import os, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tifffile import TiffFile, imwrite
from skimage import io
from skimage.segmentation import find_boundaries
from skimage.measure import regionprops_table
from scipy.ndimage import gaussian_filter, distance_transform_edt
from csbdeep.utils import normalize
from stardist.models import StarDist2D


def ensure_dir(p): os.makedirs(p, exist_ok=True)
def savefig(path, dpi=300):
    plt.tight_layout(); plt.savefig(path, dpi=dpi); plt.close()
def minmax01(x):
    x = x.astype(np.float32)
    d = x.max() - x.min()
    return np.zeros_like(x, dtype=np.float32) if d == 0 else (x - x.min()) / d


def extract_dapi_from_multitif(all_tif_path, out_dapi_path, dapi_channel=0):
    f = TiffFile(all_tif_path)
    data = f.series[0].asarray()
    axes = f.series[0].axes or ""
    f.close()
    if "C" in axes:
        dapi = data.take(dapi_channel, axis=axes.index("C"))
    elif data.ndim == 3:
        dapi = data[:, :, dapi_channel] if data.shape[-1] <= 10 \
               else data[dapi_channel, :, :]
    elif data.ndim == 2:
        dapi = data
    else:
        raise ValueError(f"Cannot infer channel axis: ndim={data.ndim}, axes='{axes}'")
    imwrite(out_dapi_path, dapi)
    return dapi, axes, data.shape


def assign_rna_by_dilated_nuclei_nearest(labels, df, H, W, gene_col,
                                         y_col="axis-1", x_col="axis-2",
                                         radius_px=10, flip_y=False):
    """Nearest-nucleus assignment with a distance threshold."""
    df2 = df.copy()
    y = np.round(df2[y_col].to_numpy()).astype(np.int32)
    x = np.round(df2[x_col].to_numpy()).astype(np.int32)
    if flip_y:
        y = (H - 1) - y
    inb = (y >= 0) & (y < H) & (x >= 0) & (x < W)
    df2 = df2.loc[inb].reset_index(drop=True)
    y, x = y[inb], x[inb]
    df2["y"] = y; df2["x"] = x

    nuc_bin = labels > 0
    dist, (iy, ix) = distance_transform_edt(~nuc_bin, return_indices=True)
    nearest_id = labels[iy, ix]
    within     = dist <= float(radius_px)

    cell_id = nearest_id[y, x].astype(np.int32)
    cell_id[~within[y, x]] = 0
    df2["cell_id"] = cell_id
    return df2, dist, nearest_id, within


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base_dir",     required=True)
    ap.add_argument("--dapi_channel", type=int,   default=0)
    ap.add_argument("--prob_thresh",  type=float, default=0.4)
    ap.add_argument("--nms_thresh",   type=float, default=0.2)
    ap.add_argument("--radius_px",    type=int,   default=10)
    ap.add_argument("--flip_y",       action="store_true")
    ap.add_argument("--sigma_rna",    type=float, default=3.0)
    args = ap.parse_args()

    base_dir = args.base_dir
    raw_tif  = os.path.join(base_dir, "01_raw_data",      "all_tif.tif")
    gene_csv = os.path.join(base_dir, "02_gene_location", "all_gene_location.csv")
    outdir   = os.path.join(base_dir, "stardist")
    ensure_dir(outdir)

    # 1) Extract DAPI
    dapi_path = os.path.join(outdir, "DAPI.tif")
    dapi, axes, full_shape = extract_dapi_from_multitif(
        raw_tif, dapi_path, dapi_channel=args.dapi_channel)
    img = io.imread(dapi_path)
    if img.ndim != 2:
        img = np.squeeze(img)
    H, W = img.shape

    # 2) StarDist nuclear segmentation
    img_norm = normalize(img, pmin=1, pmax=99.8, clip=True)
    model = StarDist2D.from_pretrained("2D_versatile_fluo")
    labels, _ = model.predict_instances(
        img_norm,
        prob_thresh=args.prob_thresh,
        nms_thresh =args.nms_thresh)
    num_nuclei = int(labels.max())
    imwrite(os.path.join(outdir, "nuclei_labels.tif"),
            labels.astype(np.uint16))

    # Segmentation QC stats
    total_area = int(np.sum(labels > 0))
    coverage   = total_area / (H * W) * 100.0
    avg_size   = total_area / num_nuclei if num_nuclei > 0 else 0
    print(f"Nuclei={num_nuclei}, coverage={coverage:.2f}%, "
          f"avg_size={avg_size:.1f}px")

    # Segmentation visualization
    fig, ax2 = plt.subplots(1, 2, figsize=(14, 6))
    ax2[0].imshow(img, cmap="gray");      ax2[0].set_title("Raw DAPI")
    ax2[1].imshow(labels, cmap="nipy_spectral")
    ax2[1].set_title(f"StarDist mask (n={num_nuclei})")
    for ax in ax2: ax.axis("off")
    savefig(os.path.join(outdir, "stardist_segmentation.png"))

    boundaries_nuc = find_boundaries(labels, mode="outer")
    dapi01 = minmax01(img)
    vmin, vmax = np.percentile(dapi01, (2, 99))
    plt.figure(figsize=(8, 8))
    plt.imshow(dapi01, cmap="gray", vmin=vmin, vmax=vmax)
    plt.imshow(boundaries_nuc, cmap="Reds", alpha=0.7)
    plt.title("Nuclei boundaries on DAPI"); plt.axis("off")
    savefig(os.path.join(outdir, "nuclei_boundaries.png"))

    # 3) RNA density map (QC)
    df = pd.read_csv(gene_csv)
    gene_col = "gene" if "gene" in df.columns else "Gene"
    density, _, _ = np.histogram2d(
        df["axis-1"].to_numpy(), df["axis-2"].to_numpy(),
        bins=[H, W], range=[[0, H], [0, W]])
    density_smooth = gaussian_filter(density, sigma=float(args.sigma_rna))

    plt.figure(figsize=(12, 6))
    plt.imshow(density_smooth, cmap="hot", origin="lower")
    plt.colorbar(label="RNA density")
    plt.title("Smoothed RNA density (all genes)"); plt.axis("off")
    savefig(os.path.join(outdir, "density_smooth.png"))

    plt.figure(figsize=(12, 6))
    plt.imshow(density_smooth, cmap="hot", origin="lower")
    plt.imshow(boundaries_nuc, cmap="gray", alpha=0.8)
    plt.title("RNA density + nuclei boundaries"); plt.axis("off")
    savefig(os.path.join(outdir, "RNA_density_with_nuclei_boundaries.png"))

    # 4) RNA → cell assignment with dilation radius
    df_assigned, dist, nearest_id, within = assign_rna_by_dilated_nuclei_nearest(
        labels=labels, df=df, H=H, W=W,
        gene_col=gene_col,
        y_col="axis-1", x_col="axis-2",
        radius_px=int(args.radius_px),
        flip_y=bool(args.flip_y))

    assigned_n = int((df_assigned["cell_id"] > 0).sum())
    print(f"Assigned RNAs within {args.radius_px}px: "
          f"{assigned_n}/{len(df_assigned)} "
          f"({assigned_n / len(df_assigned) * 100:.2f}%)")

    # QC: dilation range
    dilated_mask  = (nearest_id * within).astype(np.int32)
    boundaries_dil = find_boundaries(dilated_mask > 0, mode="outer")
    plt.figure(figsize=(8, 8))
    plt.imshow(dapi01, cmap="gray", vmin=vmin, vmax=vmax)
    plt.imshow(boundaries_nuc, cmap="cool", alpha=0.8)
    plt.imshow(boundaries_dil, cmap="Reds", alpha=0.6)
    plt.title(f"Nuclei (cool) vs dilated assign area (red), r={args.radius_px}px")
    plt.axis("off")
    savefig(os.path.join(outdir, "QC_dilated_range.png"))

    # QC: random RNA points on DAPI
    sample_n = min(800, len(df_assigned))
    df_s = df_assigned.sample(sample_n, random_state=0)
    plt.figure(figsize=(8, 8))
    plt.imshow(dapi01, cmap="gray", vmin=vmin, vmax=vmax)
    in_mask = df_s["cell_id"] > 0
    plt.scatter(df_s.loc[ in_mask, "x"], df_s.loc[ in_mask, "y"],
                s=6, c="lime", alpha=0.6, label="assigned")
    plt.scatter(df_s.loc[~in_mask, "x"], df_s.loc[~in_mask, "y"],
                s=6, c="red",  alpha=0.35, label="unassigned")
    plt.legend(loc="lower right")
    plt.title("RNA points on DAPI (assigned vs unassigned)")
    plt.axis("off")
    savefig(os.path.join(outdir, "QC_RNA_points_on_DAPI.png"))

    # 5) FINAL output 1: gene × cell long table
    gene_by_cell = (
        df_assigned.loc[df_assigned["cell_id"] > 0]
                   .groupby(["cell_id", gene_col])
                   .size().rename("count").reset_index()
                   .sort_values(["cell_id", "count"],
                                ascending=[True, False]))
    gene_by_cell.to_csv(os.path.join(outdir, "gene_by_cell_long.csv"),
                        index=False)

    # 6) FINAL output 2: nuclear centroids
    props = regionprops_table(labels, properties=("label", "centroid"))
    cent  = pd.DataFrame(props).rename(
        columns={"label": "cell_id",
                 "centroid-0": "y", "centroid-1": "x"})
    valid_cells = set(gene_by_cell["cell_id"].astype(int).unique().tolist())
    cent = cent[cent["cell_id"].isin(valid_cells)].copy()
    cent.to_csv(os.path.join(outdir, "cell_centroids.csv"), index=False)


if __name__ == "__main__":
    main()
```

---

## 3. Per-Sample Visualization Folders

These folders are populated from the StarDist outputs and the
`all_gene_location.csv` table; representative panels are listed below.

- `03_gene_plot/`: one PDF/PNG per gene plus a 32-gene-panel scatter
  (Forbidden-City palette, Apoe down-weighted to α = 0.25 to avoid
  overwhelming the panel).
- `04_cell_segment/STARdist/figures/`,
  `04_cell_segment/STARdist/out_stardist_qc/`:
  segmentation overlays (`stardist_segmentation.png`,
  `nuclei_boundaries.png`, `QC_dilated_range.png`,
  `QC_RNA_points_on_DAPI.png`).
- `05_expression/`: per-cell count matrices (`gene_by_cell_long.csv`,
  optional wide-format pivot).
- `06_dapi_gene_plot/`: DAPI + per-gene scatter overlays.
- `07_top_20_gene_plot/`: top-20 genes per Leiden cluster.
- `08_spatial_heatmap_plot/`: spatial expression heatmaps.

### 32-Gene Spatial Scatter (`03_gene_plot/`)

```python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

H, W = 3199, 5258
df   = pd.read_csv("all_gene_location.csv")
genes = list(df["gene"].unique())
assert len(genes) == 32, f"expected 32 genes, got {len(genes)}"

palette32 = [
    "#d60000","#8c3bff","#018700","#00acc6","#97ff00","#ff7ed1","#6b004f",
    "#ffa52f","#00009c","#857067","#004942","#4f2a00","#00fdcf","#bcb6ff",
    "#95b379","#bf03b8","#2466a1","#280041","#dbb3af","#fdf490","#4f445b",
    "#a37c00","#ff7066","#3f806e","#82000c","#a37bb3","#344d00","#9ae4ff",
    "#eb0077","#2d000a","#5d90ff","#00c61f",
]
gene2color = {g: palette32[i] for i, g in enumerate(genes)}
df_plot          = df.copy()
df_plot["color"] = df_plot["gene"].map(gene2color)
df_apoe  = df_plot[df_plot["gene"] == "Apoe"]
df_other = df_plot[df_plot["gene"] != "Apoe"]

fig, ax = plt.subplots(figsize=(14, 10))
fig.patch.set_facecolor("black"); ax.set_facecolor("black")

ax.scatter(df_other["axis-2"], df_other["axis-1"],
           s=1, c=df_other["color"].values, linewidths=0, alpha=0.9)
ax.scatter(df_apoe["axis-2"],  df_apoe["axis-1"],
           s=1, c=df_apoe["color"].values,  linewidths=0, alpha=0.25)

ax.invert_yaxis(); ax.set_xlim(0, W); ax.set_ylim(H, 0)
ax.set_xlabel("axis-2 (x)", color="white")
ax.set_ylabel("axis-1 (y)", color="white")
ax.set_title("32-gene panel", color="white")
ax.tick_params(colors="white")
for spine in ax.spines.values():
    spine.set_color("white")

handles = [Line2D([0], [0], marker="o", color="none",
                  markerfacecolor=gene2color[g], markersize=6, label=g)
           for g in genes]
leg = ax.legend(handles=handles, title="Genes", loc="upper left",
                bbox_to_anchor=(1.02, 1), frameon=False)
plt.setp(leg.get_texts(), color="white")
plt.setp(leg.get_title(), color="white")
plt.tight_layout()
plt.savefig("32_gene_panel_black.png", dpi=1000,
            facecolor=fig.get_facecolor(), bbox_inches="tight")
```

---

## 4. Multi-Sample Aggregation and Clustering (Seurat)

A Seurat-based alternative loads each sample's `gene_by_cell_long.csv`
into a sparse matrix, builds one Seurat object per sample, merges the
six samples, and refines clusters across iterative noise-removal
rounds. Two dilation variants (`stardist_0` = 0-px / strict nuclear,
`stardist`/`stardist_5` = 5-px) are processed independently.

```r
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

# 1) Long CSV → gene × cell sparse matrix
long_to_sparse <- function(csv_path, cell_prefix) {
  dt <- fread(csv_path)
  stopifnot(all(c("cell_id", "gene", "count") %in% names(dt)))
  dt[, cell_id := as.character(cell_id)]
  dt[, gene    := as.character(gene)]
  dt[, count   := as.numeric(count)]
  dt <- dt[!is.na(count) & count > 0]
  dt <- dt[, .(count = sum(count)), by = .(cell_id, gene)]
  dt[, cell_id := paste0(cell_prefix, "_", cell_id)]

  gene_levels <- sort(unique(dt$gene))
  cell_levels <- sort(unique(dt$cell_id))
  dt[, i := match(gene,    gene_levels)]
  dt[, j := match(cell_id, cell_levels)]

  sparseMatrix(
    i        = dt$i,
    j        = dt$j,
    x        = dt$count,
    dims     = c(length(gene_levels), length(cell_levels)),
    dimnames = list(gene_levels, cell_levels))
}

# 2) Build one Seurat object per sample
read_one_sample <- function(sample_dir, sample_name) {
  csv_path <- file.path(sample_dir, "gene_by_cell_long.csv")
  if (!file.exists(csv_path)) stop("Missing: ", csv_path)

  counts <- long_to_sparse(csv_path, cell_prefix = sample_name)
  seu <- CreateSeuratObject(counts = counts, project = sample_name,
                            min.cells = 0, min.features = 0)
  seu$genotype <- sub("^(WT|5xFAD).*", "\\1", sample_name)
  seu$sample   <- sample_name
  seu
}

# 3) Sample list (canonical 10-px dilation under "stardist/")
root    <- "/data_result/dengys/spRNA/analysis_result"
samples <- c(
  "WT/WT_1/stardist",     "WT/WT_2/stardist",     "WT/WT_3/stardist",
  "5xFAD/5xFAD_1/stardist","5xFAD/5xFAD_2/stardist","5xFAD/5xFAD_3/stardist"
)
sample_dirs  <- file.path(root, samples)
sample_names <- sub(".*/(WT_[0-9]+|5xFAD_[0-9]+)(/.*)?$", "\\1", samples)

seu_list <- lapply(seq_along(sample_dirs), function(k)
  read_one_sample(sample_dir = sample_dirs[k],
                  sample_name = sample_names[k]))
names(seu_list) <- sample_names

seu_all <- merge(seu_list[[1]], y = seu_list[-1],
                 add.cell.ids = names(seu_list))
seu_all$sample   <- sub("^((WT_[0-9]+|5xFAD_[0-9]+)).*$", "\\1",
                        colnames(seu_all))
seu_all$genotype <- sub("^((WT|5xFAD)).*$", "\\1", seu_all$sample)

saveRDS(seu_all, file = "WT_5xFAD_merged_seurat.rds")
```

### 4.1 Iterative Clustering with Noise Removal

Cluster-level noise (high-dropout / low-feature outlier clusters) is
removed across multiple rounds, with `FindClusters` re-run after each
filtering step. Below is the strict-nuclear (0-px) variant; the 5-px
variant follows the same template at lower resolution.

```r
total <- readRDS("WT_5xFAD_merged_seurat.rds")
pbmc  <- JoinLayers(total)

run_round <- function(obj, dims = 1:10, resolution = 2) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, features = rownames(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  obj <- FindNeighbors(obj, dims = dims)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = dims)
  obj
}

# Round 0: initial clustering
pbmc <- run_round(pbmc, dims = 1:10, resolution = 2)

# Round 1: drop noise clusters identified from VlnPlot of nCount/nFeature
pbmc <- subset(pbmc, subset = seurat_clusters %in% c(0:18, 20:30, 44))
pbmc <- run_round(pbmc, dims = 1:10, resolution = 2)

# Round 2
pbmc <- subset(pbmc, subset = seurat_clusters %in%
                 c(0:39, 42:48, 44, 53, 54, 58))
pbmc <- run_round(pbmc, dims = 1:10, resolution = 2)

# Round 3 (final at resolution = 1)
pbmc <- subset(pbmc, subset = seurat_clusters %in% c(0:47, 49, 50))
pbmc <- run_round(pbmc, dims = 1:10, resolution = 1)

save(pbmc, file = "umap_expansion0.Rdata")
```

For the 5-px dilation variant the sample list points to `stardist_5/`
(or the canonical `stardist/`) instead of `stardist_0/`, the same
`run_round` helper is used, and the final Leiden / Louvain resolution
is reduced to 0.5 since transcripts associated with cytoplasmic shells
yield fewer but more distinct clusters.

---

## 5. Downstream Analysis

All sections below operate on the merged, iteratively-cleaned `pbmc`
Seurat object from Section 5, after cell-type annotation
(`cellname` = fine-grained subtype; `celltype` = collapsed class).

### 5.1 Cellname / Gene-Name Cleanup

A handful of gene names and cellname labels were typoed during the
panel design and clustering rounds; they are corrected once before
downstream analysis. Cellname is then collapsed into a 4-class
`celltype` field used for cross-class comparisons.

```r
library(dplyr)

# Gene-name typos in the rownames of the count matrix
rownames(pbmc)[rownames(pbmc) == "Igfl"]   <- "Igf1"
rownames(pbmc)[rownames(pbmc) == "Igflr"]  <- "Igf1r"
rownames(pbmc)[rownames(pbmc) == "Pydn"]   <- "Pdyn"
rownames(pbmc)[rownames(pbmc) == "Rerb1"]  <- "Rreb1"
rownames(pbmc)[rownames(pbmc) == "Vgat"]   <- "Slc32a1"
rownames(pbmc)[rownames(pbmc) == "Vglut2"] <- "Slc17a6"

# Cellname label corrections
pbmc$cellname <- gsub("^Exh_",       "Glut_",    pbmc$cellname)
pbmc$cellname <- gsub("Rerb1",       "Rreb1",    pbmc$cellname)
pbmc$cellname <- gsub("^Nos1_high$", "Inh_Nos1", pbmc$cellname)

# Cellname → 4-class celltype
pbmc$celltype <- case_when(
  grepl("^Inh_",  as.character(pbmc$cellname)) ~ "GABA",
  grepl("^Glut_", as.character(pbmc$cellname)) ~ "Glut",
  as.character(pbmc$cellname) == "MG"          ~ "microglia",
  as.character(pbmc$cellname) == "MN"          ~ "ACh",
  as.character(pbmc$cellname) == "MEP"         ~ "Glut",
  TRUE                                         ~ "Other")

# Fixed 23-color palette tied to cellname
my_colors <- c(
  "Glut_Cck"        = "#EFD6BB", "Glut_Cck_Tacr1"   = "#368550",
  "Glut_Crh"        = "#9F977E", "Glut_Grpr"        = "#1F567A",
  "Glut_Maf"        = "#B6B812", "Glut_Reln"        = "#D4DA47",
  "Glut_Rorb"       = "#A46BA7", "Glut_Rorb_Rreb1"  = "#855F16",
  "Glut_Rreb1"      = "#78C2EC", "Glut_Rreb1_Reln"  = "#DFA7C7",
  "Glut_Sst"        = "#9DC3C3", "Glut_Sst_Tac1"    = "#279AD6",
  "Glut_Tac2"       = "#931355", "Inh_Gal_Pdyn"     = "#75C7CB",
  "Inh_Nfib"        = "#D4B16C", "Inh_Nos1"         = "#5E4D9A",
  "Inh_Npy_Pdyn"    = "#7BBA5F", "Inh_Pax2_Pvalb"   = "#FBBB10",
  "Inh_Pvalb"       = "#DFDFEC", "Inh_Rreb1"        = "#019FA7",
  "MEP"             = "#DF69A6", "MG"               = "#EE7A77",
  "MN"              = "#A398CB")
```

### 5.2 Per-Sample Spatial Scatter Helper

`generate_plot_df()` reads the per-sample tissue-outline polygon
(`<sample>_XY_Coordinates.csv`) plus the StarDist centroid table
(`<sample>_cell_centroids.csv`), filters cells inside the polygon,
rotates them by a per-sample angle so all sections share a common
dorsal-up orientation, and joins back to the Seurat metadata.

```r
generate_plot_df <- function(SeuratOBJ, sample_name, rotation_angle,
                             custom_colors = NULL) {
  library(sf); library(dplyr)

  fanwei <- read.csv(paste0(sample_name, "_XY_Coordinates.csv"))
  label  <- read.csv(paste0(sample_name, "_cell_centroids.csv"))

  if (!(fanwei$X[1] == fanwei$X[nrow(fanwei)] &&
        fanwei$Y[1] == fanwei$Y[nrow(fanwei)])) {
    fanwei <- rbind(fanwei, fanwei[1, ])
  }
  poly_sf <- st_sfc(st_polygon(list(as.matrix(fanwei[, c("X", "Y")]))))
  pts     <- st_as_sf(label, coords = c("x", "y"), crs = NA)
  inside  <- st_within(pts, poly_sf, sparse = FALSE)[, 1]
  label_in <- label[inside, ]

  rotate_points <- function(df, angle_deg, x_col = "x", y_col = "y") {
    theta <- angle_deg * pi / 180
    cx <- mean(df[[x_col]]); cy <- mean(df[[y_col]])
    df$x_rot <- (df[[x_col]] - cx) * cos(theta) -
                (df[[y_col]] - cy) * sin(theta) + cx
    df$y_rot <- (df[[x_col]] - cx) * sin(theta) +
                (df[[y_col]] - cy) * cos(theta) + cy
    df
  }
  label_in <- rotate_points(label_in, rotation_angle)
  label_in <- label_in %>%
    mutate(barcode = paste0(sample_name, "_", sample_name, "_", cell_id))

  meta <- SeuratOBJ@meta.data %>%
    mutate(barcode = rownames(SeuratOBJ@meta.data)) %>%
    filter(sample == sample_name) %>%
    select(barcode, cellname, seurat_clusters)

  label_in %>% left_join(meta, by = "barcode")
}

# Per-sample rotation angles
samples <- list(
  list(id = "WT_1",    angle =  147),
  list(id = "WT_2",    angle =   -7),
  list(id = "WT_3",    angle =  -97),
  list(id = "5xFAD_1", angle =    7),
  list(id = "5xFAD_2", angle =   10),
  list(id = "5xFAD_3", angle =  -86))

# Plot one section
output <- generate_plot_df(pbmc, "WT_3", -97)
ggplot(output, aes(x_rot, y_rot, color = cellname)) +
  geom_point(size = 1.5, stroke = 0) +
  scale_color_manual(values = my_colors, na.value = "grey80") +
  scale_y_reverse() + coord_fixed() + theme_classic() +
  theme(legend.title = element_blank())
```

### 5.3 Step 1 — Marker-Gene Stacked Violin and UMAP

```r
genes <- c("Cck","Tacr1","Crh","Grpr","Maf","Reln",
           "Rorb","Rreb1","Sst","Tac1","Tac2","Gal",
           "Pdyn","Nfib","Nos1","Npy","Pax2","Pvalb",
           "Zfhx3","Chat","Tmem119","Hexb","Cst7","Apoe")

Stacked_VlnPlot(pbmc, features = genes, group.by = "cellname",
                pt.size = 0, colors_use = my_colors, x_lab_rotate = 90)

plot_scatter_scRNA(
  pbmc, group.by = "cellname", cols = my_colors,
  save = TRUE, filename = "UMAP_celltype.pdf",
  reductionname = "umap", legend.position = "right",
  legend.nrow = 14, width = 8, height = 5)

# 4-class collapsed UMAP
celltype_colors <- c(
  Glut       = "#1F577B", GABA      = "#A56BA7",
  microglia  = "#E0A7C8", ACh       = "#7CBB5F",
  Other      = "#368650")
plot_scatter_scRNA(
  pbmc, group.by = "celltype",
  cols = celltype_colors[unique(pbmc$celltype)],
  save = TRUE, filename = "UMAP_celltype_new.pdf",
  reductionname = "umap", width = 6, height = 5)
```

### 5.4 Step 2 — Cross-Modality Correlation with snRNA-seq

Per-cellname / per-Family mean expression on the shared 32-gene panel
is gene-wise z-scored, then Spearman-correlated to verify that MiP-seq
clusters recover the same neuronal families called from snRNA-seq.

```r
library(Seurat)
pbmc       <- NormalizeData(pbmc,       verbose = FALSE)
Neu_rename <- NormalizeData(Neu_rename, verbose = FALSE)

pbmc_avg <- AverageExpression(pbmc, group.by = "cellname",
                              assays = "RNA", features = genes32,
                              slot = "data")$RNA
fish_avg <- AverageExpression(Neu_rename, group.by = "Neuron_name",
                              assays = "RNA", features = genes32,
                              slot = "data")$RNA
colnames(pbmc_avg) <- gsub("-", "_", colnames(pbmc_avg))
colnames(fish_avg) <- gsub("-", "_", colnames(fish_avg))

common <- intersect(rownames(pbmc_avg), rownames(fish_avg))
pb <- as.matrix(pbmc_avg[common, , drop = FALSE])
fh <- as.matrix(fish_avg[common, , drop = FALSE])

keep <- is.finite(apply(pb, 1, sd)) & is.finite(apply(fh, 1, sd)) &
        apply(pb, 1, sd) > 0 & apply(fh, 1, sd) > 0
pb_z <- t(scale(t(pb[keep, ])))
fh_z <- t(scale(t(fh[keep, ])))

cor_mat <- cor(fh_z, pb_z, method = "spearman",
               use = "pairwise.complete.obs")
pheatmap::pheatmap(cor_mat, border_color = NA)
```

### 5.5 Step 2.1 — Within-MiP-seq Cellname Self-Correlation

```r
library(ggplot2); library(dplyr); library(grid)

mat   <- as.matrix(pbmc_avg)
mat_z <- t(scale(t(mat)))
mat_z <- mat_z[apply(mat_z, 1, function(x) all(is.finite(x))), ,
               drop = FALSE]

cor_mat <- cor(mat_z, method = "spearman",
               use = "pairwise.complete.obs")
cor01   <- (cor_mat + 1) / 2
ord     <- colnames(cor01)
cor01   <- cor01[ord, ord]

df <- as.data.frame(as.table(cor01))
colnames(df) <- c("y", "x", "val")
df$row_i <- match(df$y, ord); df$col_i <- match(df$x, ord)
df <- df %>% filter(row_i >= col_i)
df$x <- factor(df$x, levels = ord)
df$y <- factor(df$y, levels = rev(ord))

ggplot(df, aes(x = x, y = y, fill = val)) +
  geom_tile() + coord_fixed() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(position = "right", expand = c(0, 0)) +
  scale_fill_gradient2(midpoint = 0.5, limits = c(0, 1),
                       breaks = c(0, 0.5, 1), name = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        legend.key.height = unit(2.2, "cm"))
```

### 5.6 Step 3 — Per-Gene Positive-Cell Distribution Across Neuronal Classes

For each panel gene, count cells with `counts > 0`, split by collapsed
neuronal class (Excitatory / Inhibitory / Cholinergic), excluding
microglia.

```r
library(Matrix); library(dplyr); library(tidyr)
library(ggplot2); library(scales)

genes    <- rownames(pbmc)
mat      <- GetAssayData(pbmc, assay = "RNA",
                         slot  = "counts")[genes, , drop = FALSE]
expr_bin <- mat > 0

idx_list  <- split(seq_len(ncol(expr_bin)), pbmc$celltype)
count_mat <- sapply(idx_list, function(idx)
                    Matrix::rowSums(expr_bin[, idx, drop = FALSE]))

df2 <- as.data.frame(count_mat)
df2$gene <- rownames(df2)
df2_long <- df2 %>%
  pivot_longer(-gene, names_to = "celltype", values_to = "n") %>%
  mutate(celltype2 = recode(celltype,
                            Glut = "Excitatory",
                            GABA = "Inhibitory",
                            ACh  = "Cholinergic")) %>%
  filter(celltype != "microglia") %>%
  group_by(gene) %>% mutate(total = sum(n)) %>% ungroup() %>%
  arrange(total) %>%
  mutate(gene = factor(gene, levels = unique(gene)),
         celltype2 = factor(celltype2,
           levels = c("Excitatory","Inhibitory","Cholinergic")))

ggplot(df2_long, aes(x = gene, y = n, fill = celltype2)) +
  geom_col(width = 0.75) +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(
    Excitatory  = "#E78B8B",
    Inhibitory  = "#6F7DC5",
    Cholinergic = "#83C29A")) +
  labs(x = NULL, y = "Cell count", fill = NULL) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
```

### 5.7 Step 4 — Igf1 / Igf1r Cell-Type Composition (WT vs 5xFAD)

Per-cell binarized expression of Igf1 and Igf1r is broken down into
cellname pies for WT vs 5xFAD, with the title showing the percentage
of positive cells out of the genotype's total.

```r
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)

DefaultAssay(pbmc) <- "RNA"
genes <- c("Igf1", "Igf1r")

mat <- GetAssayData(pbmc, assay = "RNA",
                    slot = "counts")[genes, , drop = FALSE]
det  <- t(as.matrix(mat) > 0)
meta <- pbmc@meta.data
meta$Igf1_pos  <- det[rownames(meta), "Igf1"]
meta$Igf1r_pos <- det[rownames(meta), "Igf1r"]

geno_label <- c("WT" = "WT", "5xFAD" = "AD")
cellname_order <- names(my_colors)

make_pie <- function(meta, geno, gene,
                     color_pal     = my_colors,
                     name_order    = cellname_order,
                     label_min_pct = 4) {
  pos_col <- paste0(gene, "_pos")
  sub <- meta %>%
    filter(genotype == geno, .data[[pos_col]] == TRUE,
           !is.na(cellname))
  n_pos <- nrow(sub)
  n_all <- sum(meta$genotype == geno)
  pct_pos_all <- ifelse(n_all > 0, 100 * n_pos / n_all, NA_real_)

  if (n_pos == 0)
    return(ggplot() + theme_void() +
             ggtitle(sprintf("%s | %s+  (0 cells)",
                             geno_label[geno], gene)))

  df <- sub %>% count(cellname, name = "n") %>%
    mutate(cellname = factor(cellname, levels = name_order),
           frac = n / sum(n), pct = 100 * frac) %>%
    arrange(cellname) %>%
    mutate(ypos = cumsum(frac) - 0.5 * frac,
           lab  = ifelse(pct >= label_min_pct,
                         paste0(round(pct), "%"), ""))

  cols_use <- color_pal[levels(df$cellname)[
    levels(df$cellname) %in% df$cellname]]

  ggplot(df, aes(x = 1, y = frac, fill = cellname)) +
    geom_col(width = 1, color = "white", linewidth = 0.35) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols_use, drop = TRUE,
                      name = "Cell type",
                      guide = guide_legend(ncol = 1)) +
    theme_void(base_size = 11) +
    ggtitle(sprintf("%s | %s+  (%.1f%% of %s)",
                    geno_label[geno], gene, pct_pos_all,
                    geno_label[geno]))
}

final <- (make_pie(meta, "WT", "Igf1")  | make_pie(meta, "5xFAD", "Igf1")) /
         (make_pie(meta, "WT", "Igf1r") | make_pie(meta, "5xFAD", "Igf1r")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
```

### 5.8 Step 5 — Per-Celltype Spatial Distribution (WT_3 vs 5xFAD_1)

For the two representative sections (WT_3 and 5xFAD_1), each rotated
spatial scatter is split into ACh / GABA / Glut / microglia panels and
colored by cellname.

```r
library(dplyr); library(ggplot2)

cellname2type <- pbmc@meta.data %>%
  select(cellname, celltype) %>%
  filter(!is.na(cellname), !is.na(celltype)) %>% distinct()

add_celltype <- function(df_plot, cellname2type)
  df_plot %>% left_join(cellname2type, by = "cellname")

plot_one_celltype <- function(df_plot, sample_id, target_type,
                              my_colors, pt_size = 0.35, alpha = 0.9) {
  df_sub <- df_plot %>%
    filter(celltype == target_type, !is.na(cellname)) %>%
    mutate(cellname = as.character(cellname))
  cols_use <- my_colors[names(my_colors) %in% unique(df_sub$cellname)]

  ggplot(df_sub, aes(x = x_rot, y = y_rot, color = cellname)) +
    geom_point(size = pt_size, alpha = alpha) +
    scale_color_manual(values = cols_use, drop = FALSE) +
    coord_equal() + theme_void() +
    ggtitle(paste0(sample_id, " | ", target_type,
                   " (colored by cellname)"))
}

plot_4types_and_save <- function(df_plot, sample_id, my_colors,
                                 outdir = ".") {
  for (tp in c("ACh", "GABA", "Glut", "microglia")) {
    p <- plot_one_celltype(df_plot, sample_id, tp, my_colors)
    ggsave(file.path(outdir,
                     paste0("Spatial_", sample_id, "_", tp, ".pdf")),
           p, width = 7.2, height = 6.2, useDingbats = FALSE)
    ggsave(file.path(outdir,
                     paste0("Spatial_", sample_id, "_", tp, ".png")),
           p, width = 7.2, height = 6.2, dpi = 300)
  }
}

df_WT3    <- generate_plot_df(pbmc, "WT_3",    -97)
df_5xFAD1 <- generate_plot_df(pbmc, "5xFAD_1",   7)
plot_4types_and_save(add_celltype(df_WT3,    cellname2type),
                     "WT_3",    my_colors)
plot_4types_and_save(add_celltype(df_5xFAD1, cellname2type),
                     "5xFAD_1", my_colors)
```

### 5.9 Step 6 — Igf1 / Igf1r Spatial Expression

Inject the rotated centroid coordinates back into a Seurat
`DimReducObject` so that `FeaturePlot_scCustom` can render Igf1 and
Igf1r in the original tissue geometry. Y is flipped so dorsal is up.

```r
library(Seurat); library(ggplot2)

attach_spatial <- function(pbmc, output, sample_name) {
  obj <- subset(pbmc, subset = sample == sample_name)
  out_s <- output[grepl(paste0("^", sample_name, "_"), output$barcode), ]
  coords <- out_s[match(rownames(obj@meta.data), out_s$barcode),
                  c("x_rot", "y_rot")]
  rownames(coords) <- rownames(obj@meta.data)
  colnames(coords) <- c("spatial_1", "spatial_2")
  coords$spatial_2 <- -coords$spatial_2

  keep   <- complete.cases(coords)
  obj    <- obj[, keep]; coords <- coords[keep, ]

  obj[["spatial"]] <- CreateDimReducObject(
    embeddings = as.matrix(coords),
    key        = "spatial_",
    assay      = DefaultAssay(obj))
  obj
}

output_wt3 <- generate_plot_df(pbmc, "WT_3",    -97)
output_ad1 <- generate_plot_df(pbmc, "5xFAD_1",   7)
pbmc_wt3   <- attach_spatial(pbmc, output_wt3, "WT_3")
pbmc_ad1   <- attach_spatial(pbmc, output_ad1, "5xFAD_1")

heatmap_colors <- colorRampPalette(c(
  "#547298","#8E9DBA","#DADEE7","#F1DDD6","#DA9087","#D24846"))(100)

for (g in c("Igf1", "Igf1r")) {
  pdf(paste0(g, "_spatial_WT_3.pdf"),    width = 5, height = 4)
    print(FeaturePlot_scCustom(pbmc_wt3, g, reduction = "spatial",
                               colors_use = heatmap_colors,
                               pt.size = 1) + NoAxes())
  dev.off()
  pdf(paste0(g, "_spatial_5xFAD_1.pdf"), width = 5, height = 4)
    print(FeaturePlot_scCustom(pbmc_ad1, g, reduction = "spatial",
                               colors_use = heatmap_colors,
                               pt.size = 1) + NoAxes())
  dev.off()
}

# MG Igf1+ vs Igf1- spatial highlight
pain_related <- c("Glut_Cck_Tacr1","Glut_Reln","Glut_Sst",
                  "Glut_Sst_Tac1","MEP")

flag_mg_igf1 <- function(obj) {
  d <- FetchData(obj, vars = c("Igf1", "cellname"))
  d$color_group <- "Other"
  d$color_group[d$cellname == "MG" & d$Igf1 == 0] <- "MG_Igf1_neg"
  d$color_group[d$cellname == "MG" & d$Igf1  > 0] <- "MG_Igf1_pos"
  obj$color_group <- d$color_group
  obj
}

pbmc_ad1 <- flag_mg_igf1(pbmc_ad1)
pbmc_wt3 <- flag_mg_igf1(pbmc_wt3)

for (lab in list(c("ad1", pbmc_ad1), c("wt3", pbmc_wt3))) {
  pdf(paste0("MG_Igf1_plot_", lab[[1]], ".pdf"), width = 6, height = 4)
    print(DimPlot(lab[[2]], group.by = "color_group",
                  reduction = "spatial", pt.size = 1,
                  cols = c(MG_Igf1_pos = "#DA9087",
                           MG_Igf1_neg = "#547298",
                           Other       = "grey85"),
                  order = c("MG_Igf1_pos", "MG_Igf1_neg")) +
            NoAxes() + ggtitle("MG cells: Igf1+ vs Igf1-"))
  dev.off()
}

# Pain-related Glut Igf1r+ vs Igf1r- spatial highlight (same template)
flag_neuron_igf1r <- function(obj, pain_related) {
  d <- FetchData(obj, vars = c("Igf1r", "cellname"))
  d$color_group <- "Other"
  d$color_group[d$cellname %in% pain_related & d$Igf1r == 0] <-
    "Glut_pain_related_Igf1r-"
  d$color_group[d$cellname %in% pain_related & d$Igf1r  > 0] <-
    "Glut_pain_related_Igf1r+"
  obj$color_group <- d$color_group
  obj
}
```

### 5.10 Step 7 — MG–Neuron Spatial Proximity (Multi-Radius)

For each MG cell, count Igf1r+ pain-related neurons within radii of
25 / 50 / 100 / 150 μm (pixel scale = 0.52 μm/px); compare WT vs 5xFAD
both cell-level (violin) and sample-level (3-dot summary).

```r
library(dplyr); library(ggplot2); library(rstatix)

px_per_um <- 1 / 0.52
radii <- c("25um"  =  25 * px_per_um,
           "50um"  =  50 * px_per_um,
           "100um" = 100 * px_per_um,
           "150um" = 150 * px_per_um)

# Build a single all-sample data frame with rotated coords + role flags
all_output <- bind_rows(lapply(samples, function(s) {
  df <- generate_plot_df(pbmc, s$id, s$angle)
  df$sample <- s$id
  df
}))

expr <- FetchData(pbmc, vars = c("Igf1", "Igf1r", "cellname"))
expr$barcode <- rownames(expr)
all_output <- left_join(all_output, expr, by = "barcode") %>%
  mutate(cell_role = case_when(
    cellname.x == "MG" & Igf1  > 0 ~ "MG_Igf1pos",
    cellname.x == "MG" & Igf1 == 0 ~ "MG_Igf1neg",
    cellname.x %in% pain_related & Igf1r  > 0 ~ "Neuron_Igf1r_pos",
    cellname.x %in% pain_related & Igf1r == 0 ~ "Neuron_Igf1r_neg",
    TRUE ~ "Other"))

analyze_neighbors <- function(df, radius = 150) {
  mg_pos     <- df %>% filter(cell_role == "MG_Igf1pos")
  mg_neg     <- df %>% filter(cell_role == "MG_Igf1neg")
  neuron_pos <- df %>% filter(cell_role == "Neuron_Igf1r_pos")
  if (nrow(mg_pos) == 0 || nrow(neuron_pos) == 0) return(NULL)

  calc_stats <- function(mg_df, neuron_df, mg_label) {
    bind_rows(lapply(seq_len(nrow(mg_df)), function(i) {
      d <- sqrt((neuron_df$x_rot - mg_df$x_rot[i])^2 +
                (neuron_df$y_rot - mg_df$y_rot[i])^2)
      keep <- d <= radius
      data.frame(
        mg_barcode = mg_df$barcode[i],
        mg_type    = mg_label,
        n_neurons  = sum(keep),
        total_dist = sum(d[keep]),
        mean_dist  = if (any(keep)) mean(d[keep]) else NA)
    }))
  }
  bind_rows(calc_stats(mg_pos, neuron_pos, "MG_Igf1pos"),
            calc_stats(mg_neg, neuron_pos, "MG_Igf1neg"))
}

# Run across all radii
results_all <- bind_rows(lapply(names(radii), function(rname) {
  res <- all_output %>% group_by(sample) %>%
    group_modify(~ analyze_neighbors(.x, radius = radii[rname])) %>%
    ungroup()
  res$radius_label <- rname; res
}))

results_all <- results_all %>%
  mutate(genotype = ifelse(grepl("^WT", sample), "WT", "5xFAD"),
         radius_label = factor(radius_label,
           levels = c("25um", "50um", "100um", "150um")))

# Sample-level summary (n = 3 per genotype)
sample_summary <- results_all %>%
  filter(n_neurons > 0) %>%
  group_by(radius_label, genotype, sample, mg_type) %>%
  summarise(mean_n_neurons = mean(n_neurons),
            mean_dist      = mean(mean_dist, na.rm = TRUE),
            .groups = "drop")

# Wilcoxon: WT vs 5xFAD (sample-level, then cell-level)
stat_sample <- sample_summary %>%
  group_by(radius_label, mg_type) %>%
  wilcox_test(mean_n_neurons ~ genotype) %>% add_significance()
stat_cell   <- results_all %>%
  group_by(radius_label, mg_type) %>%
  wilcox_test(n_neurons ~ genotype) %>% add_significance()

# Plot template — sample-level dots + cell-level violin per genotype
plot_dots <- function(df, mg_lab, ylab) {
  ggplot(df %>% filter(mg_type == mg_lab),
         aes(x = genotype, y = mean_n_neurons, color = genotype)) +
    geom_jitter(size = 3, width = 0.1) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.3, linewidth = 0.6) +
    facet_wrap(~ radius_label) +
    scale_color_manual(values = c(WT = "#3498DB", "5xFAD" = "#E74C3C")) +
    labs(y = ylab, x = NULL) +
    theme_bw() + theme(legend.position = "none")
}

ggsave("Igf1pos_count_sample.pdf",
       plot_dots(sample_summary, "MG_Igf1pos", "Mean neuron count"),
       width = 8, height = 5)
ggsave("Igf1neg_count_sample.pdf",
       plot_dots(sample_summary, "MG_Igf1neg", "Mean neuron count"),
       width = 8, height = 5)
```

### 5.11 Step 8 — Co-Localization Map (MG Igf1+ vs Pain-Related Igf1r+ Neurons)

A single spatial DimPlot overlays MG_Igf1+ (red) and pain-related
Igf1r+ neurons (green) on the WT_3 / 5xFAD_1 sections, with all other
cells greyed out.

```r
color_map <- c(
  "MG_Igf1+"      = "#EF7B77",
  "Neuron_Igf1r+" = "#84a494",
  "MG_Igf1-"      = "#E0E0E0",
  "Neuron_Igf1r-" = "#E0E0E0",
  "Other"         = "#E0E0E0")

flag_combined <- function(obj, pain_related) {
  d_igf1  <- FetchData(obj, vars = c("Igf1",  "cellname"))
  d_igf1r <- FetchData(obj, vars = c("Igf1r", "cellname"))

  cg <- rep("Other", ncol(obj))
  cg[d_igf1$cellname == "MG" & d_igf1$Igf1 == 0] <- "MG_Igf1-"
  cg[d_igf1$cellname == "MG" & d_igf1$Igf1  > 0] <- "MG_Igf1+"
  cg[d_igf1r$cellname %in% pain_related & d_igf1r$Igf1r == 0] <- "Neuron_Igf1r-"
  cg[d_igf1r$cellname %in% pain_related & d_igf1r$Igf1r  > 0] <- "Neuron_Igf1r+"
  obj$color_group <- cg
  obj
}

pbmc_ad1 <- flag_combined(pbmc_ad1, pain_related)
pbmc_wt3 <- flag_combined(pbmc_wt3, pain_related)

pdf("Neu_MG_distur_AD1.pdf", width = 6, height = 4)
DimPlot(pbmc_ad1, group.by = "color_group", reduction = "spatial",
        pt.size = 2, cols = color_map,
        order = c("MG_Igf1+", "Neuron_Igf1r+")) + NoAxes()
dev.off()

pdf("Neu_MG_distur_WT3.pdf", width = 6, height = 4)
DimPlot(pbmc_wt3, group.by = "color_group", reduction = "spatial",
        pt.size = 2, cols = color_map,
        order = c("MG_Igf1+", "Neuron_Igf1r+")) + NoAxes()
dev.off()
```

### 5.12 Step 9 — Per-Gene Neuron-Class Barplot

```r
library(Seurat); library(ggplot2); library(dplyr); library(tidyr)

gene_list <- setdiff(unique(c(
  "Sst","Vglut2","Pdyn","Maf","Nfib","Tac2","Vgat","Igf1r","Cst7",
  "Tac1","Npy","Reln","Gal","Cck","Nos1","Pvalb","Crh","Tmem119",
  "Tacr1","Pax2","Rorb","Grpr","Chat","Rreb1","Igf1","Hexb")),
  c("Tmem119", "Hexb", "Cst7"))

sub_obj <- subset(pbmc, celltype %in% c("Glut", "GABA", "ACh"))
sub_obj <- JoinLayers(sub_obj)
expr    <- GetAssayData(sub_obj, layer = "data")
genes_present <- gene_list[gene_list %in% rownames(expr)]

count_list <- list()
for (g in genes_present) {
  for (ct in c("Glut", "GABA", "ACh")) {
    cells <- colnames(sub_obj)[sub_obj$celltype == ct]
    count_list[[paste0(g, "_", ct)]] <- data.frame(
      Gene = g, CellType = ct, Count = sum(expr[g, cells] > 0))
  }
}
df <- do.call(rbind, count_list)

gene_total <- df %>% group_by(Gene) %>%
  summarise(Total = sum(Count)) %>% arrange(Total)
df$Gene     <- factor(df$Gene,     levels = gene_total$Gene)
df$CellType <- factor(df$CellType, levels = c("Glut", "GABA", "ACh"))

ggplot(df, aes(x = Gene, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = c(Glut = "#DA9087",
                               GABA = "#547298",
                               ACh  = "#84a494")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Cell count", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   size = 10, face = "italic"))

ggsave("neuron_gene_barplot.pdf", width = 8, height = 5)
```

