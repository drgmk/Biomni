# Single Cell RNA-seq Quality Control

---

## Metadata

**Short Description**: Best practices for quality control in single-cell RNA-seq data, including filtering and doublet detection.

**Authors**: GMK

**Affiliations**: Malaghan

**Version**: 1.0

**Last Updated**: December 2025

**License**: CC BY 4.0

**Commercial Use**: ✅ Allowed

**Source**: https://github.com/drgmk/scrna

**Citation**: GMK

---

## Overview

Quality control is essential in single-cell RNA-seq analysis, and is (once a counts matrix is generated) the first step before any downstream analysis. This guide covers a standard workflow for QC, including filtering low-quality cells and detecting doublets.

## Quality Control Steps

### 1. Filtering low-quality genes and cells
Set initial thresholds for:
- Minimum genes detected per cell (e.g., >200)
- Minimum number of cells expressing a gene (e.g., >3)
- Maximum mitochondrial gene percentage (e.g., <5-10%)
- Maximum fraction of counts from the most expressed gene (e.g., <20-30%)
- Assess and revise thresholds if necessary.
- Visualise results before and after filtering.
- Save filtering masks and thresholds in `adata.uns` for documentation.

**Tools**: scrna.functions, scanpy, rapids-singlecell
**Best for**: All datasets

### 2. Doublet Detection
- Detect and remove doublets to prevent misleading annotations.
- Visualize doublet scores and predictions to assess results.
- Save doublet results in `adata.obs` for later exclusion.

**Tools**: scanpy.pp.scrublet, rapids-singlecell.pp.scrublet
**Best for**: All datasets

### 3. Save filtered dataset
- Apply filtering and doublet masks to the `AnnData` object.
- Save the filtered dataset for downstream analysis.

## Recommended Workflow

### Step 1: Initial Filtering
- Set initial filtering thresholds based on dataset characteristics.
- Apply trimming function `trim_outliers` to remove outliers in n_genes_by_counts vs. total_counts space.
  - Retain cells within a percentile range (e.g., 1-99th percentile) of residuals.
  - Apply masks for:
    - mitochondrial percentage
    - fraction of counts from the most expressed gene
    - minimum genes per cell
    - minimum cells per gene
  - Apply filtering per sample/batch if applicable.
- Alter thresholds and repeat if necessary to retain sufficient high-quality cells.

Functions to perform these steps are in the `scrna.functions` module. Docstrings and usage are below:

```python
# docstrings for QC functions in scrna.functions
def compute_qc_metrics(adata, extra_genes=[]):
    """Calculate QC metrics for RNA data.

    Annotated data matrix `adata` is updated in-place with percent counts for mitochondrial, ribosomal,
    and MALAT genes, as well as any extra gene prefixes provided.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing RNA expression data.
    extra_genes : list, optional
        List of extra gene prefixes to calculate percent counts for.
    """


def trim_outliers(
    adata,
    x="total_counts",
    y="n_genes_by_counts",
    groupby=None,
    extra_mask=None,
    extra_mask_boolean=None,
    pct=100.0,
):
    """Function to fit a line in log space, trim outliers, and return boolean mask.

    Parameters
    ----------
    x : array-like
        Independent variable.
    y : array-like
        Dependent variable.
    groupby : str, optional
        Column name in `adata.obs` to group by before fitting line and trimming outliers.
    extra_mask : dict, optional
        Dictionary specifying additional masks to apply before trimming outliers.
        Format is {column_name: (threshold, 'min' or 'max')}.
    extra_mask_boolean : array-like, optional
        Boolean mask to apply before trimming outliers.
    pct : int, optional
        Percentile to use for trimming outliers, default is 100 (no trimming).
    """

def plot_gene_counts(
    adata,
    hue="sample",
    mask=None,
    order=None,
    show_masked=True,
    colour_by="pct_counts_in_top_1_genes",
    size_by="pct_counts_mt",
):
    """Plot gene counts and mitochondrial fraction for each sample.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing RNA expression data.
    hue : str, optional
        Column name in `adata.obs` to use for multiple panels.
    mask : array-like, optional
        Boolean mask to filter the data before plotting.
        This could come from `trim_outliers`.
    order : list, optional
        List of group names in the order to plot. If None, uses the order in `adata.obs[hue].unique()`.
    show_masked : bool, optional
        Whether to show the masked (filtered out) points in light grey.
    colour_by : str, optional
        Column name in `adata.obs` to use for coloring the points.
    size_by : str, optional
        Column name in `adata.obs` to use for the size of the points.
    """
```

```python
# use scrna.functions to filter low-quality cells
import scrna.functions as scfunc
import scanpy as scanpy

# assume rna is an AnnData object with raw counts loaded
# and per-cell metadata including sample/batch IDs in rna.obs
min_genes = 200
min_cells = 3
max_mt_pct = 10.0
max_top1_pct = 20.0
pct_outlier_cutoff = 99.0
sample_col = "sample_id"  # or batch column if applicable

mask_cells, _ = scanpy.pp.filter_cells(rna, min_genes=min_genes, inplace=False)
mask_genes, _ = scanpy.pp.filter_genes(rna, min_cells=min_cells, inplace=False)
compute_qc_metrics(rna)

mask = scfunc.trim_outliers(
    rna,
    groupby=sample_col,
    extra_mask={
        "pct_counts_mt": [max_mt_pct, "max"],
        "pct_counts_in_top_1_genes": [max_top1_pct, "max"],
    },
    extra_mask_boolean=mask_cells,
    pct=pct_outlier_cutoff,
)
```
- Visualise basic QC metrics, e.g. violin plots.
- Visualize gene counts per sample using `plot_gene_counts` to assess filtering impact.

```python
# visualize basic QC metrics
import scanpy as sc

# plot gene counts per sample
fig = plot_gene_counts(
    rna, hue=sample_col, mask=mask, show_masked=True
)
fig.savefig(str(out_path / "gene_counts_per_sample.png"), dpi=150)

sc.pl.violin(
    rna,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby=sample_col,
    jitter=0.4,
    multi_panel=True,
)
```

- Save threshold values in `adata.uns` for documentation.

```python
# save and apply masks and thresholds
rna.uns["meta_qc_thresholds"] = {
    "min_genes_per_cell": min_genes,
    "min_cells_per_gene": min_cells,
    "max_mt_pct": max_mt_pct,
    "max_top1_pct": max_top1_pct,
    "pct_outlier_cutoff": pct_outlier_cutoff,
}
rna.uns["meta_qc_mask_cells"] = mask
rna.uns["meta_qc_mask_genes"] = mask_genes
rna = rna[mask, mask_genes].copy()```
```

### Step 2: Doublet Detection

- Choose appropriate doublet detection tool based on computational resources.
  - If rapid-singlecell is available and GPU/CUDA is accessible, use `rapids_singlecell.pp.scrublet`.
  - Otherwise, use `scanpy.pp.scrublet`.
  - Visualize doublet scores and predictions to assess results.
  - Save doublet results in `adata.obs` for later exclusion.

```python
# Scanpy example
import scanpy as sc

sc.tl.scrublet(adata, batch_key='batch')  # if multiple batches
```

```python
# rapids-singlecell example
import rapids_singlecell as rsc

rsc.tl.scrublet(adata, batch_key='batch')  # if multiple batches
```

```python
# discard predicted doublets
doublet_mask = rna.obs["predicted_doublet"] == False
rna = rna[doublet_mask].copy()
```

```python
# visualize doublet scores and predictions
sc.pl.violin(
    rna,
    ["doublet_score", "predicted_doublet"],
    groupby=sample_col,
    jitter=0.4,
    multi_panel=True,
)
```

### Step 3: Save filtered dataset
- Save the filtered dataset for downstream analysis.

```python
# save filtered dataset
rna.write_h5ad(str(out_path / "rna_filtered_qc.h5ad"))
```

## Best Practices

### Do's:
1. **Use the recommended workflow for QC** - Follow the steps outlined above
2. **Check that thresholds are reasonable** - Avoid over-filtering
3. **Visualize QC metrics** - Use plots to guide threshold selection
4. **Apply doublet detection to filtered data** - Reduces false positives
5. **Keep results for downstream analysis** - Save dataset after filtering and doublet detection
6. **Document thresholds used** - Save masks and thresholds in `adata.uns` for reproducibility

### Don'ts:
1. **Don't be too aggressive** - A small number of low-quality cells is okay
2. **Don't forget to save the results** - Always save the filtered dataset and metadata for downstream analysis
3. **Don't save too many files** - Save a final filtered dataset rather than multiple intermediate files

## Common Pitfalls

### 1. Quality may vary by sample/batch
**Problem**: Uniform thresholds may not fit all samples
**Solution**: Apply filtering per sample/batch

### 2. Filtering is too strict/lenient
**Problem**: Over/under-filtering affects downstream analysis
**Solution**: Visualize QC metrics, adjust thresholds iteratively
 before transfer

## Tool Selection Guide

| Scenario | Recommended Tool | Why |
|----------|------------------|-----|
| All QC | Functions in this document | Two-dimensional filtering |
| GPU/CUDA | rapids-singlecell | Much faster processing |
| CPU | scanpy | Slower but will work |

## Quality Control Checklist

- [ ] Check dataset for per-cell metadata, e.g. batch/sample IDs in `adata.obs`
- [ ] Set initial filtering thresholds based on dataset characteristics
- [ ] Apply trimming function to remove outliers in n_genes_by_counts vs. total_counts space
- [ ] Apply masks for mitochondrial percentage, fraction of counts from most expressed gene, minimum genes per cell, minimum cells per gene
- [ ] Apply filtering per sample/batch if applicable
- [ ] Alter thresholds and repeat if necessary to retain sufficient high-quality cells
- [ ] Detect doublets in filtered data for later exclusion using appropriate tool
- [ ] Save threshold values and masks in `adata.uns` for documentation
- [ ] Check that doublet results are saved in `adata.obs`
- [ ] Save filtered dataset for downstream analysis
- [ ] Visualize QC metrics before and after filtering, both 1d distributions with thresholds and 2d plots

## Resources

### Tools:
- **above code snippets**: see code snippets above
- **scrna**: https://github.com/drgmk/scrna
- **scanpy**: https://scanpy.readthedocs.io/
- **rapids-singlecell**: https://github.com/scverse/rapids_singlecell/

## Troubleshooting

**Issue**: Not enough cells/genes retained after filtering
→ Relax thresholds, check QC metric distributions
