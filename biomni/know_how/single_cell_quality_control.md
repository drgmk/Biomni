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
Remove cells that do not meet quality thresholds. Revise thresholds if necessary.

**Tools**: scanpy, rapids-singlecell
**Best for**: All datasets

### 2. Doublet Detection
Detect (and later remove) doublets to prevent misleading annotations.

**Tools**: scanpy.pp.scrublet, rapids-singlecell.pp.scrublet
**Best for**: All datasets

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
- Save threshold values in `adata.uns` for documentation.

```python
# useful functions for QC
import scanpy as sc
import numpy as np
import scipy.stats

def compute_qc_metrics(adata, extra_genes=[]):
    """Calculate QC metrics for RNA data.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing RNA expression data.
    extra_genes : list, optional
        List of extra gene prefixes to calculate percent counts for.
    """

    pct_counts = {
        "mt": ["MT-"],
        "ribosomal": ["RPL", "RPS"],
        "malat": ["MALAT"],
    }
    for g in extra_genes:
        pct_counts[g] = [g]

    # convert everthing to lower case for matching
    for k, v in pct_counts.items():
        pct_counts[k] = [s.lower() for s in v]

    for k in pct_counts.keys():
        adata.var[k] = adata.var_names.str.lower().str.startswith(pct_counts[k][0])
        if len(pct_counts[k]) > 1:
            for s in pct_counts[k][1:]:
                adata.var[k] = np.logical_or(
                    adata.var[k], adata.var_names.str.lower().str.startswith(s)
                )

        sc.pp.calculate_qc_metrics(
            adata, qc_vars=[k], percent_top=None, log1p=False, inplace=True
        )

    # also calculate percent in top gene
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=(), percent_top=[1], log1p=False, inplace=True
    )

    # add meta to adata.uns
    adata.uns["pct_counts"] = pct_counts


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

    mask = np.ones(adata.shape[0], dtype=bool)

    if groupby is not None:
        for g in np.unique(adata.obs[groupby]):
            mask_g = adata.obs[groupby] == g
            mask[mask_g] = trim_outliers(
                adata[mask_g, :], x=x, y=y, extra_mask=extra_mask, pct=pct
            )

        adata.uns["trim_outliers_mask"] = mask
        return mask

    if extra_mask is None:
        extra_mask_ = np.ones(adata.shape[0], dtype=bool)
    else:
        extra_mask_ = np.ones(adata.shape[0], dtype=bool)
        for k, v in extra_mask.items():
            if v[1] == "max":
                extra_mask_ = np.logical_and(extra_mask_, adata.obs[k] < v[0])
            elif v[1] == "min":
                extra_mask_ = np.logical_and(extra_mask_, adata.obs[k] > v[0])
            else:
                raise ValueError(
                    f'unknown key {k} in extra_mask, must contain "min" or "max"'
                )

    if extra_mask_boolean is not None:
        extra_mask_ = np.logical_and(extra_mask_, extra_mask_boolean)

    x_ = np.log10(adata.obs[x])
    y_ = np.log10(adata.obs[y])
    fit = scipy.stats.linregress(x_[extra_mask_], y_[extra_mask_])
    y_fit = fit.intercept + fit.slope * x_
    resid = y_ - y_fit
    thresh = np.percentile(resid, [100 - pct, pct])
    mask = np.logical_and(resid > thresh[0], resid < thresh[1])
    adata.uns["trim_outliers_mask"] = mask
    return np.logical_and(mask, extra_mask_)


def plot_nxy(n):
    """Return no of x, y panels for plotting approx square panels."""
    if n < 4:
        y = 1
    elif n < 9:
        y = 2
    elif n < 16:
        y = 3
    else:
        y = 4  # ok up to n=24
    x = int(np.ceil(n / y))
    return x, y


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

    if mask is None:
        mask = np.ones(adata.shape[0], dtype=bool)

    if order is None:
        order = adata.obs[hue].unique()

    vmax = (
        np.max(adata[mask].obs[colour_by]) if colour_by in adata.obs.columns else None
    )

    nx, ny = plot_nxy(len(order))
    fig, ax = plt.subplots(ny, nx, sharey=True, sharex=True, figsize=(10, 7))

    # check on size_by, and rescale between 1 and 5
    sizes = adata.obs[size_by] / 4
    if sizes.min() == sizes.max():
        sizes = 2 * np.ones(adata.shape[0])
    else:
        sizes = 1 + 4 * (sizes - sizes.min()) / (sizes.max() - sizes.min())

    for i, s in enumerate(order):
        a = ax.flatten()[i]
        ok = (adata.obs[hue] == s) & mask
        tmp = adata[ok, :]
        _ = a.scatter(
            tmp.obs["total_counts"],
            tmp.obs["n_genes_by_counts"],
            s=sizes[ok],
            c=tmp.obs[colour_by],
            vmin=0,
            vmax=vmax,
            cmap="viridis",
        )
        if show_masked:
            ok = (adata.obs[hue] == s) & np.invert(mask)
            tmp = adata[ok, :]
            a.scatter(
                tmp.obs["total_counts"],
                tmp.obs["n_genes_by_counts"],
                s=sizes[ok],
                c="lightgrey",
                alpha=0.5,
                zorder=-1,
            )
            #   s=tmp.obs[size_by]/4, c=tmp.obs[colour_by], alpha=0.2,
            #   vmin=0, vmax=vmax, cmap='Grays', zorder=-1)

        a.set_title(s)

    [a.set_visible(False) for a in ax.flatten()[i + 1 :]]
    if ny == 1:
        ax = ax[np.newaxis, :]
    ax[ny - 1, 0].set_ylabel("n_genes_by_counts")
    ax[ny - 1, 0].set_xlabel("total_counts")
    ax[0, 0].set_xscale("log")
    ax[0, 0].set_yscale("log")

    # Place colorbar to the right of all axes, spanning full height
    # divider = make_axes_locatable(ax.flatten()[-1])
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.tight_layout()
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes((0.9, 0.15, 0.02, 0.7))
    cb = fig.colorbar(_, cax=cbar_ax, aspect=30)
    cb.set_label(colour_by)
    # turn grid on for all axes
    for a in ax.flatten():
        a.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.5)
        a.set_axisbelow(True)
    return fig
```

```python
# use the above functions
min_genes = 200
min_cells = 3
max_mt_pct = 10.0
max_top1_pct = 30.0
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

# plot gene counts per sample
fig = plot_gene_counts(
    rna, hue=sample_col, mask=mask, show_masked=True
)
fig.savefig(str(out_path / "gene_counts_per_sample.pdf"))

# save and apply mask
rna.uns["meta_qc_mask_cells"] = mask
rna.uns["meta_qc_mask_genes"] = mask_genes
rna = rna[mask, mask_genes].copy()```
```

### Step 2: Doublet Detection

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

## Best Practices

### Do's:
1. **Check that thresholds are reasonable** - Avoid over-filtering

### Don'ts:
1. **Don't be too aggressive** - A small number of low-quality cells is okay

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

- [ ] Check dataset for per-cell metadata, e.g. batch/sample IDs
- [ ] Set initial filtering thresholds based on dataset characteristics
- [ ] Apply trimming function to remove outliers in n_genes_by_counts vs. total_counts space
- [ ] Apply masks for mitochondrial percentage, fraction of counts from most expressed gene, minimum genes per cell, minimum cells per gene
- [ ] Apply filtering per sample/batch if applicable
- [ ] Alter thresholds and repeat if necessary to retain sufficient high-quality cells
- [ ] Save threshold values in `adata.uns` for documentation
- [ ] Detect doublets for later exclusion using appropriate tool
- [ ] Visualize QC metrics before and after filtering

## Resources

### Tools:
- **functions in this document**: see code snippets above
- **scanpy**: https://scanpy.readthedocs.io/
- **rapids-singlecell**: https://github.com/scverse/rapids_singlecell/

## Troubleshooting

**Issue**: Not enough cells/genes retained after filtering
→ Relax thresholds, check QC metric distributions
