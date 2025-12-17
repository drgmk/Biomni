# Single Cell RNA-seq Clustering

---

## Metadata

**Short Description**: Best practices for dimensionality reduction and clustering in single-cell RNA-seq data.

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

This guide walks through a practical scRNA-seq clustering workflow, covering dimensionality reduction, neighborhood graph, UMAP, Leiden clustering, and batch correction. It assumes the data have undergone quality control already.

## Normalization and Clustering Steps

### 1. Normalization, Log Transformation, and Scaling
Normalize counts to a target sum (e.g., 1e4) and log-transform the data. Scaling (zero-centering and variance scaling) is recommended before PCA. But zero-centering may need to be skipped for large datasets to avoid densifying the data matrix. Memory management can be enabled to handle larger datasets on GPUs.

### 2. Highly Variable Genes and PCA
Identify highly variable genes and perform PCA for dimensionality reduction on these genes.

### 3. Batch Correction
Use Harmony for batch correction if multiple samples or batches are present. This is applied on the PCA embeddings. Separation between batches can be assessed with the Average Silhouette Width (ASW) metric.

### 4. Neighborhood Graph and UMAP
Construct a neighborhood graph using the batch-corrected PCA embeddings. Compute UMAP for visualization. Several variables can be tuned; `n_neighbors` for neighborhood size, `min_dist` as for cluster structure vs spread.

### 5. Leiden Clustering
Perform Leiden clustering on the neighborhood graph. The `resolution` parameter controls cluster granularity.

### 6. Visualization
Generate UMAP plots colored by sample and cluster assignments for visual inspection.

## Recommended Workflow

Below is the basic pipeline for normalized clustering. This assumes the `AnnData` object (`rna`) is already QC‑filtered. If it is not (e.g. `rna.obs` has no doublet-related columns, `rna.uns` has no reference to QC meta/variables), refer to the [Single Cell RNA-seq Quality Control](single_cell_quality_control.md) document.

```python
import gc
import anndata as ad
import scanpy as sc
import scib.metrics
# import rapids_singlecell as rsc  # if GPU available, rsc can be used in place of sc

# Enable managed memory with RMM for rapids-singlecell (GPU)
# import rmm
# import cupy as cp
# from rmm.allocators.cupy import rmm_cupy_allocator
# rmm.reinitialize(managed_memory=True, pool_allocator=False)
# cp.cuda.set_allocator(rmm_cupy_allocator)

rna = ad.read_h5ad('rna_qc_filtered.h5ad')

# move data to GPU for rapids-singlecell
# rsc.get.anndata_to_gpu(rna)  # only for rapids-singlecell

# Normalize counts and log‑transform
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.scale(rna, zero_center=True, max_value=10)  # zero‑centered scaling

# Highly variable genes and PCA
sc.pp.highly_variable_genes(rna)
sc.tl.pca(rna)

# Batch correction with Harmony (samples in rna.obs['sample'])
# for rapids-singlecell, use: rsc.pp.harmony_integrate(...)
sc.external.pp.harmony_integrate(rna, key='sample')

# Neighborhood graph on Harmony embedding
sc.pp.neighbors(rna, n_neighbors=20, use_rep='X_pca_harmony')

# UMAP and Leiden clustering
sc.tl.umap(rna, min_dist=0.3, random_state=42)
sc.tl.leiden(rna, resolution=0.8)

# Optional: additional UMAPs for visual comparison
sc.tl.umap(rna, min_dist=0.15, random_state=42, key_added='X_umap_half_0.3_')
sc.tl.umap(rna, min_dist=0.6,  random_state=42, key_added='X_umap_twice_0.3_')

# Average Silhouette Width (ASW) for batch mixing assessment
asw_batch = scib.metrics.silhouette_batch(
    rna,
    batch_key='sample',
    label_key='leiden',
    embed='X_pca_harmony'
)
print(f'Average Silhouette Width (ASW) for batch mixing: {asw_batch:.4f}')

# get data back to CPU if using rapids-singlecell
# rsc.get.anndata_to_cpu(rna)  # only for rapids-singlecell
# gc.collect()  # tidy up to save some memory

# save the clustered data
rna.write('rna_clustered.h5ad')
```

## Visualize Outputs

Generate UMAPs for cluster inspection and save figures:

```python
import matplotlib.pyplot as plt

# Highly variable gene selection
sc.pl.highly_variable_genes(rna, save='highly_variable_genes.png')

# PCA variance ratio and first few PCs
sc.pl.pca_variance_ratio(rna, n_pcs=50, log=True, save='pca_variance_ratio.png')
sc.pl.pca(
    rna,
    color=["sample", "leiden"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
    save='pca_samples_leiden.png'
)

# Basic overview UMAPs by sample and cluster
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
sc.pl.umap(rna, color='sample', ax=ax[0], show=False)
sc.pl.umap(rna, color='leiden', ax=ax[1], show=False)
fig.tight_layout()
fig.savefig('umap_sample_leiden.png', dpi=150)

# UMAPs with different min_dist for structure vs spread
fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
sc.pl.embedding(rna, 'X_umap_half_0.3_', color='leiden', ax=ax[0], show=False)
sc.pl.embedding(rna, 'X_umap_twice_0.3_', color='leiden', ax=ax[1], show=False)
fig.tight_layout()
fig.savefig('umap_leiden_min_dist.png', dpi=150)

# UMAP with doublet annotations if available
sc.pl.umap(
    rna,
    color=["leiden", "predicted_doublet", "doublet_score"],
    wspace=0.5,
    size=3,
    save='umap_doublets.png'
)

# Combined UMAP with QC metrics
sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    save='umap_qc_metrics.png'
)
```

## Best Practices

### Do's:
1. **Use the recommended workflow** - Follow the steps outlined above
2. **Visualize the outputs** - Use plots to guide parameter selection
3. **Keep results for downstream analysis** - Save dataset after filtering and doublet detection
4. **Document parameters used** - Save parameters in `adata.uns` for reproducibility

### Don'ts:
1. **Don't be too aggressive** - A small number of low-quality cells is okay
2. **Don't forget to save the results** - Always save the processed data and metadata for downstream analysis

## Common Pitfalls

### 1. Default values may not separate clusters well
**Problem**: Cells may not be well separated with default parameters
**Solution**: Tune `n_neighbors`, `min_dist`, and `leiden` resolution based on dataset size and desired granularity


## Tool Selection Guide

| Scenario | Recommended Tool | Why |
|----------|------------------|-----|
| GPU/CUDA | rapids-singlecell | Much faster processing |
| CPU | scanpy | Slower but will work |

## Clustering Checklist

- [ ] Ensure data have undergone quality control.
- [ ] Check dataset for per-cell metadata, e.g. batch/sample IDs in `adata.obs`
- [ ] Normalize, log-transform, and scale data.
- [ ] Identify highly variable genes and perform PCA.
- [ ] Apply batch correction if multiple samples/batches
- [ ] Construct neighborhood graph using batch-corrected PCA embeddings
- [ ] Compute UMAP for visualization
- [ ] Perform Leiden clustering
- [ ] Visualize UMAPs colored by sample and cluster assignments

## Resources

### Tools:
- **scanpy**: https://scanpy.readthedocs.io/
- **rapids-singlecell**: https://github.com/scverse/rapids_singlecell/

## Troubleshooting

**Issue**: Clusters are not well separated in UMAP.
**Solution**: Adjust `n_neighbors` and `min_dist` parameters in UMAP, and try different `resolution` values in Leiden clustering.

**Issue**: Memory errors during processing.
**Solution**: Try scaling with `zero_center=False` to avoid densifying the data matrix, or enable memory management. For a 48GB GPU, up to 0.5M cells can typically be processed without issues.