import torch
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import pickle
import anndata as ad

# dataset_path = '/home/hanyuji/Workbench/ST/paste_alignment_cortex/PASTE_align/marmoset_hippocampus_cell1_27slices_1000spot_2000gene.h5ad'
dataset_path = '/home/hanyuji/Workbench/ST/paste_alignment_cortex/PASTE_align/macaque_hippocampus_cell_macaque1_30slices_1000spot_2000gene.h5ad'
adata = sc.read_h5ad(dataset_path)
adata
fig_name = 'macaque_hippocampus_cell1_30slices_pca'
batch = 'batch'

min_dist = 0.1
spread = 0.3
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_pcs=40, n_neighbors=30, use_rep='X')

sc.tl.umap(adata, min_dist=min_dist, spread=spread)
sc.pl.umap(
    adata,
    color=[batch],
    legend_fontsize=5,
    ncols=2,
    save=f'_{fig_name}_{min_dist}_{spread}',
)

print(fig_name)
