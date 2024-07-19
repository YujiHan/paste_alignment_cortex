import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
import scipy.sparse as sp
from tqdm import tqdm
from paste import pairwise_align
import ot
import torch
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import anndata as ad

dataset_path = '/home/hanyuji/Workbench/ST/ST_data_check/paste_alignment_cortex/PASTE_align/cortex_macaque1_119slice_subset_5000spot_2000gene_HVG.h5ad'
adata = sc.read_h5ad(dataset_path)

unique_batches = adata.obs['batch'].unique()
adata_list = [adata[adata.obs['batch'] == batch].copy() for batch in unique_batches]


def get_edge_index(M):
    # 将稠密矩阵M转换为稀疏矩阵，并获得边的索引edge_index
    M_ = coo_matrix(M)  # 稠密矩阵--->稀疏矩阵
    values = M_.data
    indices = np.vstack((M_.row, M_.col))  # 我们真正需要的coo形式
    return indices


torch.cuda.set_device(1)  # 将设备设置为 'cuda:0'

import pickle

dir = '/home/hanyuji/Workbench/ST/ST_data_check/paste_alignment_cortex/PASTE_align/'  # 保存的位置


all_pairwise_index = []
for i in tqdm(range(len(adata_list) - 1)):
    pi = pairwise_align(
        sliceA=adata_list[i],
        sliceB=adata_list[i + 1],
        backend=ot.backend.TorchBackend(),
        use_gpu=True,
        numItermax=50000,
    )
    edge_index = get_edge_index(pi)
    all_pairwise_index.append(edge_index)


'''保存变量'''
with open(dir + 'pairwise_index_5000spot_0715.pkl', 'wb') as file:
    pickle.dump(all_pairwise_index, file)
print('all_pairwise_index saved')


def track_sequence(start, arrays):
    # start是起始数字，arrays是数组列表
    path = [start]  # 轨迹起始
    current = start
    for arr in arrays:
        idx = np.where(arr[0] == current)[0][0]  # 找到current在第一维的位置
        current = arr[1][idx]  # 获取该位置在第二维的对应数字
        path.append(current)  # 添加到路径中
    return path


# 跟踪所有起始数字的轨迹
paths = np.array([track_sequence(i, all_pairwise_index) for i in range(5000)])


'''保存变量'''
np.save(dir + 'track_sequence_5000spot_0715.pkl', paths)
print('paths saved')


paths[42]
