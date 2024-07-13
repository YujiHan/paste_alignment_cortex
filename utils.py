import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
import scipy.sparse as sp
import glob
from tqdm import tqdm
from paste import pairwise_align
import ot
import torch
import matplotlib.pyplot as plt


def process_spatial_data(adata):
    '''
    处理obsm['spatial']为指定格式
    '''

    adata.obsm['spatial_raw'] = adata.obsm[
        'spatial'
    ].copy()  # 复制obsm['spatial']中的数据
    adata.obsm['spatial'] = np.delete(
        adata.obsm['spatial'], [0, 1], axis=1
    )  # 删除obsm['spatial']中的前两列


process_spatial_data(slice_0)
process_spatial_data(slice_1)


def random_subset(adata, cell_num=5000):
    '''
    随机选取 3000 个细胞
    '''

    random_indices = np.random.choice(
        adata.n_obs, cell_num, replace=False
    )  # 从 276593 个细胞中随机选取 3000 个细胞的索引
    adata_subset = adata[random_indices, :]  # 使用选取的索引创建新的 AnnData 对象

    return adata_subset


slice_0_subset = random_subset(slice_0)
slice_1_subset = random_subset(slice_1)


def get_slice_HVG(slice, n_top_genes=2000, need_norm_and_log1p=False, target_sum=1e4):
    '''
    计算并选择前2000个高可变基因
    '''

    slice.var_names_make_unique()  # 确保变量名唯一，避免重复
    sc.pp.filter_genes(slice, min_counts=5)  # 过滤掉计数小于5的基因
    sc.pp.highly_variable_genes(
        slice, n_top_genes=n_top_genes, flavor='seurat_v3'
    )  # 使用seurat_v3方法计算高可变基因，并选择前3000个
    slice = slice[:, slice.var.highly_variable]  # 只保留高可变基因

    if need_norm_and_log1p:
        sc.pp.normalize_total(
            slice, target_sum=target_sum
        )  # 将每个细胞的总表达量归一化为10000
        sc.pp.log1p(slice)  # 对数据进行对数变换

    return slice  # 返回处理后的数据


slice_0_subset_HVG = get_slice_HVG(slice_0_subset)
slice_1_subset_HVG = get_slice_HVG(slice_1_subset)


'''
用循环获取对齐的spot
'''
# count = 0
# for i,arr in enumerate(s1_s2_pi):
#     for j,pi in enumerate(arr):
#         if pi !=0:
#             print(f'({i}, {j}):\t\t{pi}')
#             count+=1

# print('count: {count}')


def verify_track(cur_index=42):
    '''
    验证查找轨迹的函数是否正确
    '''
    print(cur_index)
    for i in range(7):
        cur_index = all_pairwise_index[i][1, cur_index]
        print(cur_index)


verify_track()
