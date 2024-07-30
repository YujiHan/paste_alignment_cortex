import scanpy as sc
import numpy as np
import pandas as pd
import stlearn as st
import SpaGCN as spg

def load_slices_fourST(data_dir='/home/sxa/Datasets/breast_cancer_data/', slice_names=["slice1", "slice2", "slice3", "slice4"]):
    slices = []  
    for slice_name in slice_names:
        slice_i = sc.read_csv(data_dir + slice_name + ".csv")
        slice_i_coor = np.genfromtxt(data_dir + slice_name + "_coor.csv", delimiter = ',')
        slice_i.obsm['spatial_coor'] = slice_i_coor
        # Preprocess slices
        sc.pp.filter_genes(slice_i, min_counts = 15)
        sc.pp.filter_cells(slice_i, min_counts = 100)
        slices.append(slice_i)
    adata_layer_1, adata_layer_2, adata_layer_3, adata_layer_4 = slices
    common_genes = intersect(adata_layer_1.var.index, adata_layer_2.var.index)  # 筛选出共有的基因
    common_genes = intersect(common_genes, adata_layer_3.var.index)
    common_genes = intersect(common_genes, adata_layer_4.var.index)
    adata_layer_1 = adata_layer_1[:, common_genes]
    adata_layer_2 = adata_layer_2[:, common_genes]
    adata_layer_3 = adata_layer_3[:, common_genes]
    adata_layer_4 = adata_layer_4[:, common_genes]
    slices = [adata_layer_1, adata_layer_2, adata_layer_3, adata_layer_4]
    return slices

# def load_Human_DLPFC(data_dir='/home/sxa/Datasets/Human_DLPFC/', slice_name=151673):
#     # '/home/sxa/Datasets/Human_DLPFC/151673/151673_filtered_feature_bc_matrix.h5'
#     slice = sc.read_10x_h5(data_dir+str(slice_name)+'/'+str(slice_name)+'_filtered_feature_bc_matrix.h5')
#     coldata = pd.read_csv(data_dir+str(slice_name)+'/'+'col_data_'+str(slice_name)+'.csv')
#     spatial_coor = (coldata.loc[:, ['row', 'col']]).values
#     image_coor = (coldata.loc[:, ['imagerow', 'imagecol']]).values
#     spatialLIBD_cluster = (coldata.loc[:,['Cluster']]).values
#     real_label = (coldata.loc[:, ['layer_guess_reordered_short']]).values
#     slice.obsm['spatial'] = spatial_coor
#     slice.obsm['images_coordinates'] = image_coor
#     slice.obsm['real_label'] = real_label.reshape(-1)
#     slice.obsm['spatialLIBD_cluster'] = spatialLIBD_cluster.reshape(-1)
#     return slice

def load_Human_DLPFC(data_dir='/home/sxa/Datasets/Human_DLPFC/', slice_name=151673):
    slice = sc.read_visium(path=data_dir+str(slice_name)+'/', count_file=str(slice_name)+'_filtered_feature_bc_matrix.h5')
    coldata = pd.read_csv(data_dir+str(slice_name)+'/'+'col_data_'+str(slice_name)+'.csv')
    spatial_coor = (coldata.loc[:, ['row', 'col']]).values
    spatialLIBD_cluster = (coldata.loc[:,['Cluster']]).values
    real_label = (coldata.loc[:, ['layer_guess_reordered_short']]).values
    slice.obsm['spatial_coor'] = spatial_coor
    slice.obs['Ground Truth'] = real_label.reshape(-1)
    slice.obs['spatialLIBD_cluster'] = spatialLIBD_cluster.reshape(-1)
    var = list(set(slice.var.index))
    slice.var_names_make_unique()
    slice = slice[:, var]
    sc.pp.filter_genes(slice, min_counts = 5)
    sc.pp.highly_variable_genes(slice, n_top_genes=3000 ,flavor='seurat_v3')
    slice = slice[:, slice.var.highly_variable]
    # 标准化 log (验证标准化的有效性)
    sc.pp.normalize_total(slice, target_sum=1e4)
    # sc.pp.normalize_total(slice)
    sc.pp.log1p(slice)
    return slice

def load_Human_DLPFC_stLearn(data_dir='/home/sxa/Datasets/Human_DLPFC/', slice_name=151673):
    slice = st.Read10X(path=data_dir+str(slice_name)+'/', count_file=str(slice_name)+'_filtered_feature_bc_matrix.h5')
    coldata = pd.read_csv(data_dir+str(slice_name)+'/'+'col_data_'+str(slice_name)+'.csv')
    spatial_coor = (coldata.loc[:, ['row', 'col']]).values
    spatialLIBD_cluster = (coldata.loc[:,['Cluster']]).values
    real_label = (coldata.loc[:, ['layer_guess_reordered_short']]).values
    slice.obsm['spatial_coor'] = spatial_coor
    slice.obs['Ground Truth'] = real_label.reshape(-1)
    slice.obs['spatialLIBD_cluster'] = spatialLIBD_cluster.reshape(-1)
    var = list(set(slice.var.index))
    slice.var_names_make_unique()
    slice = slice[:, var]
    sc.pp.highly_variable_genes(slice, n_top_genes=3000 ,flavor='seurat_v3')
    slice = slice[:, slice.var.highly_variable]
    st.pp.filter_genes(slice,min_cells=1)
    st.pp.normalize_total(slice)
    st.pp.log1p(slice)

    return slice

def intersect(lst1, lst2): 
    """
    Gets and returns intersection of two lists.

    Args:
        lst1: List
        lst2: List
    
    Returns:
        lst3: List of common elements.
    """

    temp = set(lst2)
    print(len(temp))
    lst3 = [value for value in lst1 if value in temp]
    return lst3

# data_dir = '/home/sxa/Datasets/breast_cancer_data/'

def load(data_dir='/home/sxa/Datasets/Human_DLPFC/', slice_name=151673):
    # '/home/sxa/Datasets/Human_DLPFC/151673/151673_filtered_feature_bc_matrix.h5'
    slice = sc.read_10x_h5(data_dir+str(slice_name)+'/'+str(slice_name)+'_filtered_feature_bc_matrix.h5')
    coldata = pd.read_csv(data_dir+str(slice_name)+'/'+'col_data_'+str(slice_name)+'.csv')
    spatial_coor = (coldata.loc[:, ['row', 'col']]).values
    image_coor = (coldata.loc[:, ['imagerow', 'imagecol']]).values
    spatialLIBD_cluster = (coldata.loc[:,['Cluster']]).values
    real_label = (coldata.loc[:, ['layer_guess_reordered_short']]).values
    slice.obsm['spatial'] = spatial_coor
    slice.obsm['images_coordinates'] = image_coor
    slice.obs['Ground Truth'] = real_label.reshape(-1)
    slice.obs['spatialLIBD_cluster'] = spatialLIBD_cluster.reshape(-1)
    var = list(set(slice.var.index))
    slice.var_names_make_unique()
    slice = slice[:, var]
    sc.pp.filter_genes(slice, min_counts = 5)
    sc.pp.highly_variable_genes(slice, n_top_genes=3000 ,flavor='seurat_v3')
    slice = slice[:, slice.var.highly_variable]
    # 标准化 log (验证标准化的有效性)
    sc.pp.normalize_total(slice, target_sum=1e4)
    # sc.pp.normalize_total(slice)
    sc.pp.log1p(slice)

    return slice

def load_Human_DLPFC_SpaGCN(data_dir='/home/sxa/Datasets/Human_DLPFC/', slice_name=151673):
    slice = sc.read_10x_h5(data_dir+str(slice_name)+'/'+str(slice_name)+'_filtered_feature_bc_matrix.h5')
    coldata = pd.read_csv(data_dir+str(slice_name)+'/'+'col_data_'+str(slice_name)+'.csv')
    spatial_coor = (coldata.loc[:, ['row', 'col']]).values
    real_label = (coldata.loc[:, ['layer_guess_reordered_short']]).values
    slice.obsm['spatial_coor'] = spatial_coor
    slice.obs['Ground Truth'] = real_label.reshape(-1)
    # SpaGCN
    spatial=pd.read_csv(data_dir+str(slice_name)+'/spatial/tissue_positions_list.csv' ,sep=",",header=None,na_filter=False,index_col=0)  # 改地址
    slice.obs["x1"]=spatial[1]
    slice.obs["x2"]=spatial[2]
    slice.obs["x3"]=spatial[3]
    slice.obs["x4"]=spatial[4]
    slice.obs["x5"]=spatial[5]
    slice.obs["x_array"]=slice.obs["x2"]
    slice.obs["y_array"]=slice.obs["x3"]
    slice.obs["x_pixel"]=slice.obs["x4"]
    slice.obs["y_pixel"]=slice.obs["x5"]
    #Select captured samples
    slice=slice[slice.obs["x1"]==1]
    slice.var_names=[i.upper() for i in list(slice.var_names)]
    slice.var["genename"]=slice.var.index.astype("str")

    var = list(set(slice.var.index))
    slice.var_names_make_unique()
    slice = slice[:, var]
    spg.prefilter_genes(slice,min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(slice)
    sc.pp.highly_variable_genes(slice, n_top_genes=3000 ,flavor='seurat_v3')
    slice = slice[:, slice.var.highly_variable]
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(slice)
    sc.pp.log1p(slice)

    return slice

def load_cSCC_46(data_dir='/home/sxa/Datasets/cSCC/', P_name='P4', slice_num='1'):
    slice = sc.read_10x_mtx(data_dir+P_name+'/'+P_name+'_rep'+slice_num)
    position = sc.read_csv(data_dir+P_name+'/'+P_name+'_rep'+slice_num+'/spatial/'+P_name+'_rep'+slice_num+'_tissue_positions_list.csv')
    position = position[slice.obs.index]
    spatial_coor = np.array([row[1:3] for row in position.X.tolist()]).astype(int)
    spatial = np.array([row[3:5] for row in position.X.tolist()]).astype(int)
    slice.obsm['spatial_coor']=spatial_coor
    slice.obsm['spatial'] = spatial
    sc.pp.filter_cells(slice, min_counts = 3)
    sc.pp.highly_variable_genes(slice, n_top_genes=3000 ,flavor='seurat_v3')
    slice = slice[:, slice.var.highly_variable]
    sc.pp.normalize_total(slice, target_sum=1e4)
    sc.pp.log1p(slice)
    return slice

def load_cSCC_25910(data_dir='/home/sxa/Datasets/cSCC/', P_name='P2', slice_num='1'):
    slice = sc.read_text(data_dir+P_name+'/'+P_name+'_ST_rep'+slice_num+'_stdata.tsv')
    position = sc.read_text(data_dir+P_name+'/spot_data-selection-'+P_name+'_ST_rep'+slice_num+'.tsv')
    spatial_coor = np.array([row[:2] for row in position.X.tolist()]).astype(int)
    index = [f'{x}x{y}' for x, y in spatial_coor]
    position.obs.index = index
    index = np.intersect1d(index, slice.obs.index)
    slice = slice[index]
    position = position[index]
    spatial_coor = []
    spatial_coor = np.array([row[:2] for row in position.X.tolist()]).astype(int)
    spatial = np.array([row[4:6] for row in position.X.tolist()]).astype(int)
    slice.obsm['spatial_coor'] = spatial_coor
    slice.obsm['spatial'] = spatial
    sc.pp.filter_cells(slice, min_counts = 3)
    sc.pp.highly_variable_genes(slice, n_top_genes=3000 ,flavor='seurat_v3')
    slice = slice[:, slice.var.highly_variable]
    sc.pp.normalize_total(slice, target_sum=1e4)
    sc.pp.log1p(slice)
    return slice