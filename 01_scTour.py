import sctour as sct
import omicverse as ov
#print(f"omicverse version: {ov.__version__}")
#ov.utils.ov_plot_set()
import os
import numpy as np
import scanpy as sc
#print(f"scanpy version: {sc.__version__}")
import anndata as ad
import pandas as pd


os.getcwd()  ##查看当前路径


os.chdir('/home/luchun/scRNA/Neutrophil/01_Human_BM/scTour/')

sc.settings.set_figure_params(dpi=300, facecolor="white")


adata_all = sc.read_h5ad("../05_anno.h5ad")
adata_all
#AnnData object with n_obs × n_vars = 21861 × 2000

adata_all.uns['Neu_anno_colors'] = ["#F46D43","#E6F598","#66C2A5","#5E4FA2"]



#——————————————随机抽取1w个细胞————————————————————————
# 随机抽取 10000 个细胞
n_cells = 10000
random_indices = np.random.choice(adata_all.n_obs, n_cells, replace=False)

# 创建一个新的 AnnData 对象
adata = adata_all[random_indices, :].copy()
adata
#AnnData object with n_obs × n_vars = 10000 × 2000

sc.pl.umap(adata, color='Neu_anno',  show=False,save="_test.png")







adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 10000 × 17610
print(adata_counts.X)




ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 10000 × 17610
print(adata_counts.X)


dense_array = adata_counts.X.toarray()
float_array = dense_array.astype(np.float32)
print(float_array.dtype)

adata_counts.X = float_array


print(adata_counts.X)




adata=adata_counts

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)





#——————————————————————Train the scTour model————————————————————————
tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
tnode.train()



#——————————————————————Infer cellular dynamics——————————————————————
adata.obs['ptime'] = tnode.get_time()

sc.pl.umap(adata, color='ptime',  show=False,save="_01-1.png")


mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs


adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])


#adata = adata[np.argsort(adata.obs['ptime'].values), :]


adata.uns['leiden_res_1_colors'] = ["#ffe788","#e8743c","#6a51a3","#CE9FCA","#9C9BE9","#76afda","#9e9ac8","#2873B3","#b0d45d"]

sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='leiden_res_1',frameon=False, show=False, legend_loc='none', size=18, alpha=1,save="./figures/stream_01-2.png")


sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='Neu_anno',frameon=False, show=False, legend_loc='none', size=18, alpha=1,save="./figures/stream_01-3.png")




##reverse pseudotime
adata.obs['ptime'] = sct.train.reverse_time(adata.obs['ptime'].values)

sct.vf.plot_vector_field(adata, reverse=True, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='leiden_res_1', show=False, frameon=False, size=18, alpha=1,save="./figures/stream_01-4.png")

sct.vf.plot_vector_field(adata, reverse=True, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='Neu_anno', show=False, frameon=False, size=18, alpha=1,save="./figures/stream_01-5.png")



adata.write_h5ad('Human_BM_scTour.h5ad',compression='gzip')

#adata=sc.read_h5ad("Human_BM_scTour.h5ad")
#adata
