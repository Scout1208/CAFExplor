import scanpy as sc
import matplotlib.pyplot as plt     
import json
import pandas as pd
# 加載數據
adata = sc.read_10x_mtx(
    'filtered_gene_bc_matrices/hg19',  # 使用解壓的目錄
    var_names='gene_symbols',         # 用基因名稱而非基因編碼
    cache=True
)
print(adata)
# Ensure that all observation (cell) names and variable (gene) names are unique.
# This prevents errors in downstream analyses caused by duplicate names.
adata.obs_names_make_unique()
adata.var_names_make_unique()
adata.raw = adata.copy()
print(adata)

# 篩選細胞與基因
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(adata)
# 基本質控
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # 標記線粒體基因
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# # 視覺化 QC 結果
# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)

#Remove cells that have too many mitochondrial genes expressed or too many total counts:
#You can set your preferred criteria
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
print(adata)
# 標準化與降維
# Save the raw count data to a separate layer called 'counts' for future reference.
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_transformed"] = adata.X.copy()

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(adata)
print(adata.var)

# 進行 PCA 與 UMAP
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, metric='cosine')
sc.tl.umap(adata)

# 聚類
# sc.tl.leiden(adata)
sc.tl.leiden(adata, resolution=0.8)

# 畫出 UMAP
# sc.pl.umap(adata, color=['leiden'], title=["Leiden clustering"], frameon=False, legend_loc='right margin', save=False, size=10, legend_fontsize=12)
cell_type = {"0":'CD4 T',
    "1":'CD14 Monocytes',
    "2":'CD8 T & NK',
    "3":'B cells',
    "4":'FCGR3A+ Monocytes',
    "5":'Dendritic',
    "6":'Megakaryocytes',
    "7":'NK or Unknown'
}
#Manual Annotation for corresponding leiden clusters
adata.obs["Manual Annotation"] = adata.obs.leiden.map(cell_type)

# export.py
# 導出 UMAP 座標與聚類標籤
umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_df['leiden'] = adata.obs['leiden'].values
umap_df.to_json('umap_data.json', orient='records')

# 導出基因表達數據（熱圖數據）
# gene_list = ['CD3D', 'MS4A1', 'GNLY']  # 舉例幾個感興趣的基因
# heatmap_data = adata[:, gene_list].X.toarray()
# 將 heatmap_data.json 改為包含所有基因
heatmap_data = adata.X.toarray()  # 提取所有基因的數據
gene_list = adata.var_names.tolist()  # 獲取所有基因名稱
heatmap_df = pd.DataFrame(heatmap_data, columns=gene_list)
heatmap_df['leiden'] = adata.obs['leiden'].values
heatmap_df.to_json('heatmap_data.json', orient='records')
