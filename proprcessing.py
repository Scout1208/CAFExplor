import scanpy as sc
import json
import pandas as pd

# 加載數據
adata = sc.read_10x_mtx(
    'filtered_gene_bc_matrices/hg19',  # 使用解壓的目錄
    var_names='gene_symbols',         # 用基因名稱而非基因編碼
    cache=True
)
print(adata)

# 確保觀察（細胞）名稱和變量（基因）名稱唯一
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

# 移除表達過多線粒體基因或總計數過高的細胞
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
print(adata)

# 標準化與降維
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_transformed"] = adata.X.copy()

# 計算高變量基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(adata)
print(adata.var)

# 選取前 100 個高變量基因
top_genes = adata.var.highly_variable.sort_values(ascending=False).head(50).index.tolist()
print(f"選取的前50個高變量基因: {top_genes}")

# 進行 PCA 與 UMAP
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, metric='cosine')
sc.tl.umap(adata)

# 聚類
sc.tl.leiden(adata, resolution=0.8)

# 定義 Cell Type 名稱
cell_type = {
    "0": 'CD4 T',
    "1": 'CD14 Monocytes',
    "2": 'CD8 T & NK',
    "3": 'B cells',
    "4": 'FCGR3A+ Monocytes',
    "5": 'Dendritic',
    "6": 'Megakaryocytes',
    "7": 'NK or Unknown'
}

# 手動註釋對應的 leiden 聚類
adata.obs["Manual Annotation"] = adata.obs.leiden.map(cell_type)

# 導出 UMAP 座標與聚類標籤
umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_df['leiden'] = adata.obs['leiden'].values
umap_df.to_json('umap_data.json', orient='records')

# 導出基因表達數據（熱圖數據）
# 計算每個leiden聚類的基因表達平均值
heatmap_df = adata[:, top_genes].to_df().copy()
heatmap_df['leiden'] = adata.obs['leiden'].values

# 計算每個聚類的基因平均表達
heatmap_cluster_df = heatmap_df.groupby('leiden').mean().reset_index()

# 將聚類標籤轉換為細胞類型名稱
heatmap_cluster_df['cell_type'] = heatmap_cluster_df['leiden'].map(cell_type)

# 將 'leiden' 列移除或保留根據需求
heatmap_cluster_df = heatmap_cluster_df.drop(columns=['leiden'])

# 將 DataFrame 轉換為適合 D3.js 的格式
# 我們將轉換為一個包含基因、細胞類型及表達值的陣列
heatmap_data = []
for _, row in heatmap_cluster_df.iterrows():
    cell_type_name = row['cell_type']
    for gene in top_genes:
        heatmap_data.append({
            'gene': gene,
            'cell_type': cell_type_name,
            'expression': row[gene]
        })

# 導出為 JSON 文件
with open('heatmap_data.json', 'w') as f:
    json.dump(heatmap_data, f, indent=4)

print("UMAP 和 Heatmap 數據已成功導出。")
