import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked.h5ad") # enter h5ad file name 

## DGEs ##

sc.tl.rank_genes_groups(adata, groupby='leiden_fusion', method='wilcoxon', use_raw=False, pts=True)

dge_results = adata.uns['rank_genes_groups']


# Extracting values to DataFrames
gene_names = pd.DataFrame(dge_results['names'])
logfoldchanges = pd.DataFrame(dge_results['logfoldchanges'])
pvals_adj = pd.DataFrame(dge_results['pvals_adj'])
scores = pd.DataFrame(dge_results['scores'])
pts = pd.DataFrame(dge_results['pts'])

cluster_dfs = {}

for i in range(len(gene_names.columns)): #for loop to take the nr of columns
    
    gene_reindex = gene_names.iloc[:, i] # take the correct order of genes names
    
    pts_reindexed = pts.iloc[:, i].reindex(gene_reindex.values) # Reindex the pts
    
    cluster_df = pd.DataFrame({ # creates a df that takes all the data corresponding to the cluster in loop
        'score': scores.iloc[:, i].values,
        'logfoldchange': logfoldchanges.iloc[:, i].values,
        'pvals_adj': pvals_adj.iloc[:, i].values,
        'pts': pts_reindexed.values}, 
        index=gene_reindex.values)
    # cluster_df.index.name = 'gene'
    cluster_dfs[f'Cluster_{i}'] = cluster_df    # append it to the dictionary

#Print
for cluster, df in cluster_dfs.items():
    print(f"\n{cluster}:")
    print(df.head())


#Filter, sort, and select top 5 genes per cluster
top_genes_cluster = {}

for cluster, df in cluster_dfs.items():
    filtered_df = df.mask(df['pvals_adj'] >= 0.05) # padj < 0.05
    sorted_df = filtered_df.sort_values(by='logfoldchange', ascending=False) # orders by fold change
    top_genes = sorted_df.head(5) # top 5 genes
    
    top_genes_cluster[cluster] = top_genes # Assigns the DataFrame (top_genes) to the key cluster in the dictionary top_genes_per_cluster.

# # print 
# for cluster, df in top_genes_per_cluster.items():
#     print(f"\nTop genes in {cluster}:")
#     print(df)


# append the genes to a full list in cluster order, every 5 genes is a cluster
top_genes = []

for cluster, df in top_genes_cluster.items():
    for gene in df.index:
        top_genes.append(gene)


print("\nTop genes per cluster:")
print(top_genes)

# Complete the code to save and load adata_deg.
save_path = '/home/makowlg/Documents/Immune-CCI/adata_deg.h5ad'
adata.write_h5ad(save_path, compression='gzip')

adata_deg = sc.read_h5ad(save_path)
print(adata_deg)    # To confirm the file has the DGE data stored inside
