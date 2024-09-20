import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import os
import shutil

adata = sc.read_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy.h5ad") # enter h5ad file name 
#print(adata)

## DGEs ##

dge_fusion = adata.uns['rank_genes_groups_leiden_fusion']
# print(dge_fusion)


# Extracting values to DataFrames
gene_names = pd.DataFrame(dge_fusion['names'])
logfoldchanges = pd.DataFrame(dge_fusion['logfoldchanges'])
pvals_adj = pd.DataFrame(dge_fusion['pvals_adj'])
scores = pd.DataFrame(dge_fusion['scores'])
pts = pd.DataFrame(dge_fusion['pts'])


# print("Cluster Names:")
# for cluster in gene_names.columns:
#     print(f"Cluster {cluster}")

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
    
    mask_pts = (cluster_df['pts'] >= 0.5) 
    filtered_df = cluster_df[mask_pts]
    
    one_sorted_df = filtered_df.sort_values(by='logfoldchange', ascending=False)


    cluster_names = gene_names.columns[i]
    cluster_dfs[cluster_names] = one_sorted_df    # append it to the dictionary

# Delete the NA cluster

clusters_to_delete = [cluster for cluster in cluster_dfs if cluster.endswith("NA")]

for cluster in clusters_to_delete:
    del cluster_dfs[cluster]

# #Print
# for cluster, df in cluster_dfs.items():
#     print(f"\n{cluster}:")
#     print(df.head())

# Filter, sort, and select top 5 genes per cluster
top_genes_cluster = {}

for cluster, df in cluster_dfs.items():
    mask_pval = df['pvals_adj'] < 0.05
    filtered_df = df[mask_pval]
    top_genes = filtered_df.head(5)
    
    if (top_genes['logfoldchange'] < 0).any(): # top 5 genes

        num_nega_lfc = (top_genes['logfoldchange'] < 0).sum()

        pos_lfc = one_sorted_df.head(num_nega_lfc)

        mask1 = top_genes['logfoldchange'] > 0
        top_genes_positive = top_genes[mask1]

        replacements = []
        for gene, values in df.iterrows():
            if gene not in top_genes_positive.index:
                replacements.append(values)

            if len(replacements) == num_nega_lfc:
                break  # Stop once we have enough replacement genes

        replacements_df = pd.DataFrame(replacements)

        combined_top_genes = pd.concat([top_genes_positive,replacements_df])

        
        top_genes_cluster[cluster] = combined_top_genes

    elif top_genes.empty: 
        top_genes_new = df.head(5)

        top_genes_cluster[cluster] = top_genes_new

    else:
        top_genes_cluster[cluster] = top_genes # Assigns the DataFrame (top_genes) to the key cluster in the dictionary top_genes_per_cluster.


# print 
for cluster, df in top_genes_cluster.items():
    print(f"\nTop genes in {cluster}:")
    print(df)


# Create a dictionary to hold top gene names for each cluster
top_genes_names = {}

# Iterate over the top genes per cluster
for cluster, df in top_genes_cluster.items():
    # Collect the top 5 gene names for this cluster
    top_genes_names[cluster] = df.index.tolist()


# print("\nTop genes per cluster:")
# print(top_genes_names)


## Visualize the DotPlots of the DGE's ###

# New directory 
dotplots = "dotplots"

# Delete the existing dir if exists
if os.path.exists(dotplots):
    shutil.rmtree(dotplots)

os.makedirs(dotplots)

# Creates the dir
if not os.path.exists(dotplots):
  os.makedirs(dotplots)


# DotPlot representing the ranked genes by top 5 filter
dotplot1 = sc.pl.rank_genes_groups_dotplot(
    adata,
    var_names = top_genes_names,
    groupby='leiden_fusion',
    key = 'rank_genes_groups_leiden_fusion',
    cmap='bwr',
    vmin=-4,
    vmax=4,
    values_to_plot='logfoldchanges',
    colorbar_title='log fold change',
    use_raw=False,
    dendrogram=False,
    return_fig=True
)

output_path1 = os.path.join(dotplots, "dotplot_1.png")
dotplot1.savefig(output_path1, bbox_inches="tight")
# plt.close()  # Close the current figure to avoid overlap
plt.close()

# lista de labels do eixo do x e quando tiver essa lista substituir por uma nova - not done
# apagar o NA - done
# min de 30%, 40%, 50%
