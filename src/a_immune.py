import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import os
import shutil

adata = sc.read_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy.h5ad") # enter h5ad file name 
#print(adata)

## DGEs ##

#sc.tl.rank_genes_groups(adata, groupby='leiden_fusion', method='wilcoxon', use_raw=False, pts=True)

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
    #print(f"Gene e\n{pts_reindexed}:")

    cluster_df = pd.DataFrame({ # creates a df that takes all the data corresponding to the cluster in loop
        'score': scores.iloc[:, i].values,
        'logfoldchange': logfoldchanges.iloc[:, i].values,
        'pvals_adj': pvals_adj.iloc[:, i].values,
        'pts': pts_reindexed.values}, 
        index=gene_reindex.values)
    one_sorted_df = cluster_df.sort_values(by='logfoldchange', ascending=False)
    # cluster_df.index.name = 'gene'

    cluster_names = gene_names.columns[i]
    cluster_dfs[f'Cluster:{cluster_names}'] = one_sorted_df    # append it to the dictionary

# #Print
# for cluster, df in cluster_dfs.items():
#     print(f"\n{cluster}:")
#     print(df.head())

# Filter, sort, and select top 5 genes per cluster
top_genes_cluster = {}

for cluster, df in cluster_dfs.items():
    mask = df['pvals_adj'] < 0.05
    filtered_df = df[mask]
    #filtered_df = df.mask(df['pvals_adj'] >= 0.05) # padj < 0.05
    #sorted_df = filtered_df.sort_values(by='logfoldchange', ascending=False) # orders by fold change
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

#ordenar primeiro por logfoldchange por ordem decrescente, ver e dps quando for retirar os pvalue mais elevados ver quandos genes sobram que tenham um
#logfoldchange inferior a 0, dependendo o numero vamos buscar os top n de genes com o logfoldchange superior da lista inicial de genes ordenados

# print 
for cluster, df in top_genes_cluster.items():
    print(f"\nTop genes in {cluster}:")
    print(df)


# # append the genes to a full list in cluster order, every 5 genes is a cluster
# top_genes = []

# for cluster, df in top_genes_cluster.items():
#     for gene in df.index:
#         top_genes.append(gene)


# print("\nTop genes per cluster:")
# print(top_genes)



# ### Visualize the DotPlots of the DGE's ###

# # New directory 
# dotplots = "dotplots"

# # Delete the existing dir if exists
# if os.path.exists(dotplots):
#     shutil.rmtree(dotplots)

# os.makedirs(dotplots)

## Creates the dir
#if not os.path.exists(dotplots):
#   os.makedirs(dotplots)


# # DotPlot 1 representing the ranked genes by deafault method
# dotplot1 = sc.pl.rank_genes_groups_dotplot(
#     adata,
#     n_genes=5,
#     groupby='leiden_fusion',
#     key = 'rank_genes_groups_leiden_fusion',
#     cmap='bwr',
#     vmin=-4,
#     vmax=4,
#     values_to_plot = 'logfoldchanges',
#     #var_names = ['Arhgap15','Parp14','Diaph3'],
#     #min_logfoldchange=1,
#     colorbar_title='log fold change',
#     use_raw=False,
#     dendrogram=False,
#     return_fig=True
# )


# # DotPlot 2 representing the ranked genes by top 5 filter
# dotplot2 = sc.pl.rank_genes_groups_dotplot(
#     adata,
#     var_names=top_genes,
#     groupby='leiden_fusion',
#     key = 'rank_genes_groups_leiden_fusion',
#     cmap='bwr',
#     vmin=-4,
#     vmax=4,
#     values_to_plot='logfoldchanges',
#     #min_logfoldchange=1,
#     colorbar_title='log fold change',
#     use_raw=False,
#     dendrogram=False,
#     return_fig=True
# )


# output_path1 = os.path.join(dotplots, "dotplot_1.png")
# dotplot1.savefig(output_path1, bbox_inches="tight")
# plt.close()



# output_path2 = os.path.join(dotplots, "dotplot_2.png")
# dotplot2.savefig(output_path2, bbox_inches="tight")
# plt.close()  # Close the current figure to avoid overlap


print("hello world")