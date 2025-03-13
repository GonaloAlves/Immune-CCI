import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import os
import shutil

adata = sc.read_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy.h5ad") # enter h5ad file name 
#print(adata)

## DGEs ##

dge_fusion = adata.uns['rank_genes_groups_leiden_fusion']
# print(dge_fusion)s


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



#Unsed code from immune analysis

# DGEs #

#sc.tl.rank_genes_groups(adata, groupby='leiden_fusion', method='wilcoxon', use_raw=False, pts=True)

# ####################################################

#filter the top 5 genes:

#filtered_df = df.mask(df['pvals_adj'] >= 0.05) # padj < 0.05 
#sorted_df = filtered_df.sort_values(by='logfoldchange', ascending=False) # orders by fold change

# ####################################################

# Dotplot creation (simple data, selected by wilcoxon scores)

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

# ####################################################

# Dani explanation of comprehension of for loops

# d1 = {'a': 1, 'b': 2, 'c': 3}
# print(d1)
# d1 = {f'{k}*': v for k, v in d1.items()}
# print(d1)

# ####################################################

# Try to do a copy of the adata in order to add the * on the genes that are statiscaly unsignificant 

# adata.var['modified_gene_names'] = adata.var_names.copy()

# # Iterate over the clusters and gene data in top_genes_cluster
# for cluster, df in top_genes_cluster.items():
#     # Iterate over the genes and their data
#     for gene, row in df.iterrows():
#         # Check if p-value is not statistically significant (adj_pval >= 0.05)
#         if row['pvals_adj'] >= 0.05:
#             # Update the gene name in the 'modified_gene_names' column by appending '*'
#             adata.var.loc[gene, 'modified_gene_names'] = gene + '*'

# # Now modify adata.var_names to be the modified names
# adata.var_names = adata.var['modified_gene_names']

# # You can now directly plot with the updated names in adata
# dotplot1 = sc.pl.dotplot(
#     adata,
#     var_names=top_genes_cluster,   # The dictionary of top genes by cluster
#     groupby='leiden_fusion',       # Cluster grouping
#     cmap='bwr',                    # Colormap for logfoldchanges
#     vmin=-4,                       # Minimum value for colormap
#     vmax=4,                        # Maximum value for colormap
#     values_to_plot='logfoldchanges',  # Plot logfoldchanges
#     colorbar_title='log fold change',  # Title for color bar
#     use_raw=False,                 # Use the raw data or not
#     dendrogram=False,              # Do not display dendrogram
#     return_fig=True                # Return the plot figure object
# )

# ####################################################

# Same objective of the last code but now we are just selecting the genes and giving the * on a diferent list (not working)

# adjusted_top_genes = {}

# for cluster, df in top_genes_cluster.items():
#     # Initialize an empty list to store adjusted gene names for each cluster
#     adjusted_genes = []
    
#     for gene, row in df.iterrows():
#         # Check if p-value is not statistically significant (adj_pval >= 0.05)
#         if row['pvals_adj'] >= 0.05:
#             # Append a '*' to non-significant gene names
#             adjusted_genes.append(gene + '*')
#         else:
#             # Keep the gene name unchanged for significant genes
#             adjusted_genes.append(gene)
    
#     # Store the adjusted gene names for each cluster in the dictionary
#     adjusted_top_genes[cluster] = adjusted_genes

# dotplot1 = sc.pl.dotplot(
#     adata,
#     var_names = adjusted_top_genes,   
#     groupby='leiden_fusion',       
#     cmap='bwr',                    
#     vmin=-4,                       
#     vmax=4,                        
#     values_to_plot='logfoldchanges',  
#     colorbar_title='log fold change',  
#     use_raw=False,                 
#     dendrogram=False,              
#     return_fig=True                
# )


#Deleting NA cluster 
# remove_cluster = []

# for cluster in cluster_dfs:
#     if cluster.endswith("NA"):
#         remove_cluster.append(cluster)

# for cluster in remove_cluster:
#     del cluster_dfs[cluster]
