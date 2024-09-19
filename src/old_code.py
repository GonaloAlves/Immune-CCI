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

