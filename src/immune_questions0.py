import numpy as np
import scanpy as sc



adata = sc.read_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked.h5ad") # enter h5ad file name 

### Visualization exercises

print("Prints an overview of information on the AnnData object" "\n",  adata)
#
print("View the count matrix (rows x columns, cells x genes)" "\n", adata.X)
#
print("View a pandas data frame containing metadata on the cells" "\n", adata.obs.head())
#
first5 = adata.obsm['X_umap'][:5]
print("First 5 cells' Umap coordinates" "\n", first5)
#
print("Pandas data frame containing metadata on the genes" "\n", adata.var.head())


### How many and which samples are inside your Anndata? ###
samples = adata.obs['sample_id'].unique()
num_samples = len(samples)

print(f"{num_samples} unique samples")
print(samples)


### How many cells from each category of 'injury' annotation are there? Hint: if we want cell-associated annotations/metadata, we should look into .obs to solve the problem ###
injury = adata.obs['injury'].value_counts()
injury_total = adata.obs['injury'].count()

print(injury)
print(f"Total number of injuries: {injury_total}")



### How many cells from each category of 'day' are there? ###
day = adata.obs['day'].value_counts()
day_total = adata.obs['day'].count()

print(day)
print(f"Total number of days: {day_total}")


### How many cells from each category of 'injury_day' are there? ###
injury_days = adata.obs['injury_day'].value_counts()
print(injury_days)


### What is the mean expression of "Gapdh" in the entire data set? ###
GapdhgeneList = adata[:, 'Gapdh'].X.toarray().mean()
Gapdhgene = adata[:, 'Gapdh'].X.mean()

print(GapdhgeneList,Gapdhgene)


### In which cell is the highest expression of the gene "Gapdh" found? And the lowest? ###
Gapdhgene = adata[:, 'Gapdh'].X

max_index = np.argmax(Gapdhgene)
min_index = np.argmin(Gapdhgene)

print(max_index)
print(min_index)

max_cell = adata.obs_names[max_index]
min_cell = adata.obs_names[min_index]

print(max_cell)
print(min_cell)


### In which cluster (resolution 0.7) is the mean expression of the gene "Gpnmb" the highest? Note: later on you may plot this in a dotplot and confirm these results ###
Gpnmbgene = adata[:, ['Gpnmb','Gapdh']].X.toarray()
print(Gpnmbgene[-4:])

zeroseven = adata.obs['leiden_n15_r0.7']
print(zeroseven)

adata.obs[['Gpnmb_expression','Gapdh_expression']] = Gpnmbgene # add Gpnmbgene in obs
print(adata.obs['Gpnmb_expression'])
print(adata.obs['Gapdh_expression'])

cluster_means = adata.obs.groupby('leiden_n15_r0.7')['Gpnmb_expression'].mean() # pics que 0.7 resolution and group and combain every data from that resolution in the expression values of Gpnmb and returns the mean 

max_cluster = cluster_means.idxmax()
print(f"The cluster with the max expression is {max_cluster}")

min_cluster = cluster_means.idxmin()
print(f"The cluster with the min expression is {min_cluster}")

