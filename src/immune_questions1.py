import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import os
import shutil

adata = sc.read_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked.h5ad") # enter h5ad file name 


### Visualize clusters ###

### How many resolutions were defined? ###

different_reso = []

for reso in adata.obs.columns: 
    if reso.startswith('leiden_'):
        different_reso.append(reso)

print(len(different_reso))


### Visualize the UMAP colored by each resolution ###

# New directory 
reso_clusters = "reso_clusters"

# Delete the existing dir
if os.path.exists(reso_clusters):
    shutil.rmtree(reso_clusters)

# Creates the dir
if not os.path.exists(reso_clusters):
   os.makedirs(reso_clusters)

# Loop through the different resolutions and save the UMAP for each one
for resolution in different_reso:
    # Create the UMAP plot and return the figure
    ax = sc.pl.umap(adata, color=resolution, title=f"UMAP - {resolution}", return_fig= True)
    
    output_path = os.path.join(reso_clusters, f"umap_{resolution}.png")
    ax.savefig(output_path, bbox_inches = "tight") #Save
    
    ax.clf() #Close fig

