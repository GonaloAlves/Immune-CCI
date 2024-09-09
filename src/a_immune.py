import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

print(np.__version__)

adata = sc.read_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked.h5ad") # enter h5ad file name here
print(adata)