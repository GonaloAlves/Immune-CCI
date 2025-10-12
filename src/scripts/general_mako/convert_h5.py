import scanpy as sc

# # for 10x HDF5 file:
# adata = sc.read_10x_h5("/home/makowlg/Documents/Immune-CCI/h5ad_files/test_h5/GSM5723843_Donor1_raw_feature_bc_matrix.h5")   # reads HDF5 from 10x
# adata.var_names_make_unique()
# adata.write_h5ad("/home/makowlg/Documents/Immune-CCI/h5ad_files/test_h5/result/dataset1.h5ad")

# # if you have matrix.mtx + features.tsv + barcodes.tsv:
# adata = sc.read_10x_mtx("folder_with_mtx_and_barcodes/", var_names='gene_symbols', make_unique=True)
# adata.write_h5ad("dataset_from_mtx.h5ad")

# Load the dataset
def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    print("Loading h5ad file")
    return sc.read_h5ad(file_path)




# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/test_h5/result/dataset1.h5ad")

    print(adata.var)
    print("#############")
    print(adata.obs.head())
    print("#############")
    print(adata.X)