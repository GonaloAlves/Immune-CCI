# Decoding the immune crosstalk in spinal cord injury with snRNA-SEQ


This repository contains the code and data for my Master's thesis project, focused on analyzing gene expression in meningeal vascular cells. The project uses single-cell RNA sequencing (scRNA-seq) data to study differential gene expression, cluster annotations, and canonical gene markers.

## 📌 Features
- **Loading & Processing scRNA-seq Data**: Reads `.h5ad` files and filters unwanted clusters.
- **Clustering & Dendrograms**: Performs Leiden clustering and hierarchical clustering for visualization.
- **Gene Expression Analysis**: Extracts differentially expressed genes (DEGs) and filters by significance.
- **Custom Cluster Merging & Recovery**: Merges and unmarges clusters while preserving names.
- **Dotplot Visualizations**: Generates dotplots for gene expression across clusters.
- **Excel Export**: Exports gene expression results into structured Excel files.
- **Gene-Based Filtering**: Allows filtering of cells based on specific gene expression.

## 🛠 Installation
### Requirements
- Python 3.9+
- Conda or virtualenv recommended
- Required Python packages:
  ```bash
  pip install scanpy pandas matplotlib seaborn numpy openpyxl
  ```

### Cloning the Repository
```bash
git clone https://github.com/yourusername/meningeal-vascular-analysis.git
cd meningeal-vascular-analysis
```

## 🚀 Usage
### 1. Load Data
Modify `input_file` to point to your `.h5ad` dataset.
```python
adata = load_data("/path/to/your/adata.h5ad")
```

### 2. Remove Unwanted Clusters
```python
clusters_to_remove = ['MeV.Immune_doublets.0', 'MeV.Low_Quality.0']
adata_filtered = remove_clusters(adata, clusters_to_remove)
```

### 3. Generate Dotplots
```python
custom_cluster_order = ["MeV.Endothelial.0", "MeV.Pericytes.0"]
create_dotplots_with_thresholds(adata_filtered, genes, [0, 0.2, 0.3], custom_cluster_order)
```

### 4. Export Results to Excel
```python
export_to_excel_per_gene_group(filtered_genes, cluster_dfs)
```

### 5. Filter by Gene Expression
```python
adata_filtered = filter_cells_by_gene_expression(adata, "GeneX")
```

## 📊 Expected Outputs
- **Dotplots**: Saved in `canonical/canonical_dalila/`
- **Filtered Expression Data**: Exported to `excels/meningeal/updates/`
- **Cluster Metadata**: Updated in `.h5ad` files

## 📄 Project Structure
```
├── data/  # Contains raw .h5ad files
├── scripts/  # Python scripts for analysis
├── results/  # Outputs like dotplots and Excel files
├── README.md  # This documentation
```

## 🤝 Contributing
Feel free to submit pull requests for improvements and bug fixes.

## 📜 License
This project is open-source and licensed under the MIT License.

---
✉️ **Contact**: [Your Name] - your.email@example.com

📌 **GitHub Repository**: [GitHub Link]

