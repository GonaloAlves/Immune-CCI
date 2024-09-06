# Helper functions
#
# General function that may be used troughout the analysis scripts

# Daniel Ribeiro, 2023

from typing import Literal, Union

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse import issparse, sparray


def plot_umap(data: sc.AnnData,
              dest_dir: str,
              color: Union[str, list[str]],
              name_suffix: Union[str, None] = None,
              neighbors_key: Union[str, None] = None,
              multipanel: bool = False,
              multipanel_main_title: Union[str, None] = None,
              palette=None,
              **kwargs) -> None:
    ax = sc.pl.umap(data,
                    use_raw=False,
                    color=color,
                    palette=palette,
                    wspace=0.3,
                    ncols=3,
                    neighbors_key=neighbors_key,
                    show=False,
                    **kwargs)
    if multipanel:
        plt.suptitle(f"UMAP - {multipanel_main_title} - {name_suffix}")
        plt.savefig(dest_dir + f"/UMAP - {multipanel_main_title} - {name_suffix}.pdf", bbox_inches="tight", dpi=600)
    else:
        # Dynamic detemrine legend columns
        if not isinstance(ax, list):
            numitems = len(list(ax.get_legend_handles_labels()[0]))
            nrows = 10
            ncols = int(np.ceil(numitems / float(nrows)))
            if ncols == 0:
                ncols = 1
            ax.legend(loc='center left', ncol=ncols, bbox_to_anchor=(1.05, 0.5), frameon=False)
        #if 'projection' in kwargs.keys():
        #    if kwargs['projection'] == '3d':
        #        print('3d lines')
        #        ax.grid(True)
        plt.title(f"UMAP - {color} - {name_suffix}")
        plt.savefig(dest_dir + f"/UMAP - {color} - {name_suffix}.pdf", bbox_inches="tight", dpi=600)
    plt.close()


def plot_tsne(data: sc.AnnData,
              dest_dir: str,
              color: Union[str, list[str]],
              name_suffix: Union[str, None] = None,
              neighbors_key: Union[str, None] = None,
              multipanel: bool = False,
              multipanel_main_title: Union[str, None] = None,
              palette=None,
              **kwargs) -> None:
    ax = sc.pl.tsne(data,
                    use_raw=False,
                    color=color,
                    palette=palette,
                    alpha=0.5,
                    wspace=0.3,
                    ncols=3,
                    neighbors_key=neighbors_key,
                    show=False,
                    **kwargs)
    if multipanel:
        if name_suffix is None:
            plt.suptitle(f"tSNE - {multipanel_main_title}")
            plt.savefig(dest_dir + f"/tSNE - {multipanel_main_title}.pdf", bbox_inches="tight")
        else:
            plt.suptitle(f"tSNE - {multipanel_main_title} - {name_suffix}")
            plt.savefig(dest_dir + f"/tSNE - {multipanel_main_title} - {name_suffix}.pdf", bbox_inches="tight")
    else:
        # Dynamic determine legend columns
        numitems = len(list(ax.get_legend_handles_labels()[0]))
        nrows = 10
        ncols = int(np.ceil(numitems / float(nrows)))
        if ncols == 0:
            ncols = 1
        ax.legend(loc='center left', ncol=ncols, bbox_to_anchor=(1.05, 0.5), frameon=False)
        if name_suffix is None:
            plt.title(f"tSNE - {color}")
            plt.savefig(dest_dir + f"/tSNE - {color}.pdf", bbox_inches="tight")
        else:
            plt.title(f"tSNE - {color} - {name_suffix}")
            plt.savefig(dest_dir + f"/tSNE - {color} - {name_suffix}.pdf", bbox_inches="tight")
    plt.close()


def _plot_UMAP(adata, set_name, n_neighbors, resolution, clustering_dir, color,
               legend_where: Literal['on data', 'out data', 'both'] = 'out data',
               **kwargs) -> None:
    from src.helper_functions import random_colors
    
    neigh = n_neighbors
    plot_umap(adata, dest_dir=clustering_dir, color='sample_id', name_suffix=f'{set_name}', neighbors_key=f'neighbors_{neigh}', **kwargs)
    plot_umap(adata, dest_dir=clustering_dir, color=color, name_suffix=f'{set_name}', neighbors_key=f'neighbors_{neigh}',
              multipanel=True, multipanel_main_title='per sample', **kwargs)
    for r in resolution:
        ## PAGA
        #sc.tl.paga(adata, groups=f'leiden_n{i}_r{r}', neighbors_key=f'neighbors_{i}')
        #sc.pl.paga(adata, title=f'PAGA leiden_n{i}_r{r}', show=False)
        #plt.savefig(clustering_dir + '/' + f"PAGA - leiden_n{i}_r{r}.pdf", bbox_inches='tight')
        #plt.close()
        # UMAP
        print(f"Plot UMAP for neigh={neigh} res={r}...")
        colors = random_colors(len(adata.obs[f'leiden_n{neigh}_r{r}'].cat.categories.to_list()))
        if legend_where == 'on data' or legend_where == 'both':
            plot_umap(adata,
                      dest_dir=clustering_dir,
                      color=f'leiden_n{neigh}_r{r}',
                      name_suffix=f'{set_name}_cluster',
                      neighbors_key=f'neighbors_{neigh}',
                      legend_loc='on data',
                      legend_fontsize=5,
                      legend_fontweight='bold',
                      palette=colors,
                      **kwargs)
        if legend_where == 'out data' or legend_where == 'both':
            plot_umap(adata,
                      dest_dir=clustering_dir,
                      color=f'leiden_n{neigh}_r{r}',
                      name_suffix=f'{set_name}_cluster',
                      neighbors_key=f'neighbors_{neigh}',
                      palette=colors,
                      **kwargs)


def calculate_UMAP(adata, set_name, n_neighbors, resolution, clustering_dir, color,
                   analyze_data: bool = True,   # whether to also run analyzing functions
                   plot_data: bool = True,      # whether to also run plotting functions
                   dimensions=None,
                   legend_where: Literal['on data', 'out data', 'both'] = 'out data',
                   **kwargs) -> None:
    if analyze_data:
        print("UMAP over Leiden clustering...")
        # Compute PAGA to use in UMAP embedding
        for i in n_neighbors:
            print(f"Computing UMAP for neigh={i}...")
            sc.tl.umap(adata, min_dist=0.5, spread=1.0, alpha=1, neighbors_key=f'neighbors_{i}', init_pos='spectral', random_state=0, **kwargs)
            adata.obsm[f'X_umap_neighbors_n{i}'] = adata.obsm['X_umap'].copy()
        kwargs.pop('n_components', None)
        
        # Plot each sample individualy in the UMAP
        print("Plot each sample in UMAP, as quality control...")
        with PdfPages(clustering_dir + f"/UMAP {set_name} - each sample_id.pdf") as pdf:
            for i in n_neighbors:
                min_xy = np.min(adata.obsm[f'X_umap_neighbors_n{i}'], axis=0)  # min x, y
                max_xy = np.max(adata.obsm[f'X_umap_neighbors_n{i}'], axis=0)  # max x, y
                for sample in color:
                    # Initiate figure and plot
                    plt.figure(figsize=(7, 7), dpi=300)
                    ax = sc.pl.umap(adata[adata.obs[sample] == 1],
                                    use_raw=False, color=None, na_color='black', alpha=0.5, wspace=0.3, ncols=3, show=False, **kwargs)
                    ax.set_xlim(min_xy[0], max_xy[0])
                    ax.set_ylim(min_xy[1], max_xy[1])
                    ax.set_title(f"UMAP {set_name} - sample_id - {sample}")
                    plt.title(sample)
                    pdf.savefig(bbox_inches='tight')
                    plt.close()

    if plot_data:
        if isinstance(dimensions, list):
            for dims in dimensions:
                for i in n_neighbors:
                    _plot_UMAP(adata, set_name, i, resolution, clustering_dir, color, legend_where, dimensions=dims, **kwargs)
        else:
            for i in n_neighbors:
                _plot_UMAP(adata, set_name, i, resolution, clustering_dir, color, legend_where, dimensions=dimensions, **kwargs)


def _plot_tSNE(adata, set_name, n_neighbors, resolution, clustering_dir,
               legend_where: Literal['on data', 'out data', 'both'] = 'out data',
               **kwargs) -> None:
    from src.helper_functions import random_colors
    
    neigh = n_neighbors
    for r in resolution:
        print(f"Plot tSNE for neigh={neigh} res={r}...")
        colors = random_colors(len(adata.obs[f'leiden_n{neigh}_r{r}'].cat.categories.to_list()))
        if legend_where == 'on data' or legend_where == 'both':
            plot_tsne(adata,
                      dest_dir=clustering_dir,
                      color=f'leiden_n{neigh}_r{r}',
                      name_suffix=f'{set_name}_cluster',
                      legend_loc='on data',
                      legend_fontsize=5,
                      legend_fontweight='bold',
                      palette=colors,
                      **kwargs)
        if legend_where == 'out data' or legend_where == 'both':
            plot_tsne(adata,
                      dest_dir=clustering_dir,
                      color=f'leiden_n{neigh}_r{r}',
                      name_suffix=f'{set_name}_cluster',
                      palette=colors,
                      **kwargs)


def calculate_tSNE(adata, set_name, n_neighbors, resolution, rep, clustering_dir, color,
                   analyze_data: bool = True,   # whether to also run analyzing functions
                   plot_data: bool = True,      # whether to also run plotting functions
                   legend_where: Literal['on data', 'out data', 'both'] = 'out data',
                   **kwargs) -> None:
    if analyze_data:
        print("Compute t-SNE...")
        sc.tl.tsne(adata, use_rep=rep, metric='euclidean', early_exaggeration=12, learning_rate=200, perplexity=30, random_state=0)
    
        # Plot each sample individualy in the tSNE
        print("Plot each sample in tSNE, as quality control...")
        with PdfPages(clustering_dir + f"/tSNE {set_name} - each sample_id.pdf") as pdf:
            for i in n_neighbors:
                min_xy = np.min(adata.obsm['X_tsne'], axis=0)  # min x, y
                max_xy = np.max(adata.obsm['X_tsne'], axis=0)  # max x, y
                for sample in color:
                    # Initiate figure and plot
                    plt.figure(figsize=(7, 7), dpi=300)
                    ax = sc.pl.tsne(adata[adata.obs[sample] == 1],
                                    use_raw=False, color=None, na_color='black', alpha=0.5, wspace=0.3, ncols=3, show=False, **kwargs)
                    ax.set_xlim(min_xy[0], max_xy[0])
                    ax.set_ylim(min_xy[1], max_xy[1])
                    ax.set_title(f"tSNE - sample_id - {sample}")
                    plt.title(sample)
                    pdf.savefig(bbox_inches='tight')
                    plt.close()
    if plot_data:
        for i in n_neighbors:
            _plot_tSNE(adata, set_name, i, resolution, clustering_dir, legend_where, **kwargs)
            

def _plot_dendrogram(adata, set_name, n_neighbors, resolution, clustering_dir) -> None:
    
    neigh = n_neighbors
    res = resolution
    print(f"Plot dendrogram for neigh={neigh} res={res}...")
    old_fig_size = plt.rcParams['figure.figsize']
    fig_size = (len(adata.obs[f'leiden_n{neigh}_r{res}'].cat.categories.tolist()) + 1, 7)
    
    sc.settings.set_figure_params(figsize=fig_size)
    sc.pl.dendrogram(adata, groupby=f'leiden_n{neigh}_r{res}', dendrogram_key=f'dendrogram_leiden_n{neigh}_r{res}', show=False)
    plt.title(f"Dendrogram {set_name} - leiden_n{neigh}_r{res}")
    plt.savefig(clustering_dir + '/' + f"Dendrogram {set_name} - leiden_n{neigh}_r{res}.pdf", bbox_inches="tight")
    plt.close()
    
    plt.rcParams['figure.figsize'] = old_fig_size
    
    
def calculate_dendrogram(adata, set_name, n_neighbors, resolution, rep, clustering_dir,
                         analyze_data: bool = True,   # whether to also run analyzing functions
                         plot_data: bool = True,      # whether to also run plotting functions
                         ) -> None:
    if analyze_data:
        print("Compute dendrogram...")
        for i in n_neighbors:
            for r in resolution:
                # Dendrogram cannot be computed if there's only 1 cluster... Check for that
                if len(adata.obs[f'leiden_n{i}_r{r}'].cat.categories.tolist()) > 1:
                    sc.tl.dendrogram(adata,
                                     use_raw=False,
                                     groupby=f'leiden_n{i}_r{r}',
                                     use_rep=rep,
                                     linkage_method='ward',
                                     cor_method='spearman')
    if plot_data:
        for i in n_neighbors:
            for r in resolution:
                if len(adata.obs[f'leiden_n{i}_r{r}'].cat.categories.tolist()) > 1:
                    _plot_dendrogram(adata, set_name, i, r, clustering_dir)


def draw_clustermap(adata: sc.AnnData,
                    dataset: str,
                    cluster_name: str,
                    gene_list: list,
                    mean_expression_path: str,
                    clustering_dir: str):
    import seaborn as sns
    
    # Load cluster average non-zero expression
    df_mean = pd.read_csv(mean_expression_path, sep='\t', index_col=0, header=0)
    # Initiate figure and plot
    if len(df_mean.columns) > 1 and len(gene_list) > 1:
        
        fig_sizes = {'neurotransmitter': (len(df_mean.columns) * 0.5, len(gene_list) * 0.7),
                     'DEG': (len(df_mean.columns) * 0.5, len(gene_list) / 8)
                     }
        if dataset.find('neurotransmitter') != -1:
            fig_size = fig_sizes['neurotransmitter']
        elif dataset.find('DEG') != -1:
            fig_size = fig_sizes['DEG']
        sns.set_style("whitegrid", {'axes.grid': False})
        print("Draw clustermap...")
        clustermap = sns.clustermap(pd.DataFrame(np.expm1(df_mean), index=df_mean.index, columns=df_mean.columns),
                                    row_cluster=True,
                                    col_linkage=adata.uns[f'dendrogram_{cluster_name}']['linkage'],
                                    col_cluster=True,
                                    standard_scale=0,
                                    #z_score=0,
                                    method='ward',
                                    figsize=fig_size,
                                    #cmap='Blues',
                                    cmap="viridis",
                                    #cmap='turbo',
                                    yticklabels=True,
                                    annot_kws={"size": 80},      # size of axis labels
                                    cbar_pos=(1.05, .2, .03, .4),
                                    cbar_kws={'label': 'standard scale'}
                                    )
        clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xmajorticklabels(), fontsize=30)
        cax = clustermap.figure.axes[-1]
        cax.tick_params(labelsize=50)
        clustermap.fig.suptitle(f"Clustermap {dataset}_{cluster_name}", fontsize=30)
        
        dest = f"{clustering_dir}/Clustermap {dataset}_{cluster_name}.pdf"
        clustermap.savefig(dest)
        print(dest)
        plt.close()
        # Get gene dendrogram order
        labels = [lab.get_text() for lab in clustermap.ax_heatmap.get_yticklabels(which='both')]
        
        return labels

    
def sparsity(x: Union[np.ndarray, sparray]) -> float:
    """
    Parameters
    -----
    x: `array type`
        An array or a sparse matrix.
    
    Returns
    --------
    `float`
        Fraction of zero-valued entries
    """
    if issparse(x):
        nz = x.count_nonzero()
        return 1.0 - float(nz) / float(np.prod(x.shape))
    else:
        nz = np.count_nonzero(x)
        return 1.0 - float(nz) / float(x.size)


def get_significant_genes(adata: sc.AnnData,                 # AnnData
                          key: str,                          # Key used to store the ranking results in adata.uns
                          n_genes: int = 5,                  # Number of genes to show
                          pval_cutoff: float = 0.05,         # pvals_adj cutoff to consider
                          expression_cutoff: float = 0.4,    # Expression cutoff to consider
                          #fc_cutoff: float = 1.0             # LogFC cutoff  (TOO RESTRICTIVE)
                          ) -> list:
    '''\
        Collect n_genes that have highest logFC and pvals_adj < 'pval_cutoff'
        and expression cutoff > 'expression_cutoff'
        key: `str`
            Key used to store the ranking results in adata.uns
        n_genes: `int` = 5
            Number of genes to show
        pval_cutoff: `float` = 0.05
            pvals_adj cutoff to consider
        expression_cutoff: `float` = 0.4
            Expression cutoff to consider (fraction of cell sexpressing the gene)
        fc_cutoff: `float` = 1.0
            Log FC cutoff to consider

        Returns
        -------
        list((genes, n genes added, bool for pval)). Tuple containing list of genes, the size and if they were selected from significant pval.
        Values will be sorted by highest logFC, then expression
    '''
    import pandas as pd
    # Get size col numbers
    cols = pd.DataFrame(adata.uns[key]['names']).columns.to_list()
    gene_list = []
    for i, col in enumerate(cols):
        df = pd.DataFrame()
        df['names'] = pd.DataFrame(adata.uns[key]['names']).iloc[:, i]
        df['scores'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, i]
        df['logfoldchanges'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, i]
        df['pvals_adj'] = pd.DataFrame(adata.uns[key]['pvals_adj']).iloc[:, i]
        df['pts'] = pd.DataFrame(adata.uns[key]['pts']).loc[df['names'], str(col)].values
        df = df[df['pvals_adj'] < pval_cutoff].copy()
        df = df[df['pts'] >= expression_cutoff].copy()
        df.sort_values(by=['logfoldchanges', 'pts'], ascending=False, ignore_index=True, inplace=True)
        #df = df[df['logfoldchanges'] >= fc_cutoff].copy()
        values = df['names'].iloc[:n_genes].values
        gene_list.append([values.tolist(), len(values), True])  # Add True for significant pval
        # When no genes are found use highest score from within pts cutoff
        if gene_list[-1][1] == 0:
            df = pd.DataFrame()
            df['names'] = pd.DataFrame(adata.uns[key]['names']).iloc[:, i]
            df['scores'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, i]
            df['logfoldchanges'] = pd.DataFrame(adata.uns[key]['logfoldchanges']).iloc[:, i]
            df['pvals_adj'] = pd.DataFrame(adata.uns[key]['pvals_adj']).iloc[:, i]
            df['pts'] = pd.DataFrame(adata.uns[key]['pts']).loc[df['names'], str(col)].values
            df = df[df['pts'] >= expression_cutoff].copy()
            df.sort_values(by=['logfoldchanges', 'pts'], ascending=False, ignore_index=True, inplace=True)
            values = df['names'].iloc[:n_genes].values
            gene_list[-1] = [values.tolist(), len(values), False]   # Add False for significant pval

    return gene_list
    

def get_deg_genes(deg_path: str,
                  n_top_genes: None = None,
                  p_val_cutoff: float = 0.05,
                  expression_cutoff: float = 0.3    # Expression cutoff to consider
                  ) -> list:
    # The column selected is always top.xxx since this represents the comparison of
    # top NES score cells vs low NES score cells
    # If no n_top_genes are provided, all genes are returned,
    # else the top highest LogFC genes are returned
    columns = ['top.names', 'top.scores', 'top.logfoldchanges', 'top.pvals', 'top.pvals_adj', 'top.pts']
    deg_df = pd.read_csv(deg_path, sep='\t', header=0, index_col=0)
    deg_df = deg_df.loc[:, columns]
    deg_df = deg_df[deg_df.loc[:, 'top.pvals_adj'] < p_val_cutoff]
    deg_df = deg_df[deg_df.loc[:, 'top.pts'] > expression_cutoff]
    
    gene_list = []
    if n_top_genes is None:
        gene_list = deg_df['top.names'].values.tolist()
        return gene_list
    else:
        deg_df.sort_values(by='top.logfoldchanges', ascending=False, inplace=True)
        gene_list = deg_df['top.names'].values[:n_top_genes].tolist()
        return gene_list
    

def random_colors(n: int, seed: Union[int, None] = 0) -> list:
    import numpy as np
    from scipy.interpolate import interp1d
    from matplotlib.colors import is_color_like
    
    # n must be smaller than 256 * 256 * 256
    if not n < 256 * 256 * 256:
        raise ValueError("n must be smaller than 256 * 256 * 256")
        
    colors = []
    if seed is None:
        randomizer = np.random.default_rng()
    else:
        randomizer = np.random.default_rng(seed=seed)
    scaler = interp1d([0, 255], [0.0, 1.0])
    i = 0
    while i < n:
        rgb = randomizer.choice(a=256, size=3)
        rgb = tuple(scaler(rgb).tolist())
        while (rgb in colors) or (not is_color_like(rgb)):
            rgb = randomizer.choice(a=256, size=3)
            rgb = tuple(scaler(rgb).tolist())
        colors.append(rgb)
        i += 1
    
    return colors


def color_shades(n: int, cmap_name: str) -> list:
    """
    # Get RGB values for colormap and convert the colormap in
    # CAM02-UCS colorspace.  lab[0, :, 0] is the lightness.
    rgb = mpl.colormaps[cmap](x)[np.newaxis, :, :3]
    lab = cspace_converter("sRGB1", "CAM02-UCS")(rgb)
    """
    
    # Get cmap
    # Get lightness 50 for selected cmap
    #
    

def import_senescence_gene_sets() -> dict:
    import os
    import pandas as pd
    import mygene
    from .globals import gene_sets_senescence, gene_sets_dir
        
    paths = pd.read_csv(f"{gene_sets_dir}/import_gene_sets.csv", sep='\t', header=0, index_col=0)
    
    # Try to load final gene sets (mouse converted)
    for gs in paths.index:
        dest = f"{gene_sets_dir}/GeneSet_{gs.strip()}_final.csv"
        if os.path.exists(dest):
            with open(dest, 'r') as fp:
                # gene_sets_senescence is a global var
                gene_sets_senescence[gs.strip()] = [s.strip() for s in fp.readlines()]
    
    if any([v is None for v in gene_sets_senescence.values()]):
        for gs in paths.index:
            with open(f"{gene_sets_dir}/{paths.loc[gs.strip(), 'Path']}") as fp:
                gene_sets_senescence[gs.strip()] = [s.strip() for s in fp.readlines()]
        
        # Convert Cellage genesets form human to mouse
        cellage = ["CELLAGE", "CELLAGE_UP", "CELLAGE_DN"]
        mg = mygene.MyGeneInfo()
        for c in cellage:
            converted = mg.querymany(gene_sets_senescence[c], scopes='symbol,alias', species='mouse', fields='ensembl.gene,symbol', as_dataframe=True)
            converted.dropna(axis=0, subset='symbol', inplace=True)
            converted = converted['symbol']
            gene_sets_senescence[c] = converted.values.tolist()
        
        # Make sure no None values are present
        filtered = {k: v for k, v in gene_sets_senescence.items() if v is not None}
        gene_sets_senescence.clear()
        gene_sets_senescence.update(filtered)
        
        for k, v in gene_sets_senescence.items():
            dest = f"{gene_sets_dir}/GeneSet_{k}_final.csv"
            with open(dest, mode='w') as fp:
                fp.write('\n'.join(v))
    
    return gene_sets_senescence


def gene_sets_to_gsva() -> pd.DataFrame:
    gene_sets = import_senescence_gene_sets()
    # Create a DataFrame compatible with the GSVA algorithm
    df = {'description': [],
          'member': [],
          'name': []}
    for geneset, gene_list in gene_sets.items():
        for gene in gene_list:
            df['description'].append(pd.NA)
            df['member'].append(gene)
            df['name'].append(geneset)
    df = pd.DataFrame(df)
    
    return df


# For manipulation of CellPhoneDB results
def get_cell_types(cci: str,
                   major: bool      # Only major cell types, no clusters
                   
                   ) -> tuple[str, ...]:
    cell_types = cci.split('|')
    if major:
        cell_types[0] = cell_types[0].split('.')[0]
        cell_types[1] = cell_types[1].split('.')[0]
    
    return tuple(cell_types)


def segmented_palette(colors: tuple[str] = ("#000000", "#FFFFFF"),  # Hex colors list
                      segments: tuple = (0, 1),  # Relative position of each color
                      n: int = 4096,             # Number of quantization levels
                      name: str = "my_colormap"
                      ) -> tuple[matplotlib.colors.LinearSegmentedColormap, str]:
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import to_rgba
    
    if len(colors) != len(segments):
        raise ValueError(f"The number of colors provided ({len(colors)}) must be equal to \
                         the number of segments in the colormap ({len(segments)})")
    color_list = [(segments[i], to_rgba(colors[i])) for i in range(len(colors))]
    cmap = LinearSegmentedColormap.from_list(name, color_list, n)

    return cmap, name
