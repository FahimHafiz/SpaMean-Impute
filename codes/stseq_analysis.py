import os
import torch
import pandas as pd
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics

"""
Useful links:
1) https://scanpy-tutorials.readthedocs.io/en/multiomics/analysis-visualization-spatial.html
2) https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html


Stereo-seq Dataset Links with ground truth:
1) https://db.cngb.org/stomics/lista/download/


Directory:
1) D:/Research/spatial_transcriptomics/Data/Final Data with ground truth/xyz-seq/
"""

# scanpy settings
sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


## Read spatial transcriptomics data--> h5ad Data

# dataset_name = '151507'
# DX6_D2_stereo-seq.h5ad, DT2_D0_stereo-seq.h5ad, FB2_D1_stereo-seq.h5ad
dataset = 'GSM5009529_XYZeq_raw.h5ad'
file_fold = 'D:/Research/spatial_transcriptomics/Data/Final Data with ground truth/xyz-seq/' + str(dataset) #please replace 'file_fold' with the download path

# Load the Stereo-seq data
adata = sc.read_h5ad(file_fold)
adata.var_names_make_unique()

# Exploring the data in details!

print(f"The spatial data in AnnData format:\n {adata}\n\n")
print(f"Spot's row and column index:\n {adata.obs.head()}\n\n")  # Print first 5 rows
print(f"Tissue gene information:\n {adata.var}\n\n")  
print(f"Gene Names:\n {adata.var.index}\n\n")  # Names of genes

# View Expression Matrix
print(f"Spot-wise Gene Expression Matrix:\n {adata.X.toarray()}\n\n")

df_sge_mat = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
print(f"Spot-wise Gene Expression Matrix in dataframe format:\n {df_sge_mat.head()}")
# df_sge_mat.to_csv("V1_Human_Lymph_Node_sge.csv", index=False)

# Statistics of the data
print(f"Spot Statistics:\n {adata.obs.describe()}\n\n")

# End of Dataset Exploration

# Distribution of Gene Counts per Spot
# sns.histplot(adata.obs["n_genes_by_counts"], bins=50, kde=True)
# plt.xlabel("Number of Genes")
# plt.ylabel("Frequency")
# plt.title("Distribution of Gene Counts per Spot")
# plt.show()


# QC and preprocessing

"""
The following adata.var_names.str.startswith calculate the followings:
i) Identifies mitochondrial genes in the dataset
ii) str.startswith("MT-") checks if the gene name starts with "MT-", which marks mitochondrial genes.
iii) The result (True/False) is stored as a new column "mt" in adata.var
"""
adata.var["mt"] = adata.var_names.str.startswith("MT-")

"""
The following sc.pp.calculate_qc_metrices calculate the followings:
i) Total counts per cell (total RNA molecules detected per cell)
ii) Number of genes detected per cell
iii) Percentage of mitochondrial gene counts (pct_counts_mt)
This helps in filtering low-quality cells and 
The results are stored in adata.obs (observations metadata)
"""
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

"""
The following plot visualize the following distributions:
i) Total counts per cell (RNA molecules per cell)
ii) Total counts per cell (cells with fewer than 10,000 counts)
iii) Number of genes detected per cell
iv) Number of genes detected per cell (cells with fewer than 4,000 genes detected)

Helps identify:
i) Low-quality cells with low gene detection
ii) Cells with abnormally high counts (potential doublets)
iii) The second and fourth plots focus on lower count thresholds for better visibility.
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
"""


fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)

# Filtering Low-Quality Cells and Genes
"""
The below codes perform the followings:
i) Cells with less than 5000 total counts are removed.
ii) Cells with more than 35000 total counts are removed.
iii) Cells with >20% mitochondrial gene expression are removed.
iv) Only genes that appear in at least 10 cells are retained.
"""
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)


# Normalization and Feature Selection
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=200)

# Dimensionality Reduction and Clustering
"""
The below codes perform the followings:
i) Reduce data complexity using PCA.
ii) Construct a graph of nearest neighbors.
iii) Apply UMAP for visualization.
iv) Perform Leiden clustering to identify cell groups and stores them in "clusters".
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
"""
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", directed=False, n_iterations=2)

# ARI Calculation
ARI = metrics.adjusted_rand_score(adata.obs['clusters'], adata.obs['CellType'])
adata.uns['ARI'] = ARI
print('ARI:', ARI)

# UMAP Visualization of the spots considering "total_counts", "n_genes_by_counts" and "cluster number"

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.1)
























##################################### To Be Explore Later
# Visualization in spatial coordinates i.e. Visualizing Spatially Clustered Cells

plt.rcParams["figure.figsize"] = (8, 8)



"""
The below codes perform the followings:
i) The spatial function from Scanpy (sc.pl.spatial) overlays metadata onto a spatial tissue image.
ii) The argument img_key="hires" tells Scanpy to use a high-resolution image of the tissue.
iii) It colors the spatial spots based on:
"total_counts" → The total RNA molecules detected per spot.
"n_genes_by_counts" → The number of genes detected per spot.
iv) Overlays Leiden clustering results onto the spatial image.
"""
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

# Highlighting Specific Clusters in Spatial Coordinates
"""
The below codes perform the followings:
i) Highlights only clusters 5 and 9 in the spatial image.
ii) crop_coord=[7000, 10000, 0, 6000] → Defines a specific region of interest in spatial coordinates.
"""
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["1", "4"], # groups=["5", "9"],
    crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,
    size=1.3,
)

# Cluster marker genes i.e. Identifying Marker Genes for Each Cluster
"""
The below codes perform the followings:
i) Performs differential gene expression analysis across clusters.
ii) "clusters" → Uses cluster labels as grouping variables.
iii) method="t-test" → Uses the t-test to compare expression differences.
iv) Identifies marker genes that are uniquely expressed in each cluster and added to `.uns['rank_genes_groups']`
"""
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")

# Visualizing Marker Genes with a Heatmap
"""
The below codes perform the followings:
i) Plots a heatmap of the top 10 marker genes for cluster 9.
ii) groupby="clusters" → Organizes the heatmap by cluster labels.
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
"""
sc.pl.rank_genes_groups_heatmap(adata, groups="2", n_genes=10, groupby="clusters")

# Overlaying Specific Gene Expressions on the Spatial Image
"""
The below codes perform the followings:
i) Overlays both cluster assignments and the gene "CR2" onto the spatial tissue image.
ii) Allows to check whether "CR2" expression correlates with any cluster.
"""
sc.pl.spatial(adata, img_key="hires", color=["clusters", "COL1A2"])

# Visualizing the Expression of Specific Genes
"""
The below codes perform the followings:
i) Overlays the expression levels of genes "COL1A2" and "SYPL1".
ii) Helps in understanding whether these genes are spatially restricted to certain regions.
"""
sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)
