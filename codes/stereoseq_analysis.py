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
dataset = 'DT2_D0_stereo-seq'
output_h5ad = dataset + '_processed'
dataset_path = 'D:/Research/spatial_transcriptomics/Data/Final Data with ground truth/stereo-seq/' + str(dataset) + '.h5ad' #please replace 'file_fold' with the download path
output_path = 'D:/Research/spatial_transcriptomics/Data/Final Data with ground truth/stereo-seq/'+ str(output_h5ad) + '.h5ad'
# Load the Stereo-seq data
adata = sc.read_h5ad(dataset_path)
adata.var_names_make_unique()

# Exploring the data in details!

print(f"The spatial data in AnnData format:\n {adata}\n\n")
print(f"Spot's row and column index:\n {adata.obs.head()}\n\n")  # Print first 5 rows
print(f"Tissue gene information:\n {adata.var}\n\n")  
print(f"Gene Names:\n {adata.var.index}\n\n")  # Names of genes

# View Expression Matrix
print(f"Spot-wise Gene Expression Matrix:\n {adata.X}\n\n")

df_sge_mat = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
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
sc.pp.filter_cells(adata, min_counts=50)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=5)


# Normalization and Feature Selection
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# Save AnnData object
adata.write(output_path)