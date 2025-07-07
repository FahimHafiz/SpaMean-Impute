import os
import time
import tracemalloc
import scanpy as sc
from scipy.sparse import issparse, csr_matrix
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import (
    normalized_mutual_info_score,
    adjusted_rand_score,
    adjusted_mutual_info_score,
    homogeneity_score
)
from scipy.sparse import csr_matrix
import scipy.sparse as sp
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import spearmanr, pearsonr
from tqdm import tqdm
import pandas as pd



"""
This code is designed to impute zero values in spatial transcriptomics data 
by determining k-nearest spots (kNS) using spatial data.
Then, detects zero values in the nearest spots and using a threshold value, 
it determines whether to impute the zero values or not for a particular gene-spot combination.


It performs the following steps:
1. **Zero Indices Calculation**: Identifies genes with zero expression in each spot and stores their indices.

2. **kNN Spots Calculation**: Computes k nearest neighbors for each spots based on spatial coordinates.

3. **Merging Neighbors and Similarity Dictionaries**: Merges the kNN neighbors with the similarity metrics to create a comprehensive neighbor list for each spot.

4. **Zero Imputation**: Imputes zero values in the expression matrix using the average expression of the same gene in neighboring spots.
    - **Version 3**: Similar to Version 2, but allows for a drop threshold to skip genes with a high percentage of zeros in the neighbors.
                    This also converts sparse matrices to dense for speed, which may not be memory efficient for very large datasets.

5. **Sparsity Measurement**: Calculates the sparsity of the expression matrix before and after imputation.

7. **Clustering and Evaluation**: Performs clustering on the original and imputed data, and evaluates the clustering quality 
                                    using Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI).

8. **Output**: Saves the imputed data in h5ad format and prints evaluation metrics.

9. **Imputation Variants**: Implements multiple variants of the imputation function for flexibility and performance tuning.

"""

### Functions for zero indices and kNN neighbors

def zero_indices(adata_temp):
    """........."""
    adata_temp = adata_temp.copy()
    
    # Remove columns (genes) with all zero values
    nonzero_gene_mask = (adata_temp.X != 0).sum(axis=0).A1 if issparse(adata_temp.X) else (adata_temp.X != 0).sum(axis=0)
    # print(f"Number of non-zero genes: {nonzero_gene_mask}")
    adata_temp = adata_temp[:, nonzero_gene_mask > 0]
    # print(f"Number of genes after removing all-zero genes: {adata_temp.shape[1]}")

    # zero_dict_gene = {}
    zero_dict_spot = {}
    X = adata_temp.X
    # for idx, gene in enumerate(adata_temp.var_names):
    #     if issparse(X):
    #         zero_idx = (~X[:, idx].toarray().ravel().astype(bool)).nonzero()[0]
    #     else:
    #         zero_idx = (X[:, idx] == 0).nonzero()[0]
    #     zero_dict_gene[gene] = zero_idx.tolist()
    
    for idx, spot in enumerate(adata_temp.obs_names):
        if issparse(X):
            zero_idx = (~X[idx, :].toarray().ravel().astype(bool)).nonzero()[0]
        else:
            zero_idx = (X[idx, :] == 0).nonzero()[0]
        zero_dict_spot[spot] = zero_idx.tolist()
    
    return zero_dict_spot, adata_temp

def knn_neighbors_dict(adata_temp, k=100):
    """
    Compute k nearest neighbors for each cell based on spatial coordinates.

    Parameters:
        adata_temp: AnnData object with 'spatial' in obsm
        k: number of neighbors

    Returns:
        neighbors_dict: dict mapping cell id to list of nearest neighbor cell ids
    """
    adata_temp = adata_temp.copy()
    if 'spatial' not in adata_temp.obsm:
        raise ValueError("The AnnData object must contain spatial coordinates in 'obsm['spatial']'.")
    coords = adata_temp.obsm['spatial']
    cell_ids = np.array(adata_temp.obs_names)
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto').fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    neighbors_dict = {}
    for i, cell_id in enumerate(cell_ids):
        # Exclude the cell itself (first neighbor is always itself)
        neighbor_indices = indices[i][1:]
        neighbors_dict[cell_id] = cell_ids[neighbor_indices].tolist()
    return neighbors_dict


def impute_zeros_with_neighbors(zero_dict_spot, neighbors_dict, adata_temp, drop_threshold=0.5):
    """
    Fast imputation for dense matrices using vectorized operations.
    Converts sparse to dense for speed (if memory allows).
    """
    adata_temp = adata_temp.copy()
    X = adata_temp.X.copy()
    cell_idx_map = {cell: idx for idx, cell in enumerate(adata_temp.obs_names)}
    is_sparse = issparse(X)
    if is_sparse:
        X = X.toarray()  # Convert to dense for speed

    n_cells = len(zero_dict_spot)
    for i, (cell, zero_gene_indices) in enumerate(zero_dict_spot.items()):
        if not zero_gene_indices:
            continue
        cell_idx = cell_idx_map[cell]
        neighbors = neighbors_dict[cell]
        neighbor_indices = [cell_idx_map[n] for n in neighbors if n in cell_idx_map]
        if not neighbor_indices:
            continue
        
        neighbor_matrix = X[np.ix_(neighbor_indices, zero_gene_indices)]
        nonzero_mask = neighbor_matrix != 0

        # Calculate the percentage of zeros in each column
        zero_counts = (neighbor_matrix == 0).sum(axis=0)
        zero_percentage = zero_counts / neighbor_matrix.shape[0]

        # Keep only columns (genes) where <= 50% are zeros
        valid_gene_mask = zero_percentage <= drop_threshold
        if not np.any(valid_gene_mask):
            continue  # Skip if no genes to impute

        filtered_gene_indices = np.array(zero_gene_indices)[valid_gene_mask]
        filtered_neighbor_matrix = neighbor_matrix[:, valid_gene_mask]
        filtered_nonzero_mask = filtered_neighbor_matrix != 0

        # Compute means as before, but only for filtered genes
        sums = np.where(filtered_nonzero_mask, filtered_neighbor_matrix, 0).sum(axis=0)
        counts = filtered_nonzero_mask.sum(axis=0)
        means = np.divide(sums, counts, out=np.zeros_like(sums, dtype=float), where=counts != 0)
        X[cell_idx, filtered_gene_indices] = means
        if (i+1) % 200 == 0 or i == 0 or (i+1) == n_cells:
            print(f"Processed {i+1}/{n_cells} cells.")

    adata_temp.X = X
    return adata_temp


def calculate_sparsity(X):
        """Calculate sparsity of a matrix."""
        if issparse(X):
            X = X.toarray()
        zero_elements = np.sum(X == 0)
        total_elements = X.size
        return 100.0 * float(zero_elements) / float(total_elements)


def evaluate_imputation_and_clustering(adata, k_val, threshold_val):
    # Step 1: Select top highly variable genes (optional but recommended)
    adata_temp = adata.copy()

    # Step 2: Prepare zero value dictionary and neighbor dictionary
    zero_dict_spot, adata_filtered = zero_indices(adata_temp)
    spot_neighbors_dict = knn_neighbors_dict(adata_filtered, k=k_val)

    # Step 3: Impute
    tracemalloc.start()
    start_time = time.time()

    adata_imputed = impute_zeros_with_neighbors(zero_dict_spot, spot_neighbors_dict,
                                                   adata_filtered, drop_threshold=threshold_val)

    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    runtime = end_time - start_time
    memory = peak / (1024 ** 2)  # MB
    tracemalloc.stop()
    
    # Step 4: Measure sparsity
    sparsity_before = calculate_sparsity(adata.X)
    sparsity_after = calculate_sparsity(adata_imputed.X)

    # Step 5: Convert to sparse
    adata_imputed.X = csr_matrix(adata_imputed.X)

    # Step 6: Clustering on original
    adata_original = adata.copy()
    sc.pp.pca(adata_original)
    sc.pp.neighbors(adata_original)
    sc.tl.umap(adata_original)
    sc.tl.leiden(adata_original, key_added="clusters", directed=False, n_iterations=2)

    # Step 7: Clustering on imputed
    sc.pp.pca(adata_imputed)
    sc.pp.neighbors(adata_imputed)
    sc.tl.umap(adata_imputed)
    sc.tl.leiden(adata_imputed, key_added="clusters_imputed", directed=False, n_iterations=2)

    # Step 8: Compute metrics
    ari_original = adjusted_rand_score(adata_original.obs['clusters'], adata_original.obs['annotation'])
    ari_imputed = adjusted_rand_score(adata_imputed.obs['clusters_imputed'], adata_original.obs['annotation'])

    nmi_original = normalized_mutual_info_score(adata_original.obs['clusters'], adata_original.obs['annotation'])
    nmi_imputed = normalized_mutual_info_score(adata_imputed.obs['clusters_imputed'], adata_original.obs['annotation'])

    ami_original = adjusted_mutual_info_score(adata_original.obs['clusters'], adata_original.obs['annotation'])
    ami_imputed = adjusted_mutual_info_score(adata_imputed.obs['clusters_imputed'], adata_original.obs['annotation'])

    homo_original = homogeneity_score(adata_original.obs['clusters'], adata_original.obs['annotation'])
    homo_imputed = homogeneity_score(adata_imputed.obs['clusters_imputed'], adata_original.obs['annotation'])

    return {
        "k": k_val,
        "threshold": threshold_val,
        "sparsity_before": sparsity_before,
        "sparsity_after": sparsity_after,
        "ari_before": ari_original,
        "ari_after": ari_imputed,
        "nmi_before": nmi_original,
        "nmi_after": nmi_imputed,
        "ami_before": ami_original,
        "ami_after": ami_imputed,
        "homo_before": homo_original,
        "homo_after": homo_imputed,
        "runtime (s)": runtime,
        "memory (MB)": memory
    }

if __name__ == '__main__':
    # Main execution block

    # Parameters
    filename = 'GSE166692_scispace_processed'
    input_file = f"{filename}.h5ad"
    adata = sc.read_h5ad(input_file)
    adata_hvg = adata.copy()

    n_top = 5000  # Number of top genes to consider
    sc.pp.highly_variable_genes(adata_hvg, n_top_genes=n_top, flavor='seurat_v3')
    adata_hvg = adata_hvg[:, adata_hvg.var['highly_variable']]

    # Ranges to iterate
    k = 9
    threshold = 0.1  # 0.05 to 0.35 inclusive
    results = []

    print(f"Evaluating for {filename} where k={k}, threshold={threshold:.2f}...")
    result = evaluate_imputation_and_clustering(adata_hvg.copy(), k, threshold)
    results.append(result)

    print("SpaMean-Impute Done!")
    print(f"The results are: {results}")
