# ST_impute_benchmarking
 A benchmarking analysis of scRNA-seq and Spatial transcriptomics (ST) based imputations on different Spatial Transcriptomics Data

# Gene Selection, Imputation & Clustering Pipeline

## Overview
This repository provides a complete pipeline for gene selection, imputation, and clustering using multiple imputation methods and Leiden clustering. The workflow is designed for flexibility, allowing users to process datasets efficiently based on their computational power.

### Features
- **Gene Selection Criteria**: 
  - All genes
  - Top 2000 genes
  - Top 5000 genes
- **Imputation Methods**:
  - MAGIC
  - KNN
  - Soft Impute
  - Simple Impute
  - scVI
  - gimVI
  - Tangram
  - More imputation Methods (Easily extendable within the `Imputation Evaluator` class)
- **Clustering**: Leiden clustering is utilized for inferencing.
- **Batch Processing Modes**:
  - **Full Batch Mode**: Processes all `.h5ad` files in a directory, executing runs for all genes, top 2000 genes, and top 5000 genes before saving results in a `.csv` file.
  - **Semi Batch Mode**: Allows processing of a specified number of datasets at a time.
  - **Multi-batch Mode**: Enables parameterized batch processing, where:
    - First parameter determines the number of datasets processed at once.
    - Second parameter controls the number of gene selection criteria to apply.

The flexible batch modes ensure that different computational setups can process datasets optimally without excessive delays.

## Dataset Information
The following datasets are used:
- Visium
- Stereo-seq
- Slide-seq
- SciSpace
- XYZeq

Additionally, key evaluation metrics such as **ARI, NMI, HOMO, AMI, zero sparsity, runtime, and memory usage** are calculated for each imputation method.

### Dataset Links
You can access the datasets here:
- **Processed Datasets Link**: [Dataset Link](https://drive.google.com/drive/folders/1mNmJe9xVNpLtMlJGBOsdc9aJEleVhr36?usp=sharing)


## Usage Instructions
1. Clone the repository:
   ```sh
   git clone <repo-url>

## License
This project is licensed under the MIT License.

Enjoy using the pipeline, and investigate your imputation analysis with ease!