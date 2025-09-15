# SpaMeanImpute

**Spatial Transcriptomics Imputation Tool (SpaMean-Impute)**

SpaMeanImpute is a Python package for imputing missing gene expression values in spatial transcriptomics datasets.  
It leverages spatial neighborhood information to improve downstream clustering and analysis.

---

## ✨ Features

- **Spatially-aware imputation** using k-nearest neighbors in spatial coordinates.
- **Highly variable gene selection** with `scanpy`.
- **Evaluation pipeline** for imputation and clustering performance.
- **Command-line interface (CLI)** for easy integration into workflows.
- **Compatible with `.h5ad` AnnData files**.

---

## 📦 Installation

You can install SpaMeanImpute from PyPI:

```bash
pip install spameanimpute
```
Or install directly from the GitHub repository:
```bash
git clone https://github.com/FahimHafiz/SpaMean-Impute.git
cd SpaMean-Impute
pip install -e .
```

## 🚀 Usage

After installation, you can run:
```bash
spamean-impute input_file.h5ad \
    -k 9 \
    --threshold 0.1 \
    --n_top 5000 \
    --output_file imputed_output.h5ad
```

Arguments:

input_file — Path to input .h5ad file.

-k — Number of spatial neighbors (default: 9).

--threshold — Drop threshold for zero imputation (default: 0.1).

--n_top — Number of highly variable genes (default: 5000, use 'all' to skip filtering).

--output_file — Path to save the imputed .h5ad file.


## 🐍 Python API

```python
from spa_mean_impute.imputer import SpaMeanImpute

imputer = SpaMeanImpute(k=9, threshold=0.1, n_top=5000)
results = imputer.run("input_file.h5ad", output_file="imputed_output.h5ad")

print("Results:", results)
```


## 📁 Project Structure

spa_mean_impute/
    __init__.py
    imputer.py
    cli.py
    utils.py
tests/
    test_imputer.py
setup.py
pyproject.toml
README.md
LICENSE.txt


## ✅ Running Tests

python -m unittest discover tests


## 📄 License

This project is licensed under the MIT License — see the LICENSE.txt file for details. For any problem, feel free to contact [fahimhafiz@cse.uiu.ac.bd](fahimhafiz@cse.uiu.ac.bd) or create an issue!