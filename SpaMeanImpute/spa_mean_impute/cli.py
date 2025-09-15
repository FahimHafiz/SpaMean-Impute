import argparse
from .imputer import run_spamean_impute

def main():
    parser = argparse.ArgumentParser(description="SpaMean-Impute: Spatial Transcriptomics Imputer")
    parser.add_argument("input_file", help="Input .h5ad file")
    parser.add_argument("-k", type=int, default=9, help="Number of spatial neighbors")
    parser.add_argument("-t", "--threshold", type=float, default=0.1, help="Drop threshold for zero imputation")
    parser.add_argument("-n", "--n_top", type=int, default=5000, help="Number of highly variable genes")
    parser.add_argument("-o", "--output_file", default=None, help="Output .h5ad file")
    args = parser.parse_args()
    result = run_spamean_impute(
        input_file=args.input_file,
        k=args.k,
        threshold=args.threshold,
        n_top=args.n_top,
        output_file=args.output_file
    )
    print("SpaMean-Imputation results:", result)