import unittest
import os
from spa_mean_impute.imputer import SpaMeanImpute
import scanpy as sc

class TestSpaMeanImpute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.real_data_path = r"D:\VM Data\thesis\Analysis\proposed imputation method\151507_visium_processed.h5ad"
        if not os.path.exists(cls.real_data_path):
            raise FileNotFoundError(f"Real data file not found: {cls.real_data_path}")

        adata = sc.read_h5ad(cls.real_data_path)
        if 'spatial' not in adata.obsm:
            raise ValueError("The provided AnnData file does not contain obsm['spatial'].")

    def tearDown(self):
        # Remove output file if it exists
        if os.path.exists("real_output.h5ad"):
            os.remove("real_output.h5ad")

    def test_run_returns_list(self):
        imputer = SpaMeanImpute(k=3, threshold=0.1, n_top='all')
        results = imputer.run(self.real_data_path, output_file="real_output.h5ad")
        self.assertIsInstance(results, list)
        self.assertGreater(len(results), 0)

    def test_output_file_created(self):
        imputer = SpaMeanImpute(k=3, threshold=0.1, n_top='all')
        imputer.run(self.real_data_path, output_file="real_output.h5ad")
        self.assertTrue(os.path.exists("real_output.h5ad"))

if __name__ == "__main__":
    unittest.main()
