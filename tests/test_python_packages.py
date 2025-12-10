"""
Tests to verify that Python packages listed in env_desc.py are installed.
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path to import biomni
sys.path.insert(0, str(Path(__file__).parent.parent))

from biomni.env_desc import library_content_dict


def get_python_packages():
    """Extract Python packages from library_content_dict."""
    python_packages = {}
    for name, description in library_content_dict.items():
        if "[Python]" in description:
            python_packages[name] = description
    return python_packages


class TestPythonPackages:
    """Test that all Python packages in env_desc are installed."""

    # def test_python_packages_installed(self):
    #     """Loop through Python packages and check if they're importable."""
    #     python_packages = get_python_packages()
        
    #     failed_packages = []
        
    #     for package_name, description in python_packages.items():
    #         # Map package names to import names (handle some common cases)
    #         import_name = package_name.replace("-", "_")
            
    #         # Special mappings for packages where import name differs from package name
    #         special_mappings = {
    #             "biopython": "Bio",
    #             "biom-format": "biom",
    #             "scikit-bio": "skbio",
    #             "scikit-learn": "sklearn",
    #             "umap-learn": "umap",
    #             "faiss-cpu": "faiss",
    #             "harmony-pytorch": "harmony",
    #             # "cellxgene-census": "cellxgene_census",
    #             # "gget": "gget",
    #             # "scrublet": "scrublet",
    #             # "scvelo": "scvelo",
    #             # "pysam": "pysam",
    #             # "pyfaidx": "pyfaidx",
    #             # "pyranges": "pyranges",
    #             # "pybedtools": "pybedtools",
    #             # "statsmodels": "statsmodels",
    #             # "h5py": "h5py",
    #             # "tiledb": "tiledb",
    #             # "tiledbsoma": "tiledbsoma",
    #             # "pyscenic": "pyscenic",
    #             # "arboreto": "arboreto",
    #         }
            
    #         import_name = special_mappings.get(package_name, import_name)
            
    #         try:
    #             __import__(import_name)
    #         except ImportError as e:
    #             failed_packages.append(
    #                 f"{package_name} (import: {import_name}) - {str(e)}"
    #             )
        
    #     if failed_packages:
    #         pytest.fail(
    #             f"The following Python packages are not installed:\n"
    #             + "\n".join(f"  - {pkg}" for pkg in failed_packages)
    #         )

    @pytest.mark.parametrize(
        "package_name",
        list(get_python_packages().keys()),
        ids=lambda x: f"python_{x}"
    )
    def test_individual_python_package(self, package_name):
        """Test each Python package individually for detailed reporting."""
        import_name = package_name.replace("-", "_")
        
        special_mappings = {
            "biopython": "Bio",
            "biom-format": "biom",
            "scikit-bio": "skbio",
            "scikit-learn": "sklearn",
            "umap-learn": "umap",
            "faiss-cpu": "faiss",
            "harmony-pytorch": "harmony",
            # "cellxgene-census": "cellxgene_census",
            # "gget": "gget",
            # "scrublet": "scrublet",
            # "scvelo": "scvelo",
            # "pysam": "pysam",
            # "pyfaidx": "pyfaidx",
            # "pyranges": "pyranges",
            # "pybedtools": "pybedtools",
            # "statsmodels": "statsmodels",
            # "h5py": "h5py",
            # "tiledb": "tiledb",
            # "tiledbsoma": "tiledbsoma",
            # "pyscenic": "pyscenic",
            # "arboreto": "arboreto",
        }
        
        import_name = special_mappings.get(package_name, import_name)
        
        try:
            __import__(import_name)
        except ImportError as e:
            pytest.fail(
                f"Python package '{package_name}' (import as '{import_name}') is not installed.\n"
                f"Error: {str(e)}"
            )
