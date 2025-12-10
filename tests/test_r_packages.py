"""
Tests to verify that R packages listed in env_desc.py are installed.
"""

import pytest
import subprocess
import sys
from pathlib import Path

# Add parent directory to path to import biomni
sys.path.insert(0, str(Path(__file__).parent.parent))

from biomni.env_desc import library_content_dict


def get_r_packages():
    """Extract R packages from library_content_dict."""
    r_packages = {}
    for name, description in library_content_dict.items():
        if "[R]" in description:
            r_packages[name] = description
    return r_packages


def check_r_package_installed(package_name: str) -> tuple[bool, str]:
    """
    Check if an R package is installed using Rscript.
    
    Returns:
        (is_installed, error_message)
    """
    # R code to check if package is installed
    r_code = f'if (!require("{package_name}", quietly=TRUE)) {{ quit(status=1) }}'
    
    try:
        result = subprocess.run(
            ["Rscript", "-e", r_code],
            capture_output=True,
            timeout=10,
            text=True
        )
        
        if result.returncode == 0:
            return True, ""
        else:
            error_msg = result.stderr if result.stderr else "Package load failed"
            return False, error_msg
            
    except FileNotFoundError:
        return False, "Rscript not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "Rscript check timed out"
    except Exception as e:
        return False, str(e)


class TestRPackages:
    """Test that all R packages in env_desc are installed."""

    # def test_r_packages_installed(self):
    #     """Loop through R packages and check if they're installed."""
    #     r_packages = get_r_packages()
        
    #     failed_packages = []
    #     skipped_reason = None
        
    #     for package_name, description in r_packages.items():
    #         is_installed, error_msg = check_r_package_installed(package_name)
            
    #         # If Rscript is not available, skip all R tests with a reason
    #         if "Rscript not found" in error_msg:
    #             skipped_reason = error_msg
    #             break
            
    #         if not is_installed:
    #             failed_packages.append(f"{package_name} - {error_msg}")
        
    #     if skipped_reason:
    #         pytest.skip(f"R package tests skipped: {skipped_reason}")
        
    #     if failed_packages:
    #         pytest.fail(
    #             f"The following R packages are not installed:\n"
    #             + "\n".join(f"  - {pkg}" for pkg in failed_packages)
    #         )

    @pytest.mark.parametrize(
        "package_name",
        list(get_r_packages().keys()),
        ids=lambda x: f"r_{x}"
    )
    def test_individual_r_package(self, package_name):
        """Test each R package individually for detailed reporting."""
        is_installed, error_msg = check_r_package_installed(package_name)
        
        if "Rscript not found" in error_msg:
            pytest.skip(f"R package tests skipped: {error_msg}")
        
        if not is_installed:
            pytest.fail(
                f"R package '{package_name}' is not installed.\n"
                f"Error: {error_msg}\n\n"
                f"To install, run in R:\n"
                f"  install.packages('{package_name}')\n"
                f"or for Bioconductor packages:\n"
                f"  BiocManager::install('{package_name}')"
            )
