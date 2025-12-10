"""
Tests to verify that CLI tools listed in env_desc.py are installed.
"""

import pytest
import subprocess
import sys
from pathlib import Path

# Add parent directory to path to import biomni
sys.path.insert(0, str(Path(__file__).parent.parent))

from biomni.env_desc import library_content_dict


def get_cli_tools():
    """Extract CLI tools from library_content_dict."""
    cli_tools = {}
    for name, description in library_content_dict.items():
        if "[CLI]" in description:
            cli_tools[name] = description
    return cli_tools


def check_cli_tool_installed(tool_name: str) -> tuple[bool, str]:
    """
    Check if a CLI tool is installed using 'which' or '--version'.
    
    Returns:
        (is_installed, error_message)
    """
    # Try 'which' command first
    try:
        result = subprocess.run(
            ["which", tool_name],
            capture_output=True,
            timeout=5,
            text=True
        )
        
        if result.returncode == 0 and result.stdout.strip():
            return True, result.stdout.strip()
    except Exception:
        pass
    
    # Try running the tool with --version or -version
    for version_flag in ["--version", "-version", "-v"]:
        try:
            result = subprocess.run(
                [tool_name, version_flag],
                capture_output=True,
                timeout=5,
                text=True
            )
            
            if result.returncode == 0:
                return True, f"Found via {version_flag}"
        except FileNotFoundError:
            return False, f"{tool_name} not found in PATH"
        except subprocess.TimeoutExpired:
            # If it times out, tool exists but maybe doesn't respond to version flag
            return True, "Found (timed out on version check)"
        except Exception:
            continue
    
    return False, f"{tool_name} not found or not responding"


class TestCLITools:
    """Test that all CLI tools in env_desc are installed."""

    def test_cli_tools_installed(self):
        """Loop through CLI tools and check if they're available."""
        cli_tools = get_cli_tools()
        
        if not cli_tools:
            pytest.skip("No CLI tools found in library_content_dict")
        
        failed_tools = []
        
        for tool_name, description in cli_tools.items():
            is_installed, result = check_cli_tool_installed(tool_name)
            
            if not is_installed:
                failed_tools.append(f"{tool_name} - {result}")
        
        if failed_tools:
            pytest.fail(
                f"The following CLI tools are not installed:\n"
                + "\n".join(f"  - {tool}" for tool in failed_tools)
            )

    @pytest.mark.parametrize(
        "tool_name",
        list(get_cli_tools().keys()),
        ids=lambda x: f"cli_{x}"
    )
    def test_individual_cli_tool(self, tool_name):
        """Test each CLI tool individually for detailed reporting."""
        is_installed, result = check_cli_tool_installed(tool_name)
        
        if not is_installed:
            pytest.fail(
                f"CLI tool '{tool_name}' is not installed.\n"
                f"Error: {result}\n\n"
                f"The tool should be available in your PATH.\n"
                f"Check installation instructions for '{tool_name}'."
            )
