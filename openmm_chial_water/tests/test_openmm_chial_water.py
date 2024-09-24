"""
Unit and regression test for the openmm_chial_water package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import openmm_chial_water


def test_openmm_chial_water_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "openmm_chial_water" in sys.modules
