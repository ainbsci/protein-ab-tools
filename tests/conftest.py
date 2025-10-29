"""
Pytest configuration and fixtures.
"""
import pytest


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers",
        "requires_anarci: marks tests that require anarci to be installed"
    )


def pytest_collection_modifyitems(config, items):
    """Skip tests that require anarci if it's not installed."""
    try:
        import anarci
        anarci_available = True
    except ImportError:
        anarci_available = False

    if not anarci_available:
        skip_anarci = pytest.mark.skip(reason="anarci not installed (install via conda)")
        for item in items:
            # Skip all tests in test_numbering.py as they require anarci
            if "test_numbering" in item.nodeid:
                item.add_marker(skip_anarci)
