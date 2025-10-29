"""
Protein and Antibody Tools (PAT) - Utilities for protein and antibody sequence analysis.
"""
from .ab_analysis import (
    run_numbering,
    get_numbered_seq,
    extract_regions,
    extract_species,
    ANARCI_AVAILABLE
)
from .align import calc_percent_similarity

__all__ = [
    # Antibody numbering functions
    'run_numbering',
    'get_numbered_seq',
    'extract_regions',
    'extract_species',
    'ANARCI_AVAILABLE',
    # Sequence alignment functions
    'calc_percent_similarity',
]

__version__ = '0.0.1'
