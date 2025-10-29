# Tests for Protein-Antibody-Tools

This directory contains the test suite for the protein-ab-tools package.

## Test Organization

- `test_sequence_align.py` - Tests for sequence alignment functionality
- `test_numbering.py` - Tests for antibody numbering functionality (requires anarci)
- `conftest.py` - Pytest configuration and fixtures

## Running Tests

### Without anarci (limited testing)

If you don't have anarci installed, you can still run the sequence alignment tests:

```bash
pytest tests/
```

This will:
- Run all sequence alignment tests (13 tests)
- Skip all numbering tests (21 tests) with a message explaining that anarci is required

### With anarci (full testing)

To run all tests including the antibody numbering tests, first install anarci via conda:

```bash
# Create and activate conda environment
conda env create -f conda-env.yaml
conda activate protein-tools

# Install the package in development mode
pip install -e ".[dev]"

# Run all tests
pytest tests/
```

This will run all 34 tests.

## Test Categories

### Sequence Alignment Tests (test_sequence_align.py)
- Tests for `calc_percent_similarity()` function
- Validates alignment calculations
- Tests edge cases (empty sequences, different lengths, etc.)
- **Does not require anarci**

### Antibody Numbering Tests (test_numbering.py)
- Tests for all numbering schemes (IMGT, Kabat, Chothia, AHo)
- **Validates the Kabat FWR3 fix** (positions 66-94 for heavy chain)
- Tests region extraction and CDR/framework boundaries
- Tests species detection and germline assignment
- **Requires anarci to be installed**

## Key Tests for Bug Fixes

The following tests specifically validate the bug fixes made during code review:

1. **`test_kabat_ranges_no_overlap_heavy`** - Validates that the Kabat heavy chain FWR3 starts at position 66 (not 64), eliminating overlap with CDR2
2. **`test_kabat_ranges_no_overlap_light`** - Validates that light chain regions don't overlap

## Running Specific Tests

```bash
# Run only sequence alignment tests
pytest tests/test_sequence_align.py

# Run only numbering tests (requires anarci)
pytest tests/test_numbering.py

# Run tests with verbose output
pytest tests/ -v

# Run a specific test class
pytest tests/test_numbering.py::TestExtractRegionsKabat -v

# Run a specific test
pytest tests/test_numbering.py::TestExtractRegionsKabat::test_kabat_ranges_no_overlap_heavy -v
```

## Test Coverage

The test suite provides comprehensive coverage of:
- All public API functions
- All numbering schemes (IMGT, Kabat, Chothia, AHo)
- Both heavy and light chain sequences
- Error handling and edge cases
- Validation of region boundaries and overlaps

## CI/CD Note

When running tests in CI/CD without conda, the numbering tests will be automatically skipped. This is expected behavior. For full test coverage in CI/CD, ensure anarci is installed via conda before running pytest.
