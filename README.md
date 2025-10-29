# Protein-Antibody-Tools (PAB Tools)

Utilities for protein and antibody sequence analysis.

## Installation

### Option 1: Using Conda (Recommended)

1. Create and activate the conda environment:
    ```sh
    conda env create -f conda-env.yaml
    conda activate protein-tools
    ```

2. Install the package:
    ```sh
    pip install git+https://github.com/ainbsci/protein-ab-tools.git
    ```

### Option 2: Manual Installation

1. Install `anarci` via conda (required dependency):
    ```sh
    conda install -c bioconda anarci
    ```

2. Install the package:
    ```sh
    pip install git+https://github.com/ainbsci/protein-ab-tools.git
    ```

**Note:** `anarci` must be installed via conda from the bioconda channel due to its complex dependencies (HMMER, muscle). It cannot be installed via pip.

3. Import package

    ```python
    import protein_ab_tools as pat
    ```

## How to use

### Align two proteins

```python
seq1 = 'MALWMRLLPLLALLALWGPDPAAA'
seq2 = 'MALWMRLLPLLALSSALWGPDPAAA'

pat.calc_percent_similarity(seq1, seq2)
```

### Run Anarci numbering

```python
pat.get_numbered_seq("QVQLVESGGGVVQPGRSLRLDCKASGITFSNSGMHWVRQAPGKGLEWVAVIWYDGSKRYYADSVKGRFTISRNSKNTLFLQMNSLRAEDTAVYYCATNDDYWGQGTLVTTVSS")
```

### Extract regions

```python
pat.extract_regions("QVQLVESGGGVVQPGRSLRLDCKASGITFSNSGMHWVRQAPGKGLEWVAVIWYDGSKRYYADSVKGRFTISRNSKNTLFLQMNSLRAEDTAVYYCATNDDYWGQGTLVTTVSS")
```
