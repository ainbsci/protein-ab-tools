# Protein-Antibody-Tools (PAB Tools)

Utilities for protein and antibody sequence analysis.

## Installation

1. Install `biopython` and `anarci`
2. Run below.

    ```sh
    pip install git+https://github.com/ainbsci/protein-ab-tools.git
    ```

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
