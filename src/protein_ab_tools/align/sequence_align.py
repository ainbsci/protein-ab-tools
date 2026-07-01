from functools import lru_cache

from Bio.Align import PairwiseAligner


@lru_cache(maxsize=None)
def _get_aligner(mode):
    """Aligner setup (substitution matrix, gap penalties) is identical for a
    given mode on every call, so build it once per mode and reuse it — this
    function runs in tight loops scoring millions of sequence pairs."""
    return PairwiseAligner(mode)


def calc_percent_similarity(seq1, seq2, mode='blastp'):
    """
    Calculate the percent similarity between two sequences.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        mode (str): The alignment mode. Default is 'blastp' for protein sequences.

    Returns:
        float: The percent similarity between the two sequences.
    """
    aligner = _get_aligner(mode)
    alignments = aligner.align(seq1, seq2)
    counts = alignments[0].counts()
    identities = counts.identities
    gaps = counts.gaps
    mismatches = counts.mismatches
    return identities / (identities + gaps + mismatches) * 100


if __name__ == '__main__':
    seq1 = 'MALWMRLLPLLALLALWGPDPAAA'
    seq2 = 'MALWMRLLPLLALSSALWGPDPAAA'
    print(calc_percent_similarity(seq1, seq2))
