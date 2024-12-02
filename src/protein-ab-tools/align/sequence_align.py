from Bio.Align import PairwiseAligner


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
    aligner = PairwiseAligner(mode)
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
