"""
Tests for sequence alignment functionality.
"""
import pytest
from protein_ab_tools.align.sequence_align import calc_percent_similarity


class TestCalcPercentSimilarity:
    """Test the calc_percent_similarity function."""

    def test_identical_sequences(self):
        """Test that identical sequences have 100% similarity."""
        seq1 = 'MALWMRLLPLLALLALWGPDPAAA'
        seq2 = 'MALWMRLLPLLALLALWGPDPAAA'
        similarity = calc_percent_similarity(seq1, seq2)
        assert similarity == 100.0

    def test_example_from_readme(self):
        """Test the example from the README."""
        seq1 = 'MALWMRLLPLLALLALWGPDPAAA'
        seq2 = 'MALWMRLLPLLALSSALWGPDPAAA'
        similarity = calc_percent_similarity(seq1, seq2)
        # Should be high similarity but not 100%
        assert 80 < similarity < 100
        assert isinstance(similarity, float)

    def test_completely_different_sequences(self):
        """Test sequences with low similarity."""
        seq1 = 'AAAAAAAAAA'
        seq2 = 'GGGGGGGGGG'
        similarity = calc_percent_similarity(seq1, seq2)
        # Should have low similarity
        assert 0 <= similarity < 50
        assert isinstance(similarity, float)

    def test_different_lengths(self):
        """Test sequences of different lengths."""
        seq1 = 'MALWMRLLPLL'
        seq2 = 'MALWMRLLPLLALLALWGPDPAAA'
        similarity = calc_percent_similarity(seq1, seq2)
        # Should handle different lengths
        assert 0 <= similarity <= 100
        assert isinstance(similarity, float)

    def test_single_amino_acid(self):
        """Test with single amino acid sequences."""
        seq1 = 'M'
        seq2 = 'M'
        similarity = calc_percent_similarity(seq1, seq2)
        assert similarity == 100.0

    def test_single_vs_multiple(self):
        """Test single amino acid vs multiple."""
        seq1 = 'M'
        seq2 = 'MALW'
        similarity = calc_percent_similarity(seq1, seq2)
        assert 0 <= similarity <= 100
        assert isinstance(similarity, float)

    def test_with_gaps_removed(self):
        """Test that gaps should be removed before alignment."""
        # Note: The aligner expects ungapped sequences as input
        # Users should remove gaps before calling calc_percent_similarity
        seq1 = 'MAL-WMR'.replace('-', '')  # Remove gaps
        seq2 = 'MALWWMR'
        similarity = calc_percent_similarity(seq1, seq2)
        # Should handle sequences after gap removal
        assert 0 <= similarity <= 100
        assert isinstance(similarity, float)

    def test_empty_sequences(self):
        """Test behavior with empty sequences."""
        # This might raise an error or return 0, depending on biopython behavior
        seq1 = ''
        seq2 = ''
        try:
            similarity = calc_percent_similarity(seq1, seq2)
            # If it doesn't raise, check it returns a valid number
            assert isinstance(similarity, (int, float))
        except (ValueError, ZeroDivisionError):
            # Empty sequences might cause an error, which is acceptable
            pass

    def test_return_type(self):
        """Test that return type is always a float."""
        seq1 = 'MALWMRLLPLL'
        seq2 = 'MALWMRLLPLL'
        similarity = calc_percent_similarity(seq1, seq2)
        assert isinstance(similarity, float)

    def test_range_validation(self):
        """Test that similarity is always between 0 and 100."""
        test_cases = [
            ('MALWMRLLPLL', 'MALWMRLLPLL'),
            ('AAAA', 'GGGG'),
            ('ABC', 'ABCDEF'),
            ('QVQLVESGGG', 'QVQLVESGGG'),
        ]
        for seq1, seq2 in test_cases:
            similarity = calc_percent_similarity(seq1, seq2)
            assert 0 <= similarity <= 100, \
                f"Similarity should be 0-100, got {similarity} for {seq1} vs {seq2}"

    def test_symmetry(self):
        """Test that calc_percent_similarity(a, b) == calc_percent_similarity(b, a)."""
        seq1 = 'MALWMRLLPLLALLALWGPDPAAA'
        seq2 = 'MALWMRLLPLLALSSALWGPDPAAA'
        similarity1 = calc_percent_similarity(seq1, seq2)
        similarity2 = calc_percent_similarity(seq2, seq1)
        # Allow small floating point differences
        assert abs(similarity1 - similarity2) < 0.001

    def test_mode_parameter(self):
        """Test that mode parameter works."""
        seq1 = 'MALWMRLLPLL'
        seq2 = 'MALWMRLLPLL'
        # Default mode is 'blastp'
        similarity = calc_percent_similarity(seq1, seq2, mode='blastp')
        assert similarity == 100.0

    def test_common_antibody_sequences(self):
        """Test with real antibody sequences."""
        # Using the sequences from the numbering tests
        heavy = 'QVQLVESGGGVVQPGRSLRLDCKASGITFSNSGMHWVRQAPGKGLEWVAVIWYDGSKRYYADSVKGRFTISRNSKNTLFLQMNSLRAEDTAVYYCATNDDYWGQGTLVTTVSS'
        light = 'EIVLTQSPATLSLSPGERATLSCRASQSVSGYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQSSNWPRTFGQGTKVEIK'

        # Different chains should have low similarity
        similarity = calc_percent_similarity(heavy, light)
        assert 0 <= similarity <= 100
        assert isinstance(similarity, float)

        # Same sequence should have 100% similarity
        similarity_heavy = calc_percent_similarity(heavy, heavy)
        assert similarity_heavy == 100.0
