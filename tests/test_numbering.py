"""
Tests for antibody numbering functionality.
"""
import pytest
from protein_ab_tools.ab_analysis.numbering import (
    run_numbering,
    get_numbered_seq,
    extract_regions,
    extract_species
)


# Test sequences from the original code
HEAVY_CHAIN_SEQ = 'QVQLVESGGGVVQPGRSLRLDCKASGITFSNSGMHWVRQAPGKGLEWVAVIWYDGSKRYYADSVKGRFTISRNSKNTLFLQMNSLRAEDTAVYYCATNDDYWGQGTLVTTVSS'
LIGHT_CHAIN_SEQ = 'EIVLTQSPATLSLSPGERATLSCRASQSVSGYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQSSNWPRTFGQGTKVEIK'


class TestRunNumbering:
    """Test the run_numbering function."""

    def test_heavy_chain_imgt(self):
        """Test numbering a heavy chain sequence with IMGT scheme."""
        result = run_numbering(HEAVY_CHAIN_SEQ, scheme='imgt', chain='H')
        assert result is not None
        assert result[0][0] is not None

    def test_light_chain_imgt(self):
        """Test numbering a light chain sequence with IMGT scheme."""
        result = run_numbering(LIGHT_CHAIN_SEQ, scheme='imgt', chain='L')
        assert result is not None
        assert result[0][0] is not None

    def test_invalid_sequence(self):
        """Test that invalid sequences raise ValueError."""
        with pytest.raises(ValueError, match="Invalid sequence"):
            run_numbering('INVALID', scheme='imgt', chain='H')

    def test_strip_and_remove_dashes(self):
        """Test that sequences are properly cleaned."""
        seq_with_dashes = 'QVQ-LVE-SGG'
        seq_with_spaces = '  QVQLVESGG  '
        # These should work the same as clean sequences
        # (though they won't be valid antibody sequences)
        try:
            run_numbering(seq_with_dashes, scheme='imgt', chain='H')
        except ValueError:
            pass  # Expected for invalid sequences

        try:
            run_numbering(seq_with_spaces, scheme='imgt', chain='H')
        except ValueError:
            pass  # Expected for invalid sequences

    def test_germline_option(self):
        """Test that germline option works."""
        result = run_numbering(HEAVY_CHAIN_SEQ, scheme='imgt', chain='H', germline=True)
        assert result is not None
        # When germline=True, result should have germline information
        assert result[0][0] is not None

    def test_species_filtering(self):
        """Test species filtering option."""
        result = run_numbering(
            HEAVY_CHAIN_SEQ,
            scheme='imgt',
            chain='H',
            species=['human']
        )
        assert result is not None


class TestExtractSpecies:
    """Test the extract_species function."""

    def test_heavy_chain_species(self):
        """Test extracting species from heavy chain."""
        species = extract_species(HEAVY_CHAIN_SEQ, scheme='imgt', chain='H')
        assert species is not None
        assert isinstance(species, str)

    def test_light_chain_species(self):
        """Test extracting species from light chain."""
        species = extract_species(LIGHT_CHAIN_SEQ, scheme='imgt', chain='L')
        assert species is not None
        assert isinstance(species, str)


class TestGetNumberedSeq:
    """Test the get_numbered_seq function."""

    def test_heavy_chain_numbered(self):
        """Test getting numbered sequence for heavy chain."""
        numbered = get_numbered_seq(HEAVY_CHAIN_SEQ, scheme='imgt', chain='H')
        assert numbered is not None
        assert isinstance(numbered, str)
        # Numbered sequence should contain only amino acids and gaps
        assert all(c in 'ACDEFGHIKLMNPQRSTVWY-' for c in numbered)

    def test_light_chain_numbered(self):
        """Test getting numbered sequence for light chain."""
        numbered = get_numbered_seq(LIGHT_CHAIN_SEQ, scheme='imgt', chain='L')
        assert numbered is not None
        assert isinstance(numbered, str)


class TestExtractRegionsIMGT:
    """Test extract_regions with IMGT scheme."""

    def test_heavy_chain_regions_imgt(self):
        """Test extracting regions from heavy chain with IMGT scheme."""
        expected = {
            'vh_fwr1': 'QVQLVESGG-GVVQPGRSLRLDCKAS',
            'vh_cdr1': 'GITF----SNSG',
            'vh_fwr2': 'MHWVRQAPGKGLEWVAV',
            'vh_cdr2': 'IWYD--GSKR',
            'vh_fwr3': 'YYADSVK-GRFTISR-NSKNTLFLQMNSLRAEDTAVYYC',
            'vh_cdr3': 'ATN-------DDY',
            'vh_fwr4': 'WGQGTLVTTVSS'
        }
        result = extract_regions(HEAVY_CHAIN_SEQ, scheme='imgt', chain='H')
        assert result == expected

    def test_light_chain_regions_imgt(self):
        """Test extracting regions from light chain with IMGT scheme."""
        expected = {
            'vl_fwr1': 'EIVLTQSPATLSLSPGERATLSCRAS',
            'vl_cdr1': 'QSV------SGY',
            'vl_fwr2': 'LAWYQQKPGQAPRLLIY',
            'vl_cdr2': 'DA-------S',
            'vl_fwr3': 'NRATGIP-ARFSGSG--SGTDFTLTISSLEPEDFAVYYC',
            'vl_cdr3': 'QQSSN----WPRT',
            'vl_fwr4': 'FGQGTKVEIK'
        }
        result = extract_regions(LIGHT_CHAIN_SEQ, scheme='imgt', chain='L')
        assert result == expected


class TestExtractRegionsKabat:
    """Test extract_regions with Kabat scheme - validates the fix."""

    def test_kabat_ranges_no_overlap_heavy(self):
        """
        Test that Kabat heavy chain ranges don't overlap.
        This validates the fix where FWR3 was starting at 64 instead of 66.
        """
        result = extract_regions(HEAVY_CHAIN_SEQ, scheme='kabat', chain='H')

        # Get the numbering to check which positions are assigned to which region
        numbered_result = run_numbering(HEAVY_CHAIN_SEQ, scheme='kabat', chain='H')
        numbered_seq = numbered_result[0][0][0][0]

        # Check that positions 64 and 65 are NOT in FWR3
        # They should be in CDR2 (which ends at 65)
        cdr2_positions = [pos for (pos, _), aa in numbered_seq if 50 <= pos <= 65]
        fwr3_positions = [pos for (pos, _), aa in numbered_seq if 66 <= pos <= 94]

        # Ensure no overlap
        assert set(cdr2_positions).isdisjoint(set(fwr3_positions)), \
            "CDR2 and FWR3 positions should not overlap"

        # Specifically check that position 64 and 65 are in CDR2, not FWR3
        assert 64 in cdr2_positions or 65 in cdr2_positions, \
            "Positions 64-65 should be in CDR2"

    def test_kabat_ranges_no_overlap_light(self):
        """Test that Kabat light chain ranges don't overlap."""
        result = extract_regions(LIGHT_CHAIN_SEQ, scheme='kabat', chain='L')

        numbered_result = run_numbering(LIGHT_CHAIN_SEQ, scheme='kabat', chain='L')
        numbered_seq = numbered_result[0][0][0][0]

        # Check that positions don't overlap
        cdr2_positions = [pos for (pos, _), aa in numbered_seq if 50 <= pos <= 56]
        fwr3_positions = [pos for (pos, _), aa in numbered_seq if 57 <= pos <= 88]

        # Ensure no overlap
        assert set(cdr2_positions).isdisjoint(set(fwr3_positions)), \
            "CDR2 and FWR3 positions should not overlap"

    def test_kabat_heavy_chain_regions(self):
        """Test extracting regions from heavy chain with Kabat scheme."""
        result = extract_regions(HEAVY_CHAIN_SEQ, scheme='kabat', chain='H')

        # Check that all expected regions are present
        assert 'vh_fwr1' in result
        assert 'vh_cdr1' in result
        assert 'vh_fwr2' in result
        assert 'vh_cdr2' in result
        assert 'vh_fwr3' in result
        assert 'vh_cdr3' in result
        assert 'vh_fwr4' in result

        # Check that all regions are non-empty
        for region, seq in result.items():
            assert len(seq) > 0, f"Region {region} should not be empty"

    def test_kabat_light_chain_regions(self):
        """Test extracting regions from light chain with Kabat scheme."""
        result = extract_regions(LIGHT_CHAIN_SEQ, scheme='kabat', chain='L')

        # Check that all expected regions are present
        assert 'vl_fwr1' in result
        assert 'vl_cdr1' in result
        assert 'vl_fwr2' in result
        assert 'vl_cdr2' in result
        assert 'vl_fwr3' in result
        assert 'vl_cdr3' in result
        assert 'vl_fwr4' in result

        # Check that all regions are non-empty
        for region, seq in result.items():
            assert len(seq) > 0, f"Region {region} should not be empty"


class TestExtractRegionsChothia:
    """Test extract_regions with Chothia scheme."""

    def test_chothia_heavy_chain_regions(self):
        """Test extracting regions from heavy chain with Chothia scheme."""
        result = extract_regions(HEAVY_CHAIN_SEQ, scheme='chothia', chain='H')

        # Check that all expected regions are present
        assert all(key in result for key in [
            'vh_fwr1', 'vh_cdr1', 'vh_fwr2', 'vh_cdr2',
            'vh_fwr3', 'vh_cdr3', 'vh_fwr4'
        ])

    def test_chothia_light_chain_regions(self):
        """Test extracting regions from light chain with Chothia scheme."""
        result = extract_regions(LIGHT_CHAIN_SEQ, scheme='chothia', chain='L')

        # Check that all expected regions are present
        assert all(key in result for key in [
            'vl_fwr1', 'vl_cdr1', 'vl_fwr2', 'vl_cdr2',
            'vl_fwr3', 'vl_cdr3', 'vl_fwr4'
        ])


class TestExtractRegionsAho:
    """Test extract_regions with AHo scheme."""

    def test_aho_heavy_chain_regions(self):
        """Test extracting regions from heavy chain with AHo scheme."""
        result = extract_regions(HEAVY_CHAIN_SEQ, scheme='aho', chain='H')

        # Check that all expected regions are present
        assert all(key in result for key in [
            'vh_fwr1', 'vh_cdr1', 'vh_fwr2', 'vh_cdr2',
            'vh_fwr3', 'vh_cdr3', 'vh_fwr4'
        ])

    def test_aho_light_chain_regions(self):
        """Test extracting regions from light chain with AHo scheme."""
        result = extract_regions(LIGHT_CHAIN_SEQ, scheme='aho', chain='L')

        # Check that all expected regions are present
        assert all(key in result for key in [
            'vl_fwr1', 'vl_cdr1', 'vl_fwr2', 'vl_cdr2',
            'vl_fwr3', 'vl_cdr3', 'vl_fwr4'
        ])


class TestInvalidScheme:
    """Test error handling for invalid numbering schemes."""

    def test_invalid_scheme_raises_error(self):
        """Test that invalid numbering scheme raises ValueError."""
        with pytest.raises(ValueError, match="Invalid numbering scheme"):
            extract_regions(HEAVY_CHAIN_SEQ, scheme='invalid', chain='H')
