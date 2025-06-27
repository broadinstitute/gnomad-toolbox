"""Tests for ClinVar filtering functionality."""

import pytest
import hail as hl

from gnomad_toolbox.filtering.clinvar import filter_variants_by_clinvar


class TestClinVarFiltering:
    """Test class for ClinVar filtering functionality."""

    def test_filter_variants_by_clinvar_basic_functionality(self) -> None:
        """Test basic functionality of filter_variants_by_clinvar."""
        # Test with a small interval to avoid loading too much data.
        result = filter_variants_by_clinvar(
            # Part of BRCA1 gene interval.
            intervals="chr17:43044295-43044305",
            max_af_threshold=0.01,
            download_clinvar=False,
        )

        # Check that we get a Hail Table.
        assert isinstance(result, hl.Table)

        # Check that the result has ClinVar annotations.
        if result.count() > 0:
            assert "clinvar" in result.row, "Expected 'clinvar' annotation in result"

            # Check that clinvar annotation has expected fields.
            clinvar_fields = ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNDISDB"]
            for field in clinvar_fields:
                assert field in result.clinvar, f"Expected field {field} in clinvar annotation"

    def test_filter_variants_by_clinvar_with_genes(self) -> None:
        """Test filter_variants_by_clinvar with gene filtering."""
        # Test with a single gene.
        result = filter_variants_by_clinvar(
            genes="BRCA1",
            max_af_threshold=0.01,
            download_clinvar=False,
        )

        # Check that we get a Hail Table.
        assert isinstance(result, hl.Table)

        # Check that the result has ClinVar annotations.
        if result.count() > 0:
            assert "clinvar" in result.row, "Expected 'clinvar' annotation in result"

    def test_filter_variants_by_clinvar_with_multiple_genes(self) -> None:
        """Test filter_variants_by_clinvar with multiple genes."""
        # Test with multiple genes.
        result = filter_variants_by_clinvar(
            genes=["BRCA1", "BRCA2"],
            max_af_threshold=0.01,
            download_clinvar=False,
        )

        # Check that we get a Hail Table.
        assert isinstance(result, hl.Table)

        # Check that the result has ClinVar annotations.
        if result.count() > 0:
            assert "clinvar" in result.row, "Expected 'clinvar' annotation in result"

    def test_filter_variants_by_clinvar_af_thresholds(self) -> None:
        """Test filter_variants_by_clinvar with two AF thresholds."""
        # Test with both min and max AF thresholds.
        result = filter_variants_by_clinvar(
            intervals="chr17:43044295-43044305",
            max_af_threshold=0.001,
            min_af_threshold=0.0001,
            download_clinvar=False,
        )

        # Check that we get a Hail Table
        assert isinstance(result, hl.Table)

    def test_filter_variants_by_clinvar_custom_significances(self) -> None:
        """Test filter_variants_by_clinvar with custom ClinVar significances."""
        # Test with custom significances.
        result = filter_variants_by_clinvar(
            intervals="chr17:43044295-43044305",
            clinvar_significances=["pathogenic", "likely_pathogenic"],
            download_clinvar=False,
        )

        # Check that we get a Hail Table
        assert isinstance(result, hl.Table)

    def test_filter_variants_by_clinvar_missing_parameters(self) -> None:
        """Test that filter_variants_by_clinvar raises error with missing parameters."""
        with pytest.raises(ValueError, match="Must provide either genes, intervals, or contig/start/end parameters"):
            filter_variants_by_clinvar(
                max_af_threshold=0.01,
                download_clinvar=False,
            )

    def test_filter_variants_by_clinvar_with_contig_start_end(self) -> None:
        """Test filter_variants_by_clinvar with contig/start/end parameters."""
        result = filter_variants_by_clinvar(
            contig="chr17",
            start=43044295,
            end=43044305,
            max_af_threshold=0.01,
            download_clinvar=False,
        )

        # Check that we get a Hail Table
        assert isinstance(result, hl.Table)

        # Check that the result has ClinVar annotations
        if result.count() > 0:
            assert "clinvar" in result.row, "Expected 'clinvar' annotation in result"
