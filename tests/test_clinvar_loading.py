"""Tests for ClinVar loading functionality."""

import pytest
import hail as hl

from gnomad_toolbox.load_reference_data import _load_clinvar_resource


class TestClinVarLoading:
    """Test class for ClinVar loading functionality."""

    def test_load_clinvar_resource_invalid_build(self) -> None:
        """Test that invalid reference genome build raises ValueError."""
        with pytest.raises(ValueError, match="Unsupported reference genome build"):
            _load_clinvar_resource("invalid_build")


    @pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
    def test_clinvar_resource_basic_properties(self, build: str) -> None:
        """Test basic properties of ClinVar resources for both builds."""
        clinvar_ht = _load_clinvar_resource(build)

        # Basic type checks
        assert isinstance(clinvar_ht, hl.Table)
        assert clinvar_ht.count() > 0

        # Check that the table has the required key fields
        assert "locus" in clinvar_ht.key
        assert "alleles" in clinvar_ht.key

        # Check that the table has the required row fields
        clinvar_ht = clinvar_ht.transmute(**clinvar_ht.info)
        required_fields = ["CLNSIG", "CLNREVSTAT", "CLNDN", "CLNDISDB"]
        for field in required_fields:
            assert field in clinvar_ht.row, f"Missing required field: {field}"
