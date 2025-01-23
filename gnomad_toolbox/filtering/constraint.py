"""Functions to filter gnomAD sites HT by constraint metrics."""

import hail as hl

from gnomad_toolbox.filtering.vep import (
    filter_to_high_confidence_loftee,
    get_gene_intervals,
)
from gnomad_toolbox.load_data import (
    CONSTRAINT_DATA,
    _get_dataset,
    get_compatible_dataset_versions,
)


def get_observed_plofs_for_gene_constraint(
    gene_symbol: str,
    version: str = None,
    variant_ht: hl.Table = None,
    coverage_ht: hl.Table = None,
) -> hl.Table:
    """
    Filter to observed pLoF variants used for gene constraint metrics.

    The pLOF variant count displayed on the browser meets the following requirements:

        - PASS variant QC
        - SNV
        - Allele frequency ≤ 0.1%
        - High-confidence LOFTEE in the MANE Select or Canonical transcript
        - ≥ a specified coverage threshold (depends on the version)

    :param gene_symbol: Gene symbol.
    :param version: Optional gnomAD dataset version. If not provided, uses the gnomAD
        session version.
    :param variant_ht: Optional Hail Table with variants. If not provided, uses the
        exome variant Table for the gnomAD session version.
    :param coverage_ht: Optional Hail Table with coverage data. If not provided, uses
        the exome coverage Table for the gnomAD session version.
    :return: Table with pLoF variants.
    """
    if variant_ht is not None and coverage_ht is None:
        raise ValueError("Variant Hail Table provided without coverage Hail Table.")

    if coverage_ht is not None and variant_ht is None:
        raise ValueError("Coverage Hail Table provided without variant Hail Table.")

    # Load the variant exomes Hail Table if not provided.
    variant_ht = _get_dataset(
        dataset="variant",
        ht=variant_ht,
        data_type="exomes",
        version=version,
    )

    # Determine the coverage version compatible with the variant version.
    coverage_version = get_compatible_dataset_versions("coverage", version, "exomes")

    # Load the coverage Hail Table if not provided.
    coverage_ht = _get_dataset(
        dataset="coverage",
        ht=coverage_ht,
        data_type="exomes",
        version=coverage_version,
    )

    # Get gene intervals and filter tables.
    gencode_version = get_compatible_dataset_versions("gencode", version)
    intervals = get_gene_intervals(gene_symbol, gencode_version=gencode_version)
    variant_ht = hl.filter_intervals(variant_ht, intervals)
    coverage_ht = hl.filter_intervals(coverage_ht, intervals)

    # Determine constraint filters.
    constraint_version = get_compatible_dataset_versions("constraint", version)
    constraint_info = CONSTRAINT_DATA[constraint_version]
    cov_field = constraint_info["exome_coverage_field"]
    cov_cutoff = constraint_info["exome_coverage_cutoff"]
    af_cutoff = constraint_info["af_cutoff"]

    # Annotate the exome coverage.
    variant_ht = variant_ht.annotate(
        exome_coverage=coverage_ht[variant_ht.locus][cov_field]
    )

    # Apply constraint filters.
    variant_ht = variant_ht.filter(
        (hl.len(variant_ht.filters) == 0)
        & (hl.is_snp(variant_ht.alleles[0], variant_ht.alleles[1]))
        & (variant_ht.freq[0].AF <= af_cutoff)
        & (variant_ht.exome_coverage >= cov_cutoff)
    )

    # Filter to high-confidence LOFTEE variants.
    variant_ht = filter_to_high_confidence_loftee(
        gene_symbol=gene_symbol,
        ht=variant_ht,
        canonical_only=True,
    )

    return variant_ht


# noqa: D104, D100
