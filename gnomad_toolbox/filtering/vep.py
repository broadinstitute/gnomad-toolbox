"""Functions to filter gnomAD sites HT by VEP annotations."""

from typing import List

import hail as hl
from gnomad.resources.grch37.reference_data import gencode as grch37_gencode
from gnomad.resources.grch38.reference_data import gencode as grch38_gencode
from gnomad.utils.vep import (
    LOF_CSQ_SET,
    filter_vep_transcript_csqs,
    filter_vep_transcript_csqs_expr,
)

from gnomad_toolbox.load_data import (
    _get_gnomad_release,
    get_coverage_for_variant,
    gnomad_session,
)


# TODO: Check these csq sets, the ones in the code don't match what is listed on the
#  browser. We should make sure they are consistent.
def filter_by_consequence_category(
    plof: bool = False,
    missense: bool = False,
    synonymous: bool = False,
    other: bool = False,
    pass_filters: bool = True,
    **kwargs,
) -> hl.Table:
    """
    Filter gnomAD variants based on VEP consequence.

    https://gnomad.broadinstitute.org/help/consequence-category-filter

    The [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) consequences included in each category are:

        pLoF:

            - transcript_ablation
            - splice_acceptor_variant
            - splice_donor_variant
            - stop_gained
            - frameshift_variant

        Missense / Inframe indel:

            - stop_lost
            - start_lost
            - inframe_insertion
            - inframe_deletion
            - missense_variant

        Synonymous:

            - synonymous_variant

        Other:

            - protein_altering_variant
            - incomplete_terminal_codon_variant
            - stop_retained_variant
            - coding_sequence_variant
            - mature_miRNA_variant
            - 5_prime_UTR_variant
            - 3_prime_UTR_variant
            - non_coding_transcript_exon_variant
            - non_coding_exon_variant
            - NMD_transcript_variant
            - non_coding_transcript_variant
            - nc_transcript_variant
            - downstream_gene_variant
            - TFBS_ablation
            - TFBS_amplification
            - TF_binding_site_variant
            - regulatory_region_ablation
            - regulatory_region_amplification
            - feature_elongation
            - regulatory_region_variant
            - feature_truncation
            - intergenic_variant
            - intron_variant
            - splice_region_variant
            - upstream_gene_variant

    :param plof: Whether to include pLoF variants.
    :param missense: Whether to include missense variants.
    :param synonymous: Whether to include synonymous variants.
    :param other: Whether to include other variants.
    :param pass_filters: Boolean if the variants pass the filters.
    :param kwargs: Arguments to pass to _get_gnomad_release.
    :return: Table with variants with the specified consequences.
    """
    if not any([plof, missense, synonymous, other]):
        raise ValueError(
            "At least one of plof, missense, synonymous, or other must be True."
        )

    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    lof_csqs = list(LOF_CSQ_SET)
    missense_csqs = ["missense_variant", "inframe_insertion", "inframe_deletion"]
    synonymous_csqs = ["synonymous_variant"]
    other_csqs = lof_csqs + missense_csqs + synonymous_csqs

    csqs = (
        (lof_csqs if plof else [])
        + (missense_csqs if missense else [])
        + (synonymous_csqs if synonymous else [])
    )

    filter_expr = None

    if csqs:
        filter_expr = filter_vep_transcript_csqs_expr(ht.vep, csqs=csqs)

    if other:
        other_expr = filter_vep_transcript_csqs_expr(
            ht.vep, csqs=other_csqs, keep_csqs=False
        )
        filter_expr = other_expr if filter_expr is None else (filter_expr | other_expr)

    if pass_filters:
        pass_expr = hl.len(ht.filters) == 0
        filter_expr = pass_expr if filter_expr is None else (filter_expr & pass_expr)

    return ht.filter(filter_expr)


def get_gene_intervals(gene_symbol: str, version: str) -> List[hl.utils.Interval]:
    """
    Get the genomic intervals for a given gene symbol.

    :param gene_symbol: Gene symbol.
    :param version: Dataset version.
    :return: List of intervals for the specified gene.
    """
    gen_ht = grch37_gencode.ht() if version.startswith("2.") else grch38_gencode.ht()
    intervals = (
        gen_ht.filter(
            (gen_ht.feature == "gene")
            & (gen_ht.gene_name.lower() == gene_symbol.lower())
        )
        .select()
        .collect()
    )
    if not intervals:
        raise ValueError(f"No interval found for gene: {gene_symbol}")
    return [row["interval"] for row in intervals]


def filter_hc_variants(var_ht: hl.Table, gene_symbol: str, version: str) -> hl.Table:
    """
    Filter variants to high-confidence LOFTEE variants with optional transcript selection.

    :param var_ht: Variants Table.
    :param gene_symbol: Gene symbol.
    :param version: Dataset version.
    :return: Filtered variants Hail Table.
    """
    return filter_vep_transcript_csqs(
        var_ht,
        synonymous=False,
        mane_select=version.startswith("4."),
        genes=[gene_symbol.upper()],
        match_by_gene_symbol=True,
        additional_filtering_criteria=[
            lambda x: (x.lof == "HC")
            & (hl.is_missing(x.lof_flags) | (x.lof_flags == ""))
        ],
    )


def filter_to_plofs(
    gene_symbol: str, select_fields: bool = False, **kwargs
) -> hl.Table:
    """
    Filter to observed pLoF variants used for gene constraint metrics.

    .. note::

                pLOF variants meets the following requirements:
                - High-confidence LOFTEE variants (without any flags),
                - Only variants in the MANE Select transcript,
                - PASS variants that are SNVs with MAF ≤ 0.1%,
                - Exome median depth ≥ 30 (# TODO: This is changing in v4 constraint?)

    :param gene_symbol: Gene symbol.
    :param select_fields: Whether to limit the output to specific fields.
    :return: Table with pLoF variants.
    """
    var_version = kwargs.pop("version", gnomad_session.version)
    var_ht = _get_gnomad_release(dataset="variant", version=var_version, **kwargs)
    cov_ht = get_coverage_for_variant(var_version, **kwargs)

    # Get gene intervals and filter tables
    intervals = get_gene_intervals(gene_symbol, var_version)
    var_ht = hl.filter_intervals(var_ht, intervals)
    cov_ht = hl.filter_intervals(cov_ht, intervals)

    # Filter to high-confidence LOFTEE variants
    var_ht = filter_hc_variants(var_ht, gene_symbol, var_version)

    # Version-specific expressions
    if var_version.startswith("2."):
        allele_type_expr = var_ht.allele_type
        cov_cut_expr = cov_ht[var_ht.locus].median
    else:
        allele_type_expr = var_ht.allele_info.allele_type
        cov_cut_expr = cov_ht[var_ht.locus].median_approx

    # Apply final filters
    var_ht = var_ht.filter(
        (hl.len(var_ht.filters) == 0)
        & (allele_type_expr == "snv")
        & (var_ht.freq[0].AF <= 0.001)
        & (cov_cut_expr >= 30)
    )

    # Select specific fields if requested
    if select_fields:
        var_ht = var_ht.select(
            freq=var_ht.freq[0],
            csq=var_ht.vep.transcript_consequences[0].consequence_terms,
            coverage=cov_ht[var_ht.locus],
        )

    return var_ht
