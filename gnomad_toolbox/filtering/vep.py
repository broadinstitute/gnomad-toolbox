"""Functions to filter gnomAD sites HT by VEP annotations."""

import hail as hl
from gnomad.resources.grch37.reference_data import gencode as grch37_gencode
from gnomad.resources.grch38.reference_data import gencode as grch38_gencode
from gnomad.utils.vep import (
    LOF_CSQ_SET,
    filter_vep_transcript_csqs,
    filter_vep_transcript_csqs_expr,
)

from gnomad_toolbox.load_data import _get_gnomad_release, gnomad_session


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


def filter_to_plofs(
    gene_symbol: str, select_fields: bool = False, **kwargs
) -> hl.Table:
    """
        Filter to observed pLoF variants that we used to calculate the gene constraint metrics.

        .. note::

                pLOF variants meets the following requirements:
                - High-confidence LOFTEE variants (without any flags),
                - Only variants in the MANE Select transcript,
                - PASS variants that are SNVs with MAF ≤ 0.1%,
                - Exome median depth ≥ 30 (# TODO: This is changing in v4 constraint?)

    **Note: this number should match the number of observed pLOF SNVs on the gene page of gnomAD Browser.**

        :param gene_symbol: Gene symbol.
        :param select_fields: Boolean if the output should be limited to specific fields.
        :return: Table with pLoF variants.
    """
    var_version = kwargs.pop("version", gnomad_session.version)
    var_ht = _get_gnomad_release(dataset="variant", version=var_version, **kwargs)

    # Determine the version of the coverage table
    if var_version.startswith("4."):
        cov_ht = _get_gnomad_release(dataset="coverage", version="4.0", **kwargs)
    elif var_version.startswith("3."):
        cov_ht = _get_gnomad_release(dataset="coverage", version="3.0.1", **kwargs)
    elif var_version.startswith("2."):
        cov_ht = _get_gnomad_release(dataset="coverage", version="2.1", **kwargs)
    else:
        raise ValueError(
            f"Unrecognized version: '{var_version}'. Please specify a valid version."
        )

    # Get the gene interval from gen_ht
    gen_ht = (
        grch37_gencode.ht() if var_version.startswith("2.") else grch38_gencode.ht()
    )
    interval = (
        gen_ht.filter(
            (gen_ht.feature == "gene")
            & (gen_ht.gene_name.lower() == gene_symbol.lower())
        )
        .select()
        .collect()
    )

    if not interval:
        raise ValueError(f"No interval found for gene: {gene_symbol}")

    # Convert to a list of intervals
    interval = [row["interval"] for row in interval]
    var_ht = hl.filter_intervals(var_ht, interval)
    cov_ht = hl.filter_intervals(cov_ht, interval)

    # Filter to high-confidence LOFTEE variants
    var_ht = filter_vep_transcript_csqs(
        var_ht,
        synonymous=False,
        mane_select=True if var_version.startswith("4.") else False,
        # TODO: When this function is applied to DRD2 gene in v4.1, it will get 7 pLoF
        #  variants instead of 8 on the browser and the 4.1 constraint table,
        #  because one of them is not in mane select, nor in canonical transcript.
        genes=[gene_symbol.upper()],
        match_by_gene_symbol=True,
        additional_filtering_criteria=[
            lambda x: (x.lof == "HC")
            & (hl.is_missing(x.lof_flags) | (x.lof_flags == ""))
        ],
    )

    if var_version.startswith("2."):
        allele_type_expr = var_ht.allele_type
        cov_cut_expr = cov_ht[var_ht.locus].median
    else:
        allele_type_expr = var_ht.allele_info.allele_type
        cov_cut_expr = cov_ht[var_ht.locus].median_approx

    var_ht = var_ht.filter(
        (hl.len(var_ht.filters) == 0)
        & (allele_type_expr == "snv")
        & (var_ht.freq[0].AF <= 0.001)
        & (cov_cut_expr >= 30)
    )

    if select_fields:
        var_ht = var_ht.select(
            freq=var_ht.freq[0],
            csq=var_ht.vep.transcript_consequences[0].consequence_terms,
            coverage=cov_ht[var_ht.locus],
        )

    return var_ht
