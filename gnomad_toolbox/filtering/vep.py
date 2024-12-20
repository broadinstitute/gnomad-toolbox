"""Functions to filter gnomAD sites HT by VEP annotations."""

import hail as hl
from gnomad.utils.vep import LOF_CSQ_SET, filter_vep_transcript_csqs_expr

from gnomad_toolbox.load_data import _get_gnomad_release


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

        pLoF

            - transcript_ablation
            - splice_acceptor_variant
            - splice_donor_variant
            - stop_gained
            - frameshift_variant

        Missense / Inframe indel

            - stop_lost
            - start_lost
            - inframe_insertion
            - inframe_deletion
            - missense_variant

        Synonymous

            - synonymous_variant

        `Other

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


# TODO: The following was in one of the notebooks, and I think we should add a wrapper
#  around this function to make it much simpler instead of using it in the notebook.

# Filter to LOFTEE high-confidence variants for certain genes

# In this example, we are filtering to variants in ASH1L that are LOFTEE high-confidence
# (with no flags) in the MANE select transcript.

# from gnomad.utils.vep import filter_vep_transcript_csqs
# ht = get_gnomad_release(data_type='exomes', version='4.1')
# ht = filter_vep_transcript_csqs(
#    ht,
#    synonymous=False,
#    mane_select=True,
#    genes=["ASH1L"],
#    match_by_gene_symbol=True,
#    additional_filtering_criteria=[lambda x: (x.lof == "HC") & hl.is_missing(x.lof_flags)],
# )
# ht.show()
# ht.count()
