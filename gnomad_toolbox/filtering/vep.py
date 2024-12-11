"""Functions to filter gnomAD sites HT by VEP annotations."""

from functools import reduce

import hail as hl
from gnomad.utils.vep import CSQ_CODING, LOF_CSQ_SET

from gnomad_toolbox.load_data import _get_gnomad_release

# TODO: I haven't looked over this function yet. Is there anything in gnomad_methods
#  that could be used here? If not, is there anything here that should be moved to
#  gnomad_methods?


def filter_by_csqs(
    csqs: list[str],
    pass_filters: bool = True,
    **kwargs,
) -> hl.Table:
    """
    Filter variants by VEP transcript consequences.

    :param csqs: List of consequences to filter by. It can be specified as the
         categories on the browser: pLoF, Missense / Inframe indel, Synonymous, Other.
    :param pass_filters: Boolean if the variants pass the filters.
    :param kwargs: Arguments to pass to _get_gnomad_release.
    :return: Table with variants with the specified consequences.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    missense_inframe = ["missense_variant", "inframe_insertion", "inframe_deletion"]

    filter_expr = []
    if "lof" in csqs:
        filter_expr.append(
            hl.literal(LOF_CSQ_SET).contains(ht.vep.most_severe_consequence)
        )

    if "synonymous" in csqs:
        filter_expr.append(ht.vep.most_severe_consequence == "synonymous_variant")

    if "missense" in csqs:
        filter_expr.append(
            hl.literal(missense_inframe).contains(ht.vep.most_severe_consequence)
        )

    if "other" in csqs:
        excluded_csqs = hl.literal(
            list(LOF_CSQ_SET) + missense_inframe + ["synonymous_variant"]
        )
        filter_expr.append(~excluded_csqs.contains(ht.vep.most_severe_consequence))

    if "coding" in csqs:
        filter_expr.append(
            hl.literal(CSQ_CODING).contains(ht.vep.most_severe_consequence)
        )

    if len(filter_expr) == 0:
        raise ValueError(
            "No valid consequence specified. Choose from 'lof', 'synonymous', 'missense', 'other'."
        )

    # Combine filter expressions with logical OR
    if len(filter_expr) == 1:
        combined_filter = filter_expr[0]
    else:
        combined_filter = reduce(lambda acc, expr: acc | expr, filter_expr)

    ht = ht.filter(combined_filter)

    if pass_filters:
        ht = ht.filter(hl.len(ht.filters) == 0)

    return ht


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
