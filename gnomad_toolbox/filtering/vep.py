"""Functions to filter gnomAD sites HT by VEP annotations."""

from functools import reduce

import hail as hl
from gnomad.utils.vep import CSQ_CODING, LOF_CSQ_SET


def filter_by_csqs(
    ht: hl.Table, csqs: list[str], pass_filters: bool = True
) -> hl.Table:
    """
    Filter variants by consequences.

    :param ht: Input Table.
    :param csqs: List of consequences to filter by. It can be specified as the
         categories on the browser: pLoF, Missense / Inframe indel, Synonymous, Other.
    :param pass_filters: Boolean if the variants pass the filters.
    :return: Table with variants with the specified consequences.
    """
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
