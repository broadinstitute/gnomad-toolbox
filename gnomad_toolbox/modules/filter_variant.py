"""Small functions to filter variants in gnomAD datasets, such as by allele frequency or by variant type."""

import hail as hl


def get_variant_count(
    ht: hl.Table,
    afs: list[float] = [0.01, 0.001],
    singletons: bool = False,
    doubletons: bool = False,
) -> dict:
    """
    Count variants with frequency <1%, <0.1%, and singletons (AC == 1).

    .. note:: This function works for gnomAD exomes and genomes datasets, not yet for
              gnomAD joint dataset, since the HT schema is slightly different.

    :param ht: Input Table.
    :param afs: List of allele frequencies cutoffs.
    :param singletons: Include singletons.
    :param doubletons: Include doubletons.
    :return: Dictionary with counts.
    """
    counts = {}

    # Filter to PASS variants.
    ht = ht.filter(hl.len(ht.filters) == 0)
    if singletons:
        n_singletons = ht.aggregate(hl.agg.count_where(ht.freq[0].AC == 1))
        counts["number of singletons"] = n_singletons
    if doubletons:
        n_doubletons = ht.aggregate(hl.agg.count_where(ht.freq[0].AC == 2))
        counts["number of doubletons"] = n_doubletons

    for af in afs:
        n_variants = ht.aggregate(hl.agg.count_where(ht.freq[0].AF < af))
        counts[f"number of variants with AF < {af}"] = n_variants

    # Count variants with frequency <1%, <0.1%, and singletons (AC == 1).
    return counts
