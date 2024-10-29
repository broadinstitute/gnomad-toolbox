"""Small functions to filter variants in gnomAD datasets, such as by allele frequency or by variant type."""

import hail as hl
from gnomad.resources.grch38.gnomad import POPS_TO_REMOVE_FOR_POPMAX, coverage
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.vep import (
    CSQ_CODING,
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)


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


def filter_by_csqs(ht: hl.Table, csqs: list[str]) -> hl.Table:
    """
    Filter variants by consequence.

    :param ht: Input Table.
    :param csqs: List of consequences.
    :return: Table with variants with the given consequences.
    """
    ht = ht.filter(
        hl.any(
            hl.map(
                lambda x: (x.consequence_terms.contains(csqs)),
                ht.vep.transcript_consequences,
            )
        )
    )

    return ht


def filter_by_gene_symbol(ht: hl.Table, gene: str) -> hl.Table:
    """
    Filter variants in a gene.

    :param ht: Input Table.
    :param gene: Gene symbol or.
    :return: Table with variants in the gene.
    """
    ht = filter_vep_transcript_csqs(
        ht,
        synonymous=False,
        mane_select=True,
        genes=[gene],
        match_by_gene_symbol=True,
    )

    return ht


def filter_to_coding_variants(ht: hl.Table) -> hl.Table:
    """
    Filter to coding variants.

    :param ht: Input Table.
    :return: Table with coding variants.
    """
    ht = filter_vep_transcript_csqs(
        ht,
        synonymous=False,
        canonical=True,
    )
    ht = get_most_severe_consequence_for_summary(ht)

    filter_expr = {}
    filter_expr["coding"] = hl.any(lambda csq: ht.most_severe_csq == csq, CSQ_CODING)

    ht = ht.filter(filter_expr["coding"])

    return ht


def filter_to_lof_variants(ht: hl.Table) -> hl.Table:
    """
    Filter to loss-of-function (LoF) variants.

    :param ht: Input Table.
    :return: Table with LoF variants.
    """
    ht = filter_vep_transcript_csqs(
        ht,
        lof=True,
        canonical=True,
    )
    ht = get_most_severe_consequence_for_summary(ht)

    filter_expr = {}
    filter_expr["lof"] = hl.any(lambda csq: ht.most_severe_csq == csq, CSQ_CODING)

    ht = ht.filter(filter_expr["lof"])

    return ht
