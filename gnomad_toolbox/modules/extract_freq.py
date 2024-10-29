"""Extract callstats from 'freq' of gnomAD HTs."""

from typing import List, Optional

import hail as hl
from gnomad.resources.grch38.gnomad import POPS_TO_REMOVE_FOR_POPMAX
from gnomad.utils.filtering import filter_arrays_by_meta


def extract_callstats_for_1anc_1variant(
    ht: hl.Table, gen_anc: str, contig: str, position: int, alleles: List[str]
) -> hl.Table:
    """
    Extract callstats for a specific ancestry group and single variant.

    :param ht: Input Hail Table with variant data.
    :param gen_anc: Genetic ancestry group (e.g., 'afr', 'nfe').
    :param contig: Chromosome of the variant.
    :param position: Variant position.
    :param alleles: List of alleles for the variant (e.g., ['A', 'T']).
    :return: Filtered Table with callstats for the specified group.
    """
    # Filter to the variant of interest
    ht = ht.filter(
        (ht.locus.contig == contig)
        & (ht.locus.position == position)
        & (ht.alleles == alleles)
    )

    # Check if the variant exists
    if ht.count() == 0:
        hl.utils.warning(
            f"No variant found at {contig}:{position} with alleles {alleles}"
        )

    # Format gen_anc to lowercase and filter arrays by metadata
    items_to_filter = {"gen_anc": [gen_anc.lower()], "group": ["adj"]}
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            **{a: ht[a] for a in ["freq"]},
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=items_to_filter,
        keep=True,
        combine_operator="and",
        exact_match=True,
    )
    # Select frequency for ancestry group
    ht = ht.select(
        **{
            gen_anc: array_exprs["freq"][i]
            for i, gen_anc in enumerate([gen_anc.lower()])
        }
    )
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )
    return ht


def extract_callstats_for_multiple_ancs(
    ht: hl.Table,
    gen_ancs: List[str],
) -> hl.Table:
    """
    Extract callstats for multiple genetic ancestry groups.

    :param ht: Input Table.
    :param gen_ancs: List of Ancestry Groups (e.g., 'afr', 'amr', 'asj', 'eas', 'fin', 'nfe',
    'oth', 'sas').
    :return: Table with callstats for the given groups.
    """
    # Format the gen_ancs to lowercase if they're fed in as uppercase
    gen_ancs = [gen_anc.lower() for gen_anc in gen_ancs]
    items_to_filter = {"gen_anc": gen_ancs, "group": ["adj"]}
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            **{a: ht[a] for a in ["freq"]},
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=items_to_filter,
        keep=True,
        combine_operator="and",
        exact_match=True,
    )
    ht = ht.select(
        **{gen_anc: array_exprs["freq"][i] for i, gen_anc in enumerate(gen_ancs)}
    )
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )
    return ht
