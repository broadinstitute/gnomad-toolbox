"""Functions for filtering the gnomAD sites HT frequency data."""

from typing import List

import hail as hl
from gnomad.utils.filtering import filter_arrays_by_meta

from gnomad_toolbox.filtering.variant import get_single_variant
from gnomad_toolbox.load_data import _get_gnomad_release


def get_callstats_for_multiple_ancestries(
    gen_ancs: List[str],
    **kwargs,
) -> hl.Table:
    """
    Extract callstats for specified ancestry groups.

    :param gen_ancs: List of genetic ancestry groups (e.g., 'afr', 'amr', 'asj', 'eas',
        'fin', 'nfe', 'oth', 'sas').
    :param kwargs: Keyword arguments to pass to _get_gnomad_release.
    :return: Table with callstats for the given ancestry groups and variant.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    # Format gen_ancs to lowercase and filter arrays by metadata.
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
        "filters",
        **{gen_anc: array_exprs["freq"][i] for i, gen_anc in enumerate(gen_ancs)},
    )

    # Select a subset of the globals.
    ht = ht.select_globals(
        "date",
        "version",
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )

    return ht


def get_callstats_for_single_ancestry(
    gen_anc: str,
    **kwargs,
) -> hl.Table:
    """
    Extract callstats for a specific ancestry group and single variant.

    :param gen_anc: Genetic ancestry group (e.g., 'afr', 'amr', 'asj', 'eas', 'fin',
        'nfe', 'oth', 'sas').
    :param kwargs: Keyword arguments to pass to _get_gnomad_release.
    :return: Table with callstats for the given ancestry group.
    """
    ht = get_callstats_for_multiple_ancestries([gen_anc], **kwargs)

    # Select a subset of the globals.
    ht = ht.select_globals(
        "date",
        "version",
        sample_count=ht.freq_meta_sample_count[0],
    )

    return ht


def get_single_variant_callstats_for_multiple_ancestries(
    gen_ancs: List[str],
    **kwargs,
) -> hl.Table:
    """
    Extract callstats for specified ancestry groups and a single variant.

    :param gen_ancs: List of genetic ancestry groups (e.g., 'afr', 'amr', 'asj', 'eas',
        'fin', 'nfe', 'oth', 'sas').
    :param kwargs: Keyword arguments to pass to get_single_variant.
    :return: Table with callstats for the given ancestry groups and variant.
    """
    ht = get_single_variant(**kwargs)

    return get_callstats_for_multiple_ancestries(gen_ancs, ht=ht)


def get_single_variant_callstats_for_single_ancestry(
    gen_anc: str,
    **kwargs,
) -> hl.Table:
    """
    Extract callstats for a specific ancestry group and single variant.

    :param gen_anc: Genetic ancestry group (e.g., 'afr', 'amr', 'asj', 'eas', 'fin',
        'nfe', 'oth', 'sas').
    :param kwargs: Keyword arguments to pass to get_single_variant.
    :return: Table with callstats for the given ancestry group.
    """
    ht = get_single_variant(**kwargs)

    return get_callstats_for_single_ancestry(gen_anc, ht=ht)
