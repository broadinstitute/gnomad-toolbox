"""Functions for filtering the gnomAD sites HT frequency data."""

from typing import List, Union

import hail as hl
from gnomad.utils.filtering import filter_arrays_by_meta

from gnomad_toolbox.filtering.variant import get_single_variant
from gnomad_toolbox.load_data import _get_gnomad_release


def get_ancestry_callstats(
    gen_ancs: Union[str, List[str]],
    **kwargs,
) -> hl.Table:
    """
    Extract callstats for specified ancestry group(s).

    :param gen_ancs: Genetic ancestry group(s) (e.g., 'afr', 'amr', 'asj', 'eas',
        'fin', 'nfe', 'oth', 'sas'). Can be a single ancestry group or a list of
        ancestry groups.
    :param kwargs: Keyword arguments to pass to _get_gnomad_release.
    :return: Table with callstats for the given ancestry groups and variant.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    # Check if gen_ancs is a single ancestry group.
    one_anc = isinstance(gen_ancs, str)

    if one_anc:
        gen_ancs = [gen_ancs]

    # Format gen_ancs to lowercase and filter arrays by metadata.
    gen_ancs = [gen_anc.lower() for gen_anc in gen_ancs]
    gen_anc_label = (
        "gen_anc" if any(["gen_anc" in m for m in hl.eval(ht.freq_meta)]) else "pop"
    )
    items_to_filter = {gen_anc_label: gen_ancs, "group": ["adj"]}
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
        **{
            m[gen_anc_label]: array_exprs["freq"][i]
            for i, m in enumerate(hl.eval(freq_meta))
        },
    )

    # Select a subset of the globals.
    sample_count = array_exprs["freq_meta_sample_count"]
    if one_anc:
        sample_count = sample_count[0]
    else:
        sample_count = hl.struct(
            **{
                m[gen_anc_label]: sample_count[i]
                for i, m in enumerate(hl.eval(freq_meta))
            }
        )
    ht = ht.select_globals("date", "version", sample_count=sample_count)

    return ht


def get_single_variant_ancestry_callstats(
    gen_ancs: Union[str, List[str]],
    **kwargs,
) -> hl.Table:
    """
    Extract callstats for specified ancestry group(s) and a single variant.

    :param gen_ancs: Genetic ancestry group(s) (e.g., 'afr', 'amr', 'asj', 'eas',
        'fin', 'nfe', 'oth', 'sas'). Can be a single ancestry group or a list of
        ancestry groups.
    :param kwargs: Keyword arguments to pass to get_single_variant.
    :return: Table with callstats for the given ancestry groups and variant.
    """
    ht = get_single_variant(**kwargs)

    return get_ancestry_callstats(gen_ancs, ht=ht)
