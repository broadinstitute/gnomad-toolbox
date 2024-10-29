"""Extract callstats from 'freq' of gnomAD HTs."""

from typing import List, Optional

import hail as hl
from gnomad.resources.grch38.gnomad import POPS_TO_REMOVE_FOR_POPMAX
from gnomad.utils.filtering import filter_arrays_by_meta


def extract_callstats_for_1anc_1variant(
    ht: hl.Table, gen_anc: str, contig: str, position: int, alleles: Optional[List[str]]
) -> hl.Table:
    """
    Extract callstats for a single genetic ancestry group and a single variant.

    :param ht: Input Table.
    :param group: Ancestry Group (e.g., 'afr', 'amr', 'asj', 'eas', 'fin', 'nfe',
    'oth', 'sas').
    :param contig: Chromosome.
    :param position: Position.
    :param alleles: List of alleles.
    :return: Table with callstats for the given group.
    """
    # Filter to the variant of interest
    ht = ht.filter(
        (ht.locus.contig == contig)
        & (ht.locus.position == position)
        & (ht.alleles == alleles)
    )

    # Format the gen_anc to lowercase if it's fed in as uppercase
    gen_anc = gen_anc.lower()
    items_to_filter = {"gen_anc": [gen_anc], "group": ["adj"]}
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
    ht = ht.select(**{gen_anc: array_exprs["freq"]})
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )
    return ht
