"""Set of general functions for gnomAD analysis."""

from typing import Dict, List, Optional, Tuple, Union

import hail as hl
from gnomad.assessment.summary_stats import freq_bin_expr

from gnomad_toolbox.load_data import _get_gnomad_release


def get_variant_count_by_freq_bin(
    af_cutoffs: List[float] = [0.001, 0.01],
    singletons: bool = False,
    doubletons: bool = False,
    pass_only: bool = True,
    **kwargs,
) -> Dict[str, int]:
    """
    Count variants by frequency bin.

    By default, this function counts PASS variants that are AC0, AF < 0.01%, and
    AF 0.01% - 0.1%.

    The function can also include counts of singletons and doubletons, with or
    without passing filters.

    .. note::

        This function works for gnomAD exomes and genomes data types, not yet for gnomAD
        joint data type, since the HT schema is slightly different.

    :param af_cutoffs: List of allele frequencies cutoffs.
    :param singletons: Include singletons.
    :param doubletons: Include doubletons.
    :param pass_only: Include only PASS variants.
    :param kwargs: Keyword arguments to pass to _get_gnomad_release. Includes
        'ht', 'data_type', and 'version'.
    :return: Dictionary with counts.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    # Filter to PASS variants.
    if pass_only:
        ht = ht.filter(hl.len(ht.filters) == 0)

    # Initialize allele count cutoffs with AC0.
    ac_cutoffs = [(0, "AC0")]

    if singletons:
        ac_cutoffs.append((1, "singletons"))

    if doubletons:
        ac_cutoffs.append((2, "doubletons"))

    freq_expr = freq_bin_expr(
        ht.freq, ac_cutoffs=ac_cutoffs, af_cutoffs=af_cutoffs, upper_af=None
    )

    return ht.aggregate(hl.agg.counter(freq_expr))
