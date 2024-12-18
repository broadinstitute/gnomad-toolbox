"""Set of general functions for gnomAD analysis."""

from typing import Dict, List, Optional, Tuple, Union

import hail as hl

from gnomad_toolbox.load_data import _get_gnomad_release


# TODO: Modify this function in gnomad_methods.
def freq_bin_expr(
    freq_expr: Union[hl.expr.StructExpression, hl.expr.ArrayExpression],
    index: int = 0,
    ac_cutoffs: Optional[List[Union[int, Tuple[int, str]]]] = [
        (0, "AC0"),
        (1, "singleton"),
        (2, "doubleton"),
    ],
    af_cutoffs: Optional[List[Union[float, Tuple[float, str]]]] = [
        (1e-4, "0.01%"),
        (1e-3, "0.1%"),
        (1e-2, "1%"),
        (1e-1, "10%"),
    ],
    upper_af: Optional[Union[float, Tuple[float, str]]] = (0.95, "95%"),
) -> hl.expr.StringExpression:
    """
    Return frequency string annotations based on input AC or AF.

    .. note::

        - Default index is 0 because function assumes freq_expr was calculated with
          `annotate_freq`.
        - Frequency index 0 from `annotate_freq` is frequency for all pops calculated
          on adj genotypes only.

    :param freq_expr: Array of structs containing frequency information.
    :param index: Which index of freq_expr to use for annotation. Default is 0.
    :param ac_cutoffs: List of AC cutoffs to use for binning. 
    :param af_cutoffs: List of AF cutoffs to use for binning.
    :param upper_af: Upper AF cutoff to use for binning.
    :return: StringExpression containing bin name based on input AC or AF.
    """
    if isinstance(freq_expr, hl.expr.ArrayExpression):
        freq_expr = freq_expr[index]

    if ac_cutoffs and isinstance(ac_cutoffs[0], int):
        ac_cutoffs = [(c, f"AC{c}") for c in ac_cutoffs]

    if af_cutoffs and isinstance(af_cutoffs[0], float):
        af_cutoffs = [(f, f"{f*100}%") for f in af_cutoffs]

    if isinstance(upper_af, float):
        upper_af = (upper_af, f"{upper_af*100}%")

    freq_bin_expr = hl.case().when(hl.is_missing(freq_expr.AC), "Missing")
    prev_af = None
    for ac, name in sorted(ac_cutoffs):
        freq_bin_expr = freq_bin_expr.when(freq_expr.AC == ac, name)
        prev_af = name

    for af, name in sorted(af_cutoffs):
        prev_af = "<" if prev_af is None else f"{prev_af} - "
        freq_bin_expr = freq_bin_expr.when(freq_expr.AF < af, f"{prev_af}{name}")
        prev_af = name

    if upper_af:
        freq_bin_expr = freq_bin_expr.when(
            freq_expr.AF > upper_af[0], f">{upper_af[1]}"
        )
        default_af = "<" if prev_af is None else f"{prev_af} - "
        default_af = f"{default_af}{upper_af[1]}"
    else:
        default_af = f">{prev_af}"

    return freq_bin_expr.default(default_af)


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
