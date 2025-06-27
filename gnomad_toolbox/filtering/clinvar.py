"""Functions to filter gnomAD variants and join with ClinVar data."""

from typing import List, Optional, Union

import hail as hl
from gnomad.utils.filtering import filter_by_intervals
from gnomad.utils.reference_genome import get_reference_genome

from gnomad_toolbox.filtering.variant import filter_by_gene_symbol
from gnomad_toolbox.load_data import _get_dataset
from gnomad_toolbox.load_reference_data import (
    _import_clinvar_vcf,
    _load_clinvar_resource,
)


def filter_variants_by_clinvar(
    max_af_threshold: float = 0.01,
    min_af_threshold: Optional[float] = None,
    download_clinvar: bool = False,
    genes: Optional[Union[str, List[str]]] = None,
    intervals: Optional[Union[str, List[str], hl.Interval, List[hl.Interval]]] = None,
    contig: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    clinvar_significances: List[str] = [
        "pathogenic",
        "likely_pathogenic",
        "uncertain_significance",
    ],
    clinvar_download_path: Optional[str] = None,
    clinvar_output_path: Optional[str] = None,
    **kwargs,
) -> hl.Table:
    """
    Filter gnomAD variants by gene(s) or interval(s), join with ClinVar data, and filter by AF threshold.

    .. note::

        One of `genes`, `interval`, or all of `contig`, `start`, and `end` must be
        provided. The function will filter variants based on the provided criteria,
        then join with ClinVar data and filter to retain only variants with the
        specified ClinVar clinical significances.

    :param max_af_threshold: Maximum allele frequency threshold. Default is 0.01.
    :param min_af_threshold: Minimum allele frequency threshold. If provided, variants
        with allele frequencies below this threshold will be excluded.
    :param download_clinvar: If True, download latest ClinVar VCF. If False, use gnomAD ClinVar resource from gnomad_methods.
        Default is False.
    :param genes: Single gene name or list of gene names to filter by. If provided,
        `intervals`, `contig`, `start`, and `end` are ignored.
    :param intervals: Single interval string (e.g., "chr1:1000-2000"), list of interval
        strings, Hail interval object, or list of Hail interval objects.
        If provided, `genes`, `contig`, `start`, and `end` are ignored.
    :param contig: Chromosome/contig name. Required if using `start` and `end`.
    :param start: Start position. Required if using `contig` and `end`.
    :param end: End position. Required if using `contig` and `start`.
    :param clinvar_significances: List of ClinVar clinical significances to retain.
        Default is ["Pathogenic", "Likely_pathogenic", "Uncertain_significance"].
    :param clinvar_download_path: Path to download ClinVar VCF. Required if `download_clinvar` is True.
        Default is None.
    :param clinvar_output_path: Path to store ClinVar Hail Table. Required if `download_clinvar` is True.
        Default is None.
    :param kwargs: Additional arguments to pass to `_get_dataset`.
    :return: Hail Table containing filtered variants with ClinVar annotations.
    """
    # Load gnomAD data.
    ht = _get_dataset(dataset="variant", **kwargs)
    build = get_reference_genome(ht.locus).name

    # Apply gene/interval filtering.
    if genes is not None:
        if isinstance(genes, str):
            genes = [genes]
        ht = filter_by_gene_symbol(genes[0], **kwargs)
        for gene in genes[1:]:
            ht = ht.union(filter_by_gene_symbol(gene, **kwargs))
    elif intervals is not None:
        ht = filter_by_intervals(ht, intervals)
    elif all(x is not None for x in [contig, start, end]):
        interval = f"{contig}:{start}-{end}"
        ht = filter_by_intervals(ht, interval)
    else:
        raise ValueError(
            "Must provide either genes, intervals, or contig/start/end parameters."
        )

    # Filter by AF thresholds.
    af_expr = ht.freq[0].AF if "freq" in ht.row else ht.af[0].AF
    if min_af_threshold is not None:
        ht = ht.filter(af_expr > min_af_threshold)
    ht = ht.filter(af_expr <= max_af_threshold)

    if download_clinvar:
        # Import ClinVar data.
        clinvar_ht = _import_clinvar_vcf(
            build=build,
            clinvar_download_path=clinvar_download_path,
            output_path=clinvar_output_path,
        )
    else:
        # Load ClinVar data from gnomad_methods.
        clinvar_ht = _load_clinvar_resource(build=build)
        clinvar_ht = clinvar_ht.transmute(**clinvar_ht.info)

    # Join ClinVar data with gnomAD data.
    ht = ht.annotate(
        clinvar=clinvar_ht[ht.key].select("CLNSIG", "CLNREVSTAT", "CLNDN", "CLNDISDB")
    )

    # Filter to specified ClinVar significances (case-insensitive).
    clinvar_filter = hl.any(
        lambda sig: hl.literal(clinvar_significances).contains(hl.str(sig).lower()),
        ht.clinvar.CLNSIG,
    )
    ht = ht.filter(clinvar_filter)

    return ht
