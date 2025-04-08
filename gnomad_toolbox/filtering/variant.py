"""Functions to filter the gnomAD sites HT to a specific set of variants."""

from typing import Optional, Union

import hail as hl
from gnomad.utils.filtering import filter_by_intervals as interval_filter
from gnomad.utils.filtering import filter_gencode_ht
from gnomad.utils.parse import parse_variant
from gnomad.utils.reference_genome import get_reference_genome

from gnomad_toolbox.load_data import _get_dataset


def get_single_variant(
    variant: Optional[str] = None,
    contig: Optional[str] = None,
    position: Optional[int] = None,
    ref: Optional[str] = None,
    alt: Optional[str] = None,
    **kwargs,
) -> hl.Table:
    """
    Get a single variant from the gnomAD HT.

    .. note::

        One of `variant` or all of `contig`, `position`, `ref`, and `alt` must be
        provided. If `variant` is provided, `contig`, `position`, `ref`, and `alt` are
        ignored.

    :param variant: Variant string in the format "chr12-235245-A-C" or
        "chr12:235245:A:C". If provided, `contig`, `position`, `ref`, and `alt` are
        ignored.
    :param contig: Chromosome of the variant. Required if `variant` is not provided.
    :param position: Variant position. Required if `variant` is not provided.
    :param ref: Reference allele. Required if `variant` is not provided.
    :param alt: Alternate allele. Required if `variant` is not provided.
    :param kwargs: Additional arguments to pass to `_get_dataset`.
    :return: Table with the single variant.
    """
    if not variant and not all([contig, position, ref, alt]):
        raise ValueError(
            "Either `variant` must be provided or all of `contig`, `position`, `ref`, "
            "and `alt`."
        )

    # Load the Hail Table if not provided
    ht = _get_dataset(dataset="variant", **kwargs)

    # Determine the reference genome build for the ht.
    build = get_reference_genome(ht.locus).name

    # Filter to the Locus of the variant of interest.
    variant = parse_variant(variant, contig, position, ref, alt, build)
    ht = hl.filter_intervals(
        ht, [hl.interval(variant.locus, variant.locus, includes_end=True)]
    )

    # Filter to the variant of interest.
    ht = ht.filter(ht.alleles == variant.alleles)

    # Check if the variant exists.
    if ht.count() == 0:
        hl.utils.warning(
            f"No variant found at {hl.eval(variant.locus)} with alleles "
            f"{hl.eval(variant.alleles)}"
        )

    return ht


def filter_by_intervals(
    intervals: Union[str, list[str]],
    padding_bp: int = 0,
    **kwargs,
) -> hl.Table:
    """
    Filter variants by interval(s).

    :param intervals: Interval string or list of interval strings. The interval string
        format has to be "contig:start-end", e.g.,"1:1000-2000" (GRCh37) or
        "chr1:1000-2000" (GRCh38).
    :param padding_bp: Number of base pairs to pad the intervals. Default is 0bp.
    :param kwargs: Arguments to pass to `_get_dataset`.
    :return: Table with variants in the interval(s).
    """
    # Load the Hail Table if not provided.
    ht = _get_dataset(dataset="variant", **kwargs)

    return interval_filter(
        ht,
        intervals,
        padding_bp=padding_bp,
        reference_genome=get_reference_genome(ht.locus).name,
    )


def filter_by_gene_symbol(gene: str, exon_padding_bp: int = 75, **kwargs) -> hl.Table:
    """
    Filter variants by gene symbol.

    .. note::

           This function is to match the number of variants that you will get in the
           gnomAD browser when you search for a gene symbol. The gnomAD browser
           filters to only variants located in or within 75 base pairs of CDS or
           non-coding exons of a gene.

    :param gene: Gencode gene symbol.
    :param exon_padding_bp: Number of base pairs to pad the intervals. Default is 75bp.
    :param kwargs: Arguments to pass to `_get_dataset`.
    :return: Table with variants in the specified gene.
    """
    # Load the Hail Table if not provided.
    ht = _get_dataset(dataset="variant", **kwargs)

    # The gnomAD browser will display variants in CDS regions if present, otherwise UTR,
    # and finally exons.
    feature_order = ["CDS", "UTR", "exon"]
    gencode_ht = filter_gencode_ht(
        reference_genome=get_reference_genome(ht.locus).name,
        feature=feature_order + ["gene"],
        genes=gene,
    )

    # The 75bp padding only applies to variants in the specified gene interval
    # (without padding), so we need to filter the gencode HT to only include the gene
    # of interest first.
    ht = filter_by_intervals(
        gencode_ht.filter(gencode_ht.feature == "gene").interval, ht=ht
    )

    for f in feature_order:
        filtered_gencode_ht = gencode_ht.filter(gencode_ht.feature == f)
        filter_count = filtered_gencode_ht.count()
        if filter_count > 0:
            break

    if filter_count == 0:
        raise ValueError(f"No intervals match the gene symbol {gene}")

    ht = filter_by_intervals(
        filtered_gencode_ht.interval, ht=ht, padding_bp=exon_padding_bp
    )

    return ht


def get_age_distribution(
    variant: Optional[str] = None,
    contig: Optional[str] = None,
    position: Optional[int] = None,
    ref: Optional[str] = None,
    alt: Optional[str] = None,
    **kwargs,
) -> hl.Table:
    """
    Get the age distribution of a variant.

    :param variant: Variant string in the format "chr12-235245-A-C" or
        "chr12:235245:A:C". If provided, `contig`, `position`, `ref`, and `alt` are
        ignored.
    :param contig: Chromosome of the variant. Required if `variant` is not provided.
    :param position: Variant position. Required if `variant` is not provided.
    :param ref: Reference allele. Required if `variant` is not provided.
    :param alt: Alternate allele. Required if `variant` is not provided.
    :param kwargs: Additional arguments to pass to `_get_dataset`.
    :return: Table with the age distribution of the variant.
    """
    # Load the Hail Table if not provided
    ht = get_single_variant(variant, contig, position, ref, alt, **kwargs)

    # Age distribution is stored in different structure in different releases
    # Check to see if 'histograms' annotation exists (structure for v4)
    if "histograms" in ht.row:
        ht = ht.select(age_distributions=ht.histograms.age_hists)
    elif "exomes" in ht.row and "genomes" in ht.row:
        # If 'histograms' annotation does not exist, check `data_type`
        ht = ht.select(
            exomes_age_distributions=ht.exomes.histograms.age_hists,
            genomes_age_distributions=ht.genomes.histograms.age_hists,
        )
    else:
        raise ValueError(
            "The age distribution is not available for this dataset. Please check the "
            "selected dataset and try again."
        )

    return ht
