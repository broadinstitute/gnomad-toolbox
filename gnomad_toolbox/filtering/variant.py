"""Functions to filter the gnomAD sites HT to a specific set of variants."""

from typing import Optional, Union

import hail as hl
from gnomad.utils.filtering import filter_to_gencode_cds
from gnomad.utils.reference_genome import get_reference_genome

from gnomad_toolbox.load_data import _get_gnomad_release


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
    :param kwargs: Additional arguments to pass to `_get_gnomad_release`.
    :return: Table with the single variant.
    """
    if not variant and not all([contig, position, ref, alt]):
        raise ValueError(
            "Either `variant` must be provided or all of `contig`, `position`, `ref`, "
            "and `alt`."
        )

    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    # Determine the reference genome build for the ht.
    build = get_reference_genome(ht.locus).name

    # TODO: Move this to gnomad_methods.
    # Parse the variant string if provided.
    try:
        if variant and ":" not in variant:
            contig, position, ref, alt = variant.split("-")
        if all([contig, position, ref, alt]):
            variant = f"{contig}:{position}:{ref}:{alt}"
        variant = hl.eval(hl.parse_variant(variant, reference_genome=build))
    except ValueError:
        raise ValueError(
            f"Invalid variant format: {variant}. Expected format: chr12-235245-A-C "
            f"or chr12:235245:A:C"
        )

    # Filter to the Locus of the variant of interest.
    ht = hl.filter_intervals(
        ht, [hl.interval(variant.locus, variant.locus, includes_end=True)]
    )

    # Filter to the variant of interest.
    ht = ht.filter(ht.alleles == variant.alleles)

    # Check if the variant exists.
    if ht.count() == 0:
        hl.utils.warning(
            f"No variant found at {variant.locus} with alleles {variant.alleles}"
        )

    return ht


def filter_by_intervals(
    intervals: Union[str, list[str]],
    **kwargs,
) -> hl.Table:
    """
    Filter variants by interval(s).

    :param intervals: Interval string or list of interval strings. The interval string
        format has to be "contig:start-end", e.g.,"1:1000-2000" (GRCh37) or
        "chr1:1000-2000" (GRCh38).
    :param kwargs: Arguments to pass to `_get_gnomad_release`.
    :return: Table with variants in the interval(s).
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    # Determine the reference genome build for the ht.
    build = get_reference_genome(ht.locus).name

    if isinstance(intervals, str):
        intervals = [intervals]

    if build == "GRCh38" and any([not i.startswith("chr") for i in intervals]):
        raise ValueError("Interval must start with 'chr' for GRCh38 reference genome.")

    ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval(i, reference_genome=build) for i in intervals]
    )

    return ht


def filter_by_gene_symbol(gene: str, exon_padding_bp: int = 75, **kwargs) -> hl.Table:
    """
    Filter variants in a gene.

    .. note::

           This function is to match the number of variants that you will get in the
           gnomAD browser, which only focus on variants in "CDS" regions plus
           75bp (default of `exon_padding_bp`) up- and downstream.

           However, gnomAD browser used a preprocessed Gencode file which excluded
           46 genes on chrY that share the same gene id as chrX. For example,
           if you use this function to filter "ASMT" gene, you will get more variants
           than shown in the gnomAD browser.

    :param gene: Gene symbol.
    :param exon_padding_bp: Number of base pairs to pad the CDS intervals. Default is
        75bp.
    :param kwargs: Arguments to pass to `_get_gnomad_release`.
    :return: Table with variants in the gene.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)
    ht = filter_to_gencode_cds(ht, genes=gene, padding=exon_padding_bp)

    return ht
