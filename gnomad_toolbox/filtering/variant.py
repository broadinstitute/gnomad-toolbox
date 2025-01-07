"""Functions to filter the gnomAD sites HT to a specific set of variants."""

from typing import Optional, Union

import hail as hl
from gnomad.resources.grch37.gnomad import browser_gene as browser_gene_grch37
from gnomad.resources.grch38.gnomad import browser_gene as browser_gene_grch38
from gnomad.utils.filtering import filter_to_gencode_cds
from gnomad.utils.parse import parse_variant
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vep import filter_vep_transcript_csqs

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
    Filter variants in a gene by gene symbol.

    .. note::
       This is to match the browser display, which includes variants in the CDS +
       75bp padding for protein-coding genes, and UTR/exons + 75bp padding for
       non-protein-coding genes.

    :param gene: Gene symbol.
    :param exon_padding_bp: Number of base pairs to pad the CDS intervals. Default is 75bp.
    :param kwargs: Arguments to pass to `_get_gnomad_release`.
    :return: Table with variants in the gene.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    # Load gene information
    gene_ht = (
        browser_gene_grch37().ht()
        if ht.locus.dtype.reference_genome.name == "GRCh37"
        else browser_gene_grch38().ht()
    )

    # Pre-filter to the specified gene (case-insensitive)
    gene = gene.upper()
    gene_ht = gene_ht.filter(gene_ht.gencode_symbol == gene)

    # First get the variants in the gene region
    ht = hl.filter_intervals(ht, hl.array(gene_ht.interval.take(1)))

    # Get intervals based on feature type (CDS > UTR > Exons) with padding
    def get_intervals(feature_type: str) -> hl.expr.ArrayExpression:
        return hl.array(
            gene_ht.exons.filter(lambda exon: exon.feature_type == feature_type)
        ).map(
            lambda exon: hl.locus_interval(
                hl.if_else(
                    gene_ht.interval.start.dtype.reference_genome.name == "GRCh38",
                    "chr" + gene_ht.chrom,
                    gene_ht.chrom,
                ),
                exon.start - exon_padding_bp,
                exon.stop + exon_padding_bp,
                reference_genome=gene_ht.interval.start.dtype.reference_genome,
                includes_start=True,
                includes_end=True,
            )
        )

    cds_intervals = get_intervals("CDS")
    utr_intervals = get_intervals("UTR")
    exon_intervals = get_intervals("exon")

    # Determine which intervals to use
    gene_ht = gene_ht.annotate(
        intervals=hl.if_else(
            hl.len(cds_intervals) > 0,
            cds_intervals,
            hl.if_else(hl.len(utr_intervals) > 0, utr_intervals, exon_intervals),
        )
    )

    intervals = gene_ht.intervals.take(1)[0]

    ht = hl.filter_intervals(ht, intervals)

    # Additional filtering (e.g., with VEP consequences)
    ht = filter_vep_transcript_csqs(
        ht, genes=[gene], synonymous=False, match_by_gene_symbol=True
    )

    return ht
