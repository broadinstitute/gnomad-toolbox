"""Functions to filter the gnomAD sites HT to a specific set of variants."""

from typing import Optional, Union

import hail as hl
from gnomad.utils.parse import parse_variant
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
    Filter variants in a gene.

    .. note::

           This function is to match the number of variants that you will get in the
           gnomAD browser, which only focus on variants in "CDS" regions plus
           75bp (default of `exon_padding_bp`) up- and downstream.

    :param gene: Gene symbol.
    :param exon_padding_bp: Number of base pairs to pad the CDS intervals. Default is
        75bp.
    :param kwargs: Arguments to pass to `_get_gnomad_release`.
    :return: Table with variants in the gene.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

    # Determine the reference genome build for the ht.
    build = get_reference_genome(ht.locus).name

    # Make gene symbol uppercase
    gene = gene.upper()

    # TODO: Create a resource for this in gnomad_methods (is it different from our
    #  current gencode resources?
    # gene_ht = hl.read_table(
    #    f"gs://gcp-public-data--gnomad/resources/{build.lower()}/browser/gnomad"
    #    f".genes.{build}.GENCODEv{'19' if build == 'GRCh37' else '39'}.ht"
    # )

    # TODO: This actually takes a while to run locally for a single gene. Is there a
    #  way to speed this up?

    # Filter to the gene of interest.
    # gene_ht = gene_ht.filter(gene_ht.gencode_symbol == gene)

    # Get the CDS intervals for the gene.
    # chrom_expr = hl.if_else(build == "GRCh38", "chr" + gene_ht.chrom, gene_ht.chrom)
    # intervals = gene_ht.aggregate(
    #    hl.agg.explode(
    #        lambda exon: hl.agg.collect(
    #            hl.locus_interval(
    #                chrom_expr,
    #                exon.start - 75,
    #                exon.stop + 75,
    #                reference_genome=build,
    #                includes_end=True,
    #            )
    #        ),
    #        gene_ht.exons.filter(lambda exon: exon.feature_type == "CDS"),
    #    )
    # )

    # TODO: Consider this alternative approach to get the intervals from gencode. That
    #  is not too bad time wise

    from gnomad.resources.grch38.reference_data import gencode

    gencode_ht = gencode.ht()
    gencode_ht = gencode_ht.filter(
        (gencode_ht.gene_name == gene) & ((gencode_ht.feature == "CDS"))
    )
    intervals = hl.locus_interval(
        gencode_ht.interval.start.contig,
        gencode_ht.interval.start.position - exon_padding_bp,
        gencode_ht.interval.end.position + exon_padding_bp,
        includes_start=gencode_ht.interval.includes_start,
        includes_end=gencode_ht.interval.includes_end,
        reference_genome=build,
    ).collect()

    ht = hl.filter_intervals(ht, intervals)

    return ht
