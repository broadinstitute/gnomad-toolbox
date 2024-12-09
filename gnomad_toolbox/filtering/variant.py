"""Functions to filter the gnomAD sites HT to a specific set of variants."""

from typing import Optional

import hail as hl
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


def filter_by_interval(ht: hl.Table, interval: str) -> hl.Table:
    """
    Filter variants by interval.

    :param ht: Input Table.
    :param interval: Interval string. Format has to be "chr:start-end", e.g., "1:1000-2000".
    :return: Table with variants in the interval.
    """
    if ht.locus.dtype.reference_genome.name == "GRCh38":
        interval = "chr" + interval
    ht = hl.filter_intervals(
        ht,
        [
            hl.parse_locus_interval(
                interval,
                reference_genome=(
                    "GRCh38"
                    if ht.locus.dtype.reference_genome.name == "GRCh38"
                    else "GRCh37"
                ),
            )
        ],
    )
    return ht


def filter_by_gene_symbol(ht: hl.Table, gene: str) -> hl.Table:
    """
    Filter variants in a gene.

    .. note::
           This function is to match the number of variants that you will get in the
           gnomAD browser, which only focus on variants in "CDS" regions plus 75bp
           up- and downstream. This is not the same as filtering by gene symbol with
           our `filter_vep_transcript_csqs` function, which will include all variants.

    :param ht: Input Table.
    :param gene: Gene symbol.
    :return: Table with variants in the gene.
    """
    # Make gene symbol uppercase
    gene = gene.upper()

    if ht.locus.dtype.reference_genome.name == "GRCh37":
        gene_ht = hl.read_table(
            "gs://gcp-public-data--gnomad/resources/grch37/browser/gnomad"
            ".genes.GRCh37.GENCODEv19.ht"
        )
    else:
        gene_ht = hl.read_table(
            "gs://gcp-public-data--gnomad/resources/grch38/browser/gnomad"
            ".genes.GRCh38.GENCODEv39.ht"
        )

    gene_ht = gene_ht.annotate(
        cds_intervals=hl.array(
            gene_ht.exons.filter(lambda exon: exon.feature_type == "CDS")
        ).map(
            lambda exon: hl.locus_interval(
                hl.if_else(
                    gene_ht.interval.start.dtype.reference_genome.name == "GRCh38",
                    "chr" + gene_ht.chrom,
                    gene_ht.chrom,
                ),
                exon.start - 75,
                exon.stop + 75,
                reference_genome=gene_ht.interval.start.dtype.reference_genome,
                includes_end=True,
            )
        )
    )

    intervals = gene_ht.filter(gene_ht.gencode_symbol == gene).cds_intervals.collect()[
        0
    ]

    ht = hl.filter_intervals(ht, intervals)

    return ht
