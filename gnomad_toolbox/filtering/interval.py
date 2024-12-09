"""Functions to filter the gnmoAD sites HT by interval."""

import hail as hl


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
