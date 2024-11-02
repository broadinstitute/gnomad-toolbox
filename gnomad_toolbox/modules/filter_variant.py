"""Small functions to filter variants in gnomAD datasets, such as by allele frequency or by variant type."""

from functools import reduce

import hail as hl
from gnomad.utils.vep import CSQ_CODING, LOF_CSQ_SET


def get_variant_count(
    ht: hl.Table,
    afs: list[float] = [0.01, 0.001],
    singletons: bool = False,
    doubletons: bool = False,
) -> dict:
    """
    Count variants with frequency <1%, <0.1%, and singletons (AC == 1).

    .. note:: This function works for gnomAD exomes and genomes datasets, not yet for
              gnomAD joint dataset, since the HT schema is slightly different.

    :param ht: Input Table.
    :param afs: List of allele frequencies cutoffs.
    :param singletons: Include singletons.
    :param doubletons: Include doubletons.
    :return: Dictionary with counts.
    """
    counts = {}

    # Filter to PASS variants.
    ht = ht.filter(hl.len(ht.filters) == 0)
    if singletons:
        n_singletons = ht.aggregate(hl.agg.count_where(ht.freq[0].AC == 1))
        counts["number of singletons"] = n_singletons
    if doubletons:
        n_doubletons = ht.aggregate(hl.agg.count_where(ht.freq[0].AC == 2))
        counts["number of doubletons"] = n_doubletons

    for af in afs:
        n_variants = ht.aggregate(hl.agg.count_where(ht.freq[0].AF < af))
        counts[f"number of variants with AF < {af}"] = n_variants

    # Count variants with frequency <1%, <0.1%, and singletons (AC == 1).
    return counts


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


def filter_by_csqs(
    ht: hl.Table, csqs: list[str], pass_filters: bool = True
) -> hl.Table:
    """
    Filter variants by consequences.

    :param ht: Input Table.
    :param csqs: List of consequences to filter by. It can be specified as the
         categories on the browser: pLoF, Missense / Inframe indel, Synonymous, Other.
    :param pass_filters: Boolean if the variants pass the filters.
    :return: Table with variants with the specified consequences.
    """
    missense_inframe = ["missense_variant", "inframe_insertion", "inframe_deletion"]

    filter_expr = []
    if "lof" in csqs:
        filter_expr.append(
            hl.literal(LOF_CSQ_SET).contains(ht.vep.most_severe_consequence)
        )

    if "synonymous" in csqs:
        filter_expr.append(ht.vep.most_severe_consequence == "synonymous_variant")

    if "missense" in csqs:
        filter_expr.append(
            hl.literal(missense_inframe).contains(ht.vep.most_severe_consequence)
        )

    if "other" in csqs:
        excluded_csqs = hl.literal(
            list(LOF_CSQ_SET) + missense_inframe + ["synonymous_variant"]
        )
        filter_expr.append(~excluded_csqs.contains(ht.vep.most_severe_consequence))

    if "coding" in csqs:
        filter_expr.append(
            hl.literal(CSQ_CODING).contains(ht.vep.most_severe_consequence)
        )

    if len(filter_expr) == 0:
        raise ValueError(
            "No valid consequence specified. Choose from 'lof', 'synonymous', 'missense', 'other'."
        )

    # Combine filter expressions with logical OR
    if len(filter_expr) == 1:
        combined_filter = filter_expr[0]
    else:
        combined_filter = reduce(lambda acc, expr: acc | expr, filter_expr)

    ht = ht.filter(combined_filter)

    if pass_filters:
        ht = ht.filter(hl.len(ht.filters) == 0)

    return ht
