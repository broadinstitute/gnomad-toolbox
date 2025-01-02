"""Functions to filter gnomAD sites HT by VEP annotations."""

from functools import reduce

import hail as hl
from gnomad.resources.grch37.reference_data import gencode as grch37_gencode
from gnomad.resources.grch38.reference_data import gencode as grch38_gencode
from gnomad.utils.vep import CSQ_CODING, LOF_CSQ_SET, filter_vep_transcript_csqs

from gnomad_toolbox.load_data import _get_gnomad_release, gnomad_session

# TODO: I haven't looked over this function yet. Is there anything in gnomad_methods
#  that could be used here? If not, is there anything here that should be moved to
#  gnomad_methods?


def filter_by_csqs(
    csqs: list[str],
    pass_filters: bool = True,
    **kwargs,
) -> hl.Table:
    """
    Filter variants by VEP transcript consequences.

    :param csqs: List of consequences to filter by. It can be specified as the
         categories on the browser: pLoF, Missense / Inframe indel, Synonymous, Other.
    :param pass_filters: Boolean if the variants pass the filters.
    :param kwargs: Arguments to pass to _get_gnomad_release.
    :return: Table with variants with the specified consequences.
    """
    # Load the Hail Table if not provided
    ht = _get_gnomad_release(dataset="variant", **kwargs)

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


def filter_to_plofs(
    gene_symbol: str, select_fields: bool = False, **kwargs
) -> hl.Table:
    """
        Filter to observed pLoF variants that we used to calculate the gene constraint metrics.

        .. note::

                pLOF variants meets the following requirements:
                - High-confidence LOFTEE variants (without any flags),
                - Only variants in the MANE Select transcript,
                - PASS variants that are SNVs with MAF ≤ 0.1%,
                - Exome median depth ≥ 30 (# TODO: This is changing in v4 constraint?)

    **Note: this number should match the number of observed pLOF SNVs on the gene page of gnomAD Browser.**

        :param gene_symbol: Gene symbol.
        :param select_fields: Boolean if the output should be limited to specific fields.
        :return: Table with pLoF variants.
    """
    var_version = kwargs.pop("version", gnomad_session.version)
    var_ht = _get_gnomad_release(dataset="variant", version=var_version, **kwargs)

    # Determine the version of the coverage table
    if var_version.startswith("4."):
        cov_ht = _get_gnomad_release(dataset="coverage", version="4.0", **kwargs)
    elif var_version.startswith("3."):
        cov_ht = _get_gnomad_release(dataset="coverage", version="3.0.1", **kwargs)
    elif var_version.startswith("2."):
        cov_ht = _get_gnomad_release(dataset="coverage", version="2.1", **kwargs)
    else:
        raise ValueError(
            f"Unrecognized version: '{var_version}'. Please specify a valid version."
        )

    # Get the gene interval from gen_ht
    gen_ht = (
        grch37_gencode.ht() if var_version.startswith("2.") else grch38_gencode.ht()
    )
    interval = (
        gen_ht.filter(
            (gen_ht.feature == "gene")
            & (gen_ht.gene_name.lower() == gene_symbol.lower())
        )
        .select()
        .collect()
    )

    if not interval:
        raise ValueError(f"No interval found for gene: {gene_symbol}")

    # Convert to a list of intervals
    interval = [row["interval"] for row in interval]
    var_ht = hl.filter_intervals(var_ht, interval)
    cov_ht = hl.filter_intervals(cov_ht, interval)

    # Filter to high-confidence LOFTEE variants
    var_ht = filter_vep_transcript_csqs(
        var_ht,
        synonymous=False,
        mane_select=True if var_version.startswith("4.") else False,
        # TODO: When this function is applied to DRD2 gene in v4.1, it will get 7 pLoF
        #  variants instead of 8 on the browser and the 4.1 constraint table,
        #  because one of them is not in mane select, nor in canonical transcript.
        genes=[gene_symbol.upper()],
        match_by_gene_symbol=True,
        additional_filtering_criteria=[
            lambda x: (x.lof == "HC")
            & (hl.is_missing(x.lof_flags) | (x.lof_flags == ""))
        ],
    )

    if var_version.startswith("2."):
        allele_type_expr = var_ht.allele_type
        cov_cut_expr = cov_ht[var_ht.locus].median
    else:
        allele_type_expr = var_ht.allele_info.allele_type
        cov_cut_expr = cov_ht[var_ht.locus].median_approx

    var_ht = var_ht.filter(
        (hl.len(var_ht.filters) == 0)
        & (allele_type_expr == "snv")
        & (var_ht.freq[0].AF <= 0.001)
        & (cov_cut_expr >= 30)
    )

    if select_fields:
        var_ht = var_ht.select(
            freq=var_ht.freq[0],
            csq=var_ht.vep.transcript_consequences[0].consequence_terms,
            coverage=cov_ht[var_ht.locus],
        )

    return var_ht
