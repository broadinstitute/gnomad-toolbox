"""
Functions to import non-gnomAD reference data.

This module provides functions to load and import various reference datasets.

Currently, this module includes functions to load ClinVar data from gnomad_methods or by downloading from NCBI.
"""

import hail as hl
import requests
from gnomad.resources.resource_utils import (
    NO_CHR_TO_CHR_CONTIG_RECODING,
    import_sites_vcf,
)

CLINVAR_FTP_URL = {
    "GRCh37": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar.vcf.gz",
    "GRCh38": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz",
}
"""
Path to the latest weekly ClinVar VCF release for GRCh37 and GRCh38.
"""


def _load_clinvar_resource(build: str) -> hl.Table:
    """
    Load ClinVar resource from gnomad_methods based on reference genome build.

    :param build: Reference genome build ('GRCh37' or 'GRCh38').
    :return: Hail Table with ClinVar data containing fields like `CLNSIG`, `CLNREVSTAT`,
        `CLNDN`, `CLNDISDB`, etc.
    """
    if build == "GRCh37":
        import gnomad.resources.grch37 as grch37_res

        return grch37_res.reference_data.clinvar.ht()
    elif build == "GRCh38":
        import gnomad.resources.grch38 as grch38_res

        return grch38_res.reference_data.clinvar.ht()
    else:
        raise ValueError(
            f"Unsupported reference genome build: {build}. Must be 'GRCh37' or 'GRCh38'."
        )


def _import_clinvar_vcf(
    build: str,
    clinvar_download_path: str,
    output_path: str,
    overwrite: bool = False,
) -> hl.Table:
    """
    Import latest weekly ClinVar VCF release to Hail Table.

    :param build: Reference genome build ('GRCh37' or 'GRCh38').
    :param clinvar_download_path: Path to download the ClinVar VCF.
    :param output_path: Path to the output Hail Table.
    :param overwrite: If True, overwrite existing ClinVar Hail Table.
    :return: Hail Table with ClinVar data.
    """
    response = requests.get(CLINVAR_FTP_URL[build], stream=True)
    response.raise_for_status()

    with open(clinvar_download_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    import_args = {
        "path": clinvar_download_path,
        "force_bgz": True,
        "min_partitions": 100,
        "reference_genome": build,
        "skip_invalid_loci": True,
    }

    # ClinVar GRCh38 VCFs do not have "chr" prefix in contig names
    if build == "GRCh38":
        import_args["contig_recoding"] = NO_CHR_TO_CHR_CONTIG_RECODING

    clinvar_ht = import_sites_vcf(**import_args)

    # Remove rows with only a single allele, write table out, and return
    clinvar_ht = clinvar_ht.filter(hl.len(clinvar_ht.alleles) > 1)
    return clinvar_ht.checkpoint(output_path, overwrite=overwrite)
