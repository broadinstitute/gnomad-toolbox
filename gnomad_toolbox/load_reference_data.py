"""Functions to import non-gnomAD reference data."""

import hail as hl
import wget
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
    import_args = {
        "path": wget.download(CLINVAR_FTP_URL[build], out=clinvar_download_path),
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
