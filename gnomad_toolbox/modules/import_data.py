"""Functions to import gnomAD data."""

import hail as hl
from gnomad.resources.grch37.gnomad import EXOME_RELEASES as GRCh37_EXOME_RELEASES
from gnomad.resources.grch37.gnomad import GENOME_RELEASES as GRCh37_GENOME_RELEASES
from gnomad.resources.grch37.gnomad import public_release as grch37_public_release
from gnomad.resources.grch38.gnomad import EXOME_RELEASES as GRCh38_EXOME_RELEASES
from gnomad.resources.grch38.gnomad import GENOME_RELEASES as GRCh38_GENOME_RELEASES
from gnomad.resources.grch38.gnomad import JOINT_RELEASES as GRCh38_JOINT_RELEASES
from gnomad.resources.grch38.gnomad import public_release as grch38_public_release


def get_ht_by_datatype_and_version(
    data_type: str = "exomes", version: str = "4.1"
) -> hl.Table:
    """
    Get gnomAD HT by data type and version.

    .. note::

       Available versions for each data type are (as of 2024-10-29):

       ::

           | Data Type       | GRCh38 Versions                  | GRCh37 Versions      |
           |-----------------|----------------------------------|----------------------|
           | exomes          | 4.0, 4.1                         | 2.1, 2.1.1           |
           | genomes         | 3.0, 3.1, 3.1.1, 3.1.2, 4.0, 4.1 | 2.1, 2.1.1           |
           | joint           | 4.1                              | N/A                  |

    :param data_type: Data type (exomes, genomes, or joint).
    :param version: gnomAD version.
    :return: Hail Table.
    """
    # Mapping data types to version sets for GRCh38 and GRCh37
    versions_by_type = {
        "exomes": (GRCh38_EXOME_RELEASES, GRCh37_EXOME_RELEASES),
        "genomes": (GRCh38_GENOME_RELEASES, GRCh37_GENOME_RELEASES),
        "joint": (GRCh38_JOINT_RELEASES, []),
    }

    # Validate data type
    if data_type not in versions_by_type:
        raise ValueError(
            f"Data type {data_type} is invalid. Choose from 'exomes', 'genomes', or 'joint'."
        )

    # Get GRCh38 and GRCh37 versions for the given data type
    grch38_versions, grch37_versions = versions_by_type[data_type]

    # Check version availability for GRCh38 and GRCh37
    if version in grch38_versions:
        return grch38_public_release(data_type).ht()
    elif version in grch37_versions:
        return grch37_public_release(data_type).ht()
    else:
        raise ValueError(
            f"Version {version} is not available for {data_type}. "
            f"Available versions: GRCh38 - {grch38_versions}, GRCh37 - {grch37_versions}."
        )
