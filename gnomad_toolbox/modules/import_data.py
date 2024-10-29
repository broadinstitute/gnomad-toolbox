"""Functions to import gnomAD data."""

import hail as hl
from gnomad.resources.grch37.gnomad import public_release as grch37_public_release
from gnomad.resources.grch38.gnomad import public_release as grch38_public_release


def get_ht_by_datatype_and_version(
    data_type: str = "exomes", version: str = "4.1"
) -> hl.Table:
    """
    Get gnomAD HT by data type and version.

    .. note: Available versions for each data type are:

    | Data Type       | GRCh38 Versions                  | GRCh37 Versions      |
    |-----------------|----------------------------------|----------------------|
    | exomes          | 4.0, 4.1                         | 2.1, 2.1.1           |
    | genomes         | 3.0, 3.1, 3.1.1, 3.1.2, 4.0, 4.1 | 2.1, 2.1.1           |
    | joint           | 4.1                              | N/A                  |


    :param data_type: Data type (exomes or genomes or joint).
    :param version: gnomAD version.
    :return: Hail Table.
    """
    if version in ["2.1", "2.1.1"]:
        return grch37_public_release(data_type).ht()
    elif version in ["3.0", "3.1", "3.1.1", "3.1.2", "4.0", "4.1"]:
        return grch38_public_release(data_type).ht()
    else:
        raise ValueError(f"Version {version} not found for data type {data_type}.")
