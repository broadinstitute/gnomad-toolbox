"""Functions to import gnomAD data."""

import functools
from typing import Optional

import gnomad.resources.grch37.gnomad as grch37_gnomad
import gnomad.resources.grch38.gnomad as grch38_gnomad
import hail as hl

GNOMAD_BY_BUILD = {
    "GRCh37": grch37_gnomad,
    "GRCh38": grch38_gnomad,
}
DATASETS = {
    "variant": "public_release",
    "all_sites_an": "all_sites_an",
    "coverage": "coverage",
}
DATA_TYPES = ["exomes", "genomes", "joint"]
RELEASES_GLOBAL = {
    "variant": {
        "exomes": "EXOME_RELEASES",
        "genomes": "GENOME_RELEASES",
        "joint": "JOINT_RELEASES",
    },
    "all_sites_an": {
        "exomes": "EXOME_AN_RELEASES",
        "genomes": "GENOME_AN_RELEASES",
    },
    "coverage": {
        "exomes": "EXOME_COVERAGE_RELEASES",
        "genomes": "GENOME_COVERAGE_RELEASES",
    },
}
RELEASES = {
    dataset: {
        data_type: {
            build: (
                None
                if release_global.get(data_type) is None
                else getattr(res, release_global.get(data_type), None)
            )
            for build, res in GNOMAD_BY_BUILD.items()
        }
        for data_type in DATA_TYPES
    }
    for dataset, release_global in RELEASES_GLOBAL.items()
}


class GnomADSession:
    """Class to manage the default data type and version for a gnomAD session."""

    def __init__(self) -> None:
        """
        Initialize a gnomAD session.

        The default data type is exomes and the default version is the current exome
        release.

        :return: None.
        """
        self.data_type = "exomes"
        self.version = grch38_gnomad.CURRENT_EXOME_RELEASE

    def set_default_data(
        self,
        data_type: Optional[str] = None,
        version: Optional[str] = None,
    ) -> None:
        """
        Set default data type and version.

        :param data_type: Data type (exomes, genomes, or joint).
        :param version: gnomAD version.
        :return: None.
        """
        data_type = data_type or self.data_type
        version = version or self.version

        # Validate data type.
        if data_type and data_type not in DATA_TYPES:
            raise ValueError(
                f"Data type {data_type} is invalid. Choose from 'exomes', 'genomes', "
                f"or 'joint'."
            )

        # Get all possible versions.
        possible_versions = functools.reduce(
            lambda x, y: (x or []) + (y or []),
            [
                ds[dt][r]
                for ds in RELEASES.values()
                for dt in ([data_type] if data_type else DATA_TYPES)
                for r in GNOMAD_BY_BUILD
            ],
        )

        # Check version availability.
        if version not in possible_versions:
            raise ValueError(
                f"Version {version} is not available"
                f"{'' if data_type else f' for {data_type}'}. "
            )

        self.data_type = data_type
        self.version = version


# Global gnomad session object
gnomad_session = GnomADSession()


def _get_gnomad_release(
    ht: hl.Table = None,
    dataset: str = "variant",
    data_type: str = None,
    version: str = None,
) -> hl.Table:
    """
    Get gnomAD HT using a Hail Table, specific parameters, or session defaults.

    :param ht: Pre-loaded Hail Table. If provided, other parameters are ignored.
    :param dataset: Dataset type. One of "variant", "all_sites_an", "coverage". Default
        is variant.
    :param data_type: Data type (exomes, genomes, or joint). Default is session value.
    :param version: gnomAD version. Default is session value.
    :return: Hail Table for requested dataset, data type, and version.
    """
    # If a pre-loaded Hail Table is provided, return it directly.
    if ht is not None:
        return ht

    # Use session defaults if parameters are not provided.
    data_type = data_type or gnomad_session.data_type
    version = version or gnomad_session.version

    # Get all releases for the given dataset.
    releases = RELEASES.get(dataset)

    # Validate dataset.
    if releases is None:
        raise ValueError(f"{dataset} is invalid. Choose from {RELEASES.keys()}")

    # Get all releases for the given dataset and data_type.
    data_type_releases = releases.get(data_type)

    # Validate data type.
    if data_type_releases is None:
        raise ValueError(
            f"Data type {data_type} is invalid. Choose from 'exomes', 'genomes', or "
            "'joint'."
        )

    # Check version availability for GRCh38 and GRCh37.
    if data_type_releases["GRCh38"] and version in data_type_releases["GRCh38"]:
        return (
            getattr(grch38_gnomad, DATASETS[dataset])(data_type).versions[version].ht()
        )
    elif data_type_releases["GRCh37"] and version in data_type_releases["GRCh37"]:
        return (
            getattr(grch37_gnomad, DATASETS[dataset])(data_type).versions[version].ht()
        )
    else:
        raise ValueError(
            f"Version {version} is not available for {
                data_type} in the {dataset} dataset. "
            f"Available versions: GRCh38 - {data_type_releases['GRCh38']}, "
            f"GRCh37 - {data_type_releases['GRCh37']}."
        )


def get_gnomad_release(
    dataset: str = "variant",
    data_type: Optional[str] = None,
    version: Optional[str] = None,
) -> hl.Table:
    """
    Get gnomAD HT by dataset, data type,  and version.

    .. table:: Available versions for each dataset and data type are (as of 2024-10-29)
        :widths: auto

        +--------------+-----------------+----------------------------------+----------------------+
        | Dataset      | Data Type       | GRCh38 Versions                  | GRCh37 Versions      |
        +==============+=================+==================================+======================+
        | variant      | exomes          | 4.0, 4.1                         | 2.1, 2.1.1           |
        |              +-----------------+----------------------------------+----------------------+
        |              | genomes         | 3.0, 3.1, 3.1.1, 3.1.2, 4.0, 4.1 | 2.1, 2.1.1           |
        |              +-----------------+----------------------------------+----------------------+
        |              | joint           | 4.1                              | N/A                  |
        +--------------+-----------------+----------------------------------+----------------------+
        | coverage     | exomes          | 4.0                              | 2.1                  |
        |              +-----------------+----------------------------------+----------------------+
        |              | genomes         | 3.0.1                            | 2.1                  |
        +--------------+-----------------+----------------------------------+----------------------+
        | all_sites_an | exomes          | 4.1                              | N/A                  |
        |              +-----------------+----------------------------------+----------------------+
        |              | genomes         | 4.1                              | N/A                  |
        +--------------+-----------------+----------------------------------+----------------------+

    :param data_type: Data type (exomes, genomes, or joint). Default is "exomes".
    :param version: gnomAD version. Default is the current exome release.
    :param dataset: Dataset type. One of "variant", "all_sites_an", "coverage". Default
        is "variant".
    :return: Hail Table for requested dataset, data type, and version.
    """
    return _get_gnomad_release(dataset=dataset, data_type=data_type, version=version)
