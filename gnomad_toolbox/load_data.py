"""Functions to import gnomAD data."""

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
            build: getattr(res, release_global, None)
            for build, res in GNOMAD_BY_BUILD.items()
        }
        for data_type, release_global in data_types.items()
    }
    for dataset, data_types in RELEASES_GLOBAL.items()
}


def get_gnomad_release(
    data_type: str = "exomes",
    version: str = grch38_gnomad.CURRENT_EXOME_RELEASE,
    dataset: str = "variant",
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
