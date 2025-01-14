"""Functions to import gnomAD data."""

from typing import Callable, Optional, Union

import gnomad.resources.grch37 as grch37_res
import gnomad.resources.grch38 as grch38_res
import hail as hl
from gnomad.resources.resource_utils import VersionedTableResource

DATA_TYPES = ["exomes", "genomes", "joint"]
PEXT_DATA_TYPES = ["base_level", "annotation_level"]
RESOURCES_BY_BUILD = {
    "GRCh37": grch37_res,
    "GRCh38": grch38_res,
}
VARIANT_DATA = {
    "2.1.1": {
        "reference_genome": "GRCh37",
        "data_types": ["exomes", "genomes"],
        "dataset_versions": {
            "vep": "85",
            "gencode": "v19",
            "coverage": {
                "exomes": "2.1",
                "genomes": "2.1",
            },
            "constraint": "2.1.1",
            "pext": "v7",
            "liftover": "2.1.1",
        },
    },
    "4.1": {
        "reference_genome": "GRCh38",
        "data_types": ["exomes", "genomes", "joint"],
        "dataset_versions": {
            "vep": "105",
            "gencode": "v39",
            "coverage": {"exomes": "4.0", "genomes": "3.0.1"},
            "all_sites_an": "4.1",
            "constraint": "4.1",
            "pext": "v7",
            "browser": "4.1",
        },
    },
}
COVERAGE_DATA = {
    "2.1": {"reference_genome": "GRCh37", "data_types": ["exomes", "genomes"]},
    "3.0.1": {"reference_genome": "GRCh38", "data_types": ["genomes"]},
    "4.0": {
        "reference_genome": "GRCh38",
        "data_types": ["exomes"],
    },
}
ALL_SITES_AN_DATA = {
    "4.1": {
        "reference_genome": "GRCh38",
        "data_types": ["exomes", "genomes"],
    }
}
CONSTRAINT_DATA = {
    "2.1.1": {
        "reference_genome": "GRCh37",
        "exome_coverage_field": "median",
        "exome_coverage_cutoff": 30,
        "af_cutoff": 0.001,
        "mane_select": False,
        "canonical": True,
    },
    "4.1": {
        "reference_genome": "GRCh38",
        "exome_coverage_field": "median_approx",
        "exome_coverage_cutoff": 30,
        "af_cutoff": 0.001,
        "mane_select": True,
        "canonical": False,
    },
}
LIFTOVER_DATA = {
    "2.1.1": {
        "reference_genome": "GRCh37",
        "data_types": ["exomes", "genomes"],
    }
}
PEXT_DATA = {
    "v7": {"reference_genome": "GRCh37", "data_types": PEXT_DATA_TYPES},
    "v10": {"reference_genome": "GRCh38", "data_types": PEXT_DATA_TYPES},
}
BROWSER_DATA = {"4.1": {"reference_genome": "GRCh38"}}
SUPPORTED_DATASETS = {
    "variant": {"resource": "public_release", "versions": VARIANT_DATA},
    "coverage": {"resource": "coverage", "versions": COVERAGE_DATA},
    "all_sites_an": {"resource": "all_sites_an", "versions": ALL_SITES_AN_DATA},
    "constraint": {"resource": "constraint", "versions": CONSTRAINT_DATA},
    "liftover": {"resource": "liftover", "versions": LIFTOVER_DATA},
    "pext": {"resource": "pext", "versions": PEXT_DATA},
    "browser": {"resource": "browser", "versions": BROWSER_DATA},
}
SUPPORTED_REFERENCE_DATA = {
    "vep": {
        "resource": "vep_context",
        "versions": {
            "85": {"reference_genome": "GRCh37"},
            "105": {"reference_genome": "GRCh38"},
        },
    },
    "gencode": {
        "resource": "gencode",
        "versions": {
            "v19": {"reference_genome": "GRCh37"},
            "v39": {"reference_genome": "GRCh38"},
        },
    },
}


# TODO: Are there other things we want to store in a session? Like a PASS variants only
#  variable?
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
        self.version = "4.1"
        self.set_default_data()

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

        # Validate version.
        if version not in VARIANT_DATA:
            raise ValueError(f"Version {version} is not a supported gnomAD version.")

        # Validate data type for the version.
        version_info = VARIANT_DATA[version]
        if data_type not in version_info["data_types"]:
            raise ValueError(f"Version {version} for {data_type} is not available.")

        self.data_type = data_type
        self.version = version
        self.reference_genome = version_info["reference_genome"]
        self.compatible_datasets = version_info["dataset_versions"]


# Global gnomad session object
gnomad_session = GnomADSession()


def _get_dataset(
    ht: hl.Table = None,
    dataset: str = "variant",
    data_type: str = None,
    version: str = None,
) -> hl.Table:
    """
    Get gnomAD HT using a Hail Table, specific parameters, or session defaults.

    :param ht: Pre-loaded Hail Table. If provided, other parameters are ignored.
    :param dataset: Dataset type. One of "variant", "all_sites_an", "coverage",
        "gencode". Default is variant.
    :param data_type: Data type (exomes, genomes, or joint). Default is session value.
    :param version: gnomAD version. Default is session value.
    :return: Hail Table for requested dataset, data type, and version.
    """
    # If a pre-loaded Hail Table is provided, return it directly.
    if ht is not None:
        return ht

    # Validate dataset.
    dataset_info = SUPPORTED_DATASETS.get(dataset)
    if dataset_info is None:
        dataset_info = SUPPORTED_REFERENCE_DATA.get(dataset)

    if dataset_info is None:
        raise ValueError(
            f"{dataset} is invalid. Choose from:\n"
            f"\tgnomAD datasets:{list(SUPPORTED_DATASETS.keys())}\n"
            f"\treference datasets: {list(SUPPORTED_REFERENCE_DATA.keys())}"
        )

    # If version or data_type are not provided, use the session information.
    data_type = data_type or gnomad_session.data_type
    if dataset == "variant":
        version = version or gnomad_session.version
    else:
        version = version or gnomad_session.compatible_datasets[dataset]

    # Validate version.
    versions = dataset_info["versions"]
    version_info = versions.get(version)
    if version_info is None:
        version_format = "\t".join(
            f"{v} ({versions[v]['reference_genome']})" for v in versions
        )
        raise ValueError(
            f"Version {version} is not in the supported versions for {dataset}. "
            f"Supported versions: {version_format}"
        )

    # Validate data type.
    data_types = version_info.get("data_types")
    if data_types and data_type and data_type not in data_types:
        raise ValueError(
            f"Version {version} is not available for {data_type} in the {dataset} "
            f"dataset. Available data types: {data_types}."
        )

    # Get the resource for the given build.
    build = version_info["reference_genome"]
    res_build = RESOURCES_BY_BUILD[build]
    res_name = dataset_info["resource"]
    res = getattr(res_build.gnomad, res_name, None) or getattr(
        res_build.reference_data, res_name
    )

    if isinstance(res, Callable):
        res = res(data_type) if data_type else res()

    if isinstance(res, VersionedTableResource):
        res = res.versions[version]

    return res.ht()


def get_gnomad_release(
    dataset: str = "variant",
    data_type: Optional[str] = None,
    version: Optional[str] = None,
) -> hl.Table:
    """
    Get gnomAD HT by dataset, data type, and version.

    Not all combinations of dataset, data type, and version are available and/or
    supported by the toolbox. The table below shows what is supported.

    .. table:: Available versions for each dataset and data type are (as of 2025-1-13)
        :widths: auto

        +--------------+--------------+---------+------------------------------+
        | Dataset      | Genome Build | Version | Data Types                   |
        +==============+==============+=========+==============================+
        | variant      | GRCh37       | 2.1.1   | exomes, genomes              |
        |              +--------------+---------+------------------------------+
        |              | GRCh38       | 4.1     | exomes, genomes, joint       |
        +--------------+--------------+---------+------------------------------+
        | coverage     | GRCh37       | 2.1     | exomes, genomes              |
        |              +--------------+---------+------------------------------+
        |              | GRCh38       | 3.0.1   | genomes                      |
        |              |              +---------+------------------------------+
        |              |              | 4.0     | exomes                       |
        +--------------+--------------+---------+------------------------------+
        | all_sites_an | GRCh38       | 4.1     | exomes, genomes              |
        +--------------+--------------+---------+------------------------------+
        | constraint   | GRCh37       | 2.1.1   | N/A                          |
        |              +--------------+---------+------------------------------+
        |              | GRCh38       | 4.1     | N/A                          |
        +--------------+--------------+---------+------------------------------+
        | liftover     | GRCh37       | 2.1.1   | exomes, genomes              |
        +--------------+--------------+---------+------------------------------+
        | pext         | GRCh37       | v7      | base_level, annotation_level |
        |              +--------------+---------+------------------------------+
        |              | GRCh38       | v10     | base_level, annotation_level |
        +--------------+--------------+---------+------------------------------+
        | browser      | GRCh38       | 4.1     | N/A (joint, but doesn't need |
        |              |              |         | to be specified)             |
        +--------------+--------------+---------+------------------------------+

    :param dataset: Dataset type. One of "variant", "coverage", "all_sites_an",
        "constraint", "liftover", "pext", "browser". Default is "variant".
    :param data_type: Data type. One of "exomes", "genomes", "joint" for all datasets
        except "pext" where it is one of "base_level", "annotation_level". Default is
        the current session data type.
    :param version: gnomAD dataset version. Default is the current session version.
    :return: Hail Table for requested dataset, data type, and version.
    """
    return _get_dataset(dataset=dataset, data_type=data_type, version=version)


def get_compatible_dataset_versions(
    dataset: str, variant_version: Optional[str] = None, data_type: Optional[str] = None
) -> Union[str, dict]:
    """
    Get the compatible version of another datasets for a given gnomAD variant data version.

    :param dataset: Dataset to get the compatible version for.
    :param variant_version: Optional gnomAD variant data version. If not provided, the
        current session version is used.
    :param data_type: Optional data type for the dataset if applicable.
    :return: Compatible version of the dataset for the given variant version.
    """
    # Get the dictionary of compatible versions for the given variant version or
    # the current session version.
    if variant_version is None:
        versions = gnomad_session.compatible_datasets
    else:
        versions = VARIANT_DATA[variant_version]["dataset_versions"]

    # Validate dataset.
    if dataset not in versions:
        raise ValueError(
            f"{dataset} is not available for {variant_version}."
            f"Available datasets: {list(versions.keys())}"
        )

    # If the dataset has multiple data types and a data type is provided, return the
    # version for the data type.
    dataset_version = versions[dataset]
    if data_type and isinstance(dataset_version, dict):
        if data_type not in dataset_version:
            raise ValueError(
                f"{data_type} is not available for {variant_version} {dataset}."
                f"Available data types: {list(dataset_version.keys())}"
            )
        return dataset_version[data_type]

    return dataset_version
