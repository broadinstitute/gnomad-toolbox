#!/usr/bin/env python3

"""Generate API reference documentation for a package."""

import importlib
import os
import pathlib
import pkgutil
import re
import sys

import setuptools

REPOSITORY_ROOT_PATH = str(pathlib.Path(os.path.abspath(__file__)).parent.parent)

sys.path.insert(0, REPOSITORY_ROOT_PATH)

ROOT_PACKAGE_PATH = os.path.join(REPOSITORY_ROOT_PATH, "gnomad_toolbox")

DOCS_DIRECTORY = os.path.join(REPOSITORY_ROOT_PATH, "docs")

EXCLUDE_PACKAGES = ["tests"]

PACKAGE_DOC_TEMPLATE = """{title}

{package_doc}

.. toctree::
    :maxdepth: 2

    {module_links}
"""

MODULE_DOC_TEMPLATE = """{title}

{module_doc}

{argparse_doc}

Module Functions
****************

.. automodulesummary:: {module_name}

.. automodule:: {module_name}
    :exclude-members: get_script_argument_parser
"""

ARGPARSE_TEMPLATE = """
.. argparse::
   :ref: {module_name}.get_script_argument_parser
   :prog: {module_name}.py

"""


def module_doc_path(module):
    """Get the path for a module's documentation file."""
    return os.path.join(
        DOCS_DIRECTORY,
        "api_reference",
        re.sub(r"\.py$", ".rst", os.path.relpath(module.__file__, ROOT_PACKAGE_PATH)),
    )

def package_doc_path(package):
    """Get the path for a package's documentation file."""
    return os.path.join(
        DOCS_DIRECTORY,
        "api_reference",
        os.path.relpath(package.__path__[0], ROOT_PACKAGE_PATH),
        "index.rst",
    )


def write_file(path, contents):
    """Write a file after ensuring that the target directory exists."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as out:
        out.write(contents)


def format_title(title):
    """
    Format title for reST.

    reST requires header text to have an underline at least as long as the text.
    """
    underline = "=" * len(title)
    return f"{title}\n{underline}"


def write_module_doc(module_name):
    """Write API reference documentation file for a module."""
    module = importlib.import_module(module_name)

    if hasattr(module, "get_script_argument_parser"):
        argparse_doc = ARGPARSE_TEMPLATE.format(module_name=module_name)
    else:
        argparse_doc = ""

    doc = MODULE_DOC_TEMPLATE.format(
        module_name=module_name,
        title=format_title(module_name),
        module_doc=module.__doc__ or "",
        argparse_doc=argparse_doc,
    )

    doc_path = module_doc_path(module)
    write_file(doc_path, doc)


def write_package_doc(
    package_name,
    package_doc = None,
    doc_path = None,
):
    """Write API reference documentation file for a package."""
    package = importlib.import_module(package_name)

    module_links = []

    for module in pkgutil.iter_modules(package.__path__):
        if module.name in EXCLUDE_PACKAGES:
            continue

        full_module_name = f"{package_name}.{module.name}"
        if module.ispkg:
            write_package_doc(full_module_name)
            module_links.append(f"{module.name} <{module.name}/index>")
        else:
            write_module_doc(full_module_name)
            module_links.append(f"{module.name} <{module.name}>")

    doc = PACKAGE_DOC_TEMPLATE.format(
        title=format_title(package_name),
        package_doc= package_doc or package.__doc__ or "",
        module_links="\n    ".join(module_links),
    )

    doc_path = doc_path or package_doc_path(package)
    write_file(doc_path, doc)


if __name__ == "__main__":
    write_package_doc(
        "gnomad_toolbox",
        package_doc=(
            "This repository provides a set of Python functions to simplify working "
            "with gnomAD Hail Tables. It includes tools for data access, filtering, "
            "and analysis."
        ),
        doc_path=os.path.join(DOCS_DIRECTORY, "api_reference", "index.rst"),
    )
