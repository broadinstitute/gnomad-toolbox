"""Setup script."""

import setuptools

with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

install_requires = []
with open("requirements.txt", "r") as requirements_file:
    for req in (line.strip() for line in requirements_file):
        if req != "hail":
            install_requires.append(req)


setuptools.setup(
    name="gnomad_toolbox",
    version="0.0.1",
    author="The Genome Aggregation Database",
    author_email="gnomad@broadinstitute.org",
    description="Toolbox to help users process gnomAD release files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/broadinstitute/gnomad-toolbox",
    packages=setuptools.find_namespace_packages(include=["gnomad_toolbox.*"]),
    project_urls={
        "Documentation": "https://broadinstitute.github.io/gnomad-toolbox/",
        "Source Code": "https://github.com/broadinstitute/gnomad-toolbox",
        "Issues": "https://github.com/broadinstitute/gnomad-toolbox/issues",
        "Change Log": "https://github.com/broadinstitute/gnomad-toolbox/releases",
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD 3-Clause License",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.9",
    install_requires=install_requires,
)
