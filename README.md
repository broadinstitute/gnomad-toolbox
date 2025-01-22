# gnomad-toolbox
This repository provides a set of Python functions to simplify working with gnomAD Hail Tables. It includes tools for
data access, filtering, and analysis.

**Disclaimer: This package is in its early stages of development, and we are actively working on improving it. There may
be bugs, and the API may change. We welcome feedback and contributions.**

## Repository structure
```
gnomad_toolbox/
│
├── load_data.py         # Functions to load gnomAD release Hail Tables.
│
├── filtering/
│   ├── constraint.py    # Functions to filter constraint metrics (e.g., observed/expected ratios).
│   ├── coverage.py      # Functions to filter variants or regions based on coverage thresholds.
│   ├── frequency.py     # Functions to filter variants by allele frequency thresholds.
│   ├── pext.py          # Functions to filter variants using predicted expression (pext) scores.
|   ├── variant.py       # Functions to filter to a specific variant or set of variants.
│   ├── vep.py           # Functions to filter variants based on VEP (Variant Effect Predictor) annotations.
│
├── analysis/
│   ├── general.py       # General analysis functions, such as summarizing variant statistics.
│
├── notebooks/
│   ├── explore_release_data.ipynb    # Jupyter notebook introducing the loading of gnomAD release data.
|   |── intro_to_filtering_variant_data.ipynb    # Jupyter notebook introducing the filtering of gnomAD variant data.
|   |── dive_into_secondary_analyses.ipynb    # Jupyter notebook introducing secondary analyses using gnomAD data.
```

---

## Setting Up Your Environment for Hail and gnomAD Toolbox

This guide provides step-by-step instructions to set up a working environment for
using Hail and the gnomAD Toolbox.

**Disclaimer: We provide this guide to help you set up your environment, but we cannot
guarantee that it will work on all systems. If you encounter any issues, you can
reach out to us on the [gnomAD Forum](https://discuss.gnomad.broadinstitute.org), and
if it is something that we have come across before, we will try to help you out.**

### Prerequisites

Before starting, ensure you have the following:
* Administrator access to your system to install software.
* Internet connection for downloading dependencies.
* **Note: Hail 0.2.127+ requires Java 8 or Java 11 and jupyter labs requires Java
11+, and if you have an Apple M1 or M2, you must have arm64 Java installed, you
can first check your Java version by running:**
  ```commandline
  java -version
  ```
  and check if you have the arm64 Java by running:
  ```commandline
  file $JAVA_HOME/bin/java
  ```
  If you don't have the arm64 Java, you can find it [here](https://www.azul.com/downloads/?os=macos&architecture=arm-64-bit&package=jre#zulu)

### Install Miniconda
Miniconda is a lightweight distribution of Conda that includes just Conda and its dependencies.
1. Download Miniconda for your system from the [official website](https://docs.anaconda.com/miniconda/install/).
2. Follow the installation instructions for your operating system:
   * Linux/macOS: Run the installer script in your terminal:
      ```commandline
      bash Miniconda3-latest-Linux-x86_64.sh
      ```
   * Windows: Run the installer executable and follow the installation wizard.
   **Note: Our team has not tested the gnomad-toolbox on Windows, so we cannot guarantee that it will
   work as expected or support any issues that may arise.**
3. Confirm the installation by running:
   ```commandline
   conda --version
   ```

### Create a Conda Environment
1. Create a new Conda environment with a specific version of Python:
   ```commandline
   conda create -n gnomad-toolbox python=3.11
   ```
2. Activate the Conda environment:
   ```commandline
   conda activate gnomad-toolbox
   ```
3. Install the gnomad-toolbox package and its dependencies:
   * To install the latest version from PyPI:
      ```commandline
      pip install gnomad-toolbox
      ```
   * To install the most up-to-date version from GitHub:
      ```commandline
      pip install git+https://github.com/broadinstitute/gnomad-toolbox@main
      ```
   Note: If you encounter an error like: `Error: pg_config executable not found.`, you may need to install the `postgresql` package:
   ```commandline
   conda install postgresql
   ```

### Verify the Setup
Start a Python shell and test if Hail and gnomad_toolbox are working:
```python
import hail as hl
import gnomad_toolbox
hl.init()
print("Hail and gnomad_toolbox setup is complete!")
```

---

## Accessing gnomAD Data Locally with example notebooks
If you already have experience with gcloud and have no problem running these notebooks,
you can skip this section.

### Install the Cloud Storage Connector
Hail uses the Google Cloud Storage Connector to read and write data from Google Cloud Storage. The easiest way to
install the connector is to use the `install-gcs-connector` script provided by the Broad Institute:
```commandline
curl -sSL https://broad.io/install-gcs-connector | python3 - --auth-type UNAUTHENTICATED
```

### Using the Example Notebooks
The gnomAD tool-box package includes example notebooks to help you get started with
loading and filtering gnomAD data.

1. Run the `copy-gnomad-toolbox-notebooks` command to copy the example notebooks to a new directory:
   ```commandline
   copy-gnomad-toolbox-notebooks /path/to/copy/notebooks
   ```
   Note: If the specified directory already exists, you will need to provide a different path, or if you want to
   overwrite the existing directory, you will need to add the `--overwrite` flag:
   ```commandline
   copy-gnomad-toolbox-notebooks /path/to/copy/notebooks --overwrite
   ```

2. Use the `gnomad-toolbox-jupyter` command to start a Jupyter server:
   * To start jupyter Notebook:
      ```commandline
      gnomad-toolbox-jupyter notebook
      ```
   * To start jupyter Lab:
      ```commandline
      gnomad-toolbox-jupyter lab
      ```

   These commands will start a Jupyter notebook/lab server and open a new tab in your default web browser. The notebook
   directory containing the example notebooks will be displayed.

3. Open the `explore_release_data.ipynb` notebook to learn how to load gnomAD release data. You can run all cells by
clicking on the >> button in the toolbar (shown in the image below) or by selecting "Run All" from the "Cell" menu.
   ![jupyter notebook -- run all cells](images/jupyter_run_all.png)

4. Explore the other notebooks to learn about additional functionalities and analyses you can perform with gnomAD data.

5. Try adding your own queries to the notebooks to explore the data further.
**WARNING: you should avoid running queries on the full dataset as it may take a long time.**

---

## Resources

### gnomAD:
   * [gnomAD Toolbox Documentation](https://broadinstitute.github.io/gnomad-toolbox/)
   * [gnomAD Browser](https://gnomad.broadinstitute.org/)
   * [gnomAD Download Page](https://gnomad.broadinstitute.org/downloads)
   * [gnomAD Forum](https://discuss.gnomad.broadinstitute.org)

### Hail:
   * [Hail Documentation](https://hail.is/docs/0.2/index.html)
   * [Hail Discussion Forum](https://discuss.hail.is/)
