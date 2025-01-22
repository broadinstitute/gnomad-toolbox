# gnomad-toolbox

![License](https://img.shields.io/github/license/broadinstitute/gnomad-toolbox)

This repository provides a set of Python functions to simplify working with gnomAD Hail Tables. It includes tools for data access, filtering, and analysis.

> **Disclaimer:** This package is in its early stages of development, and we are actively working on improving it. There may be bugs, and the API may change. Feedback and contributions are welcome.

---

## Repository Structure
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
│   ├── variant.py       # Functions to filter to a specific variant or set of variants.
│   ├── vep.py           # Functions to filter variants based on VEP (Variant Effect Predictor) annotations.
│
├── analysis/
│   ├── general.py       # General analysis functions, such as summarizing variant statistics.
│
├── notebooks/
│   ├── explore_release_data.ipynb      # Guide to loading gnomAD release data.
│   ├── intro_to_filtering_variant_data.ipynb  # Introduction to filtering gnomAD variants.
│   ├── dive_into_secondary_analyses.ipynb  # Analyses with gnomAD data.
```

---

## Set Up Your Environment for Hail and gnomAD Toolbox

This guide provides step-by-step instructions to set up a working environment for using Hail and the gnomAD Toolbox.

> **Note:** We provide this guide to help you set up your environment, but we cannot guarantee that it will work on all systems. If you encounter any issues, you can reach out to us on the [gnomAD Forum](https://discuss.gnomad.broadinstitute.org), and if it is something that we have come across before, we will try to help you out.

### Prerequisites

Ensure you have the following:
- Administrator access to your system to install software.
- Internet connection for downloading dependencies.

> **Note:** Hail 0.2.127+ requires Java 8 or Java 11. If you use an Apple M1/M2 chip, you must have arm64 Java installed. Verify your setup:
>   ```commandline
>   java -version
>   file $JAVA_HOME/bin/java
>   ```
If you don’t have arm64 Java, download it [here](https://www.azul.com/downloads/?os=macos&architecture=arm-64-bit&package=jre#zulu).

---

### Install Miniconda

Miniconda is a lightweight distribution of Conda.

1. Download Miniconda from the [official website](https://docs.anaconda.com/miniconda/install/).
2. Follow the installation instructions described on the download page for your operating system.
3. Verify installation:
   ```commandline
   conda --version
   ```

---

### Create a Conda Environment

1. Create a new environment with a specified Python version:
   ```commandline
   conda create -n gnomad-toolbox python=3.11
   ```
2. Activate the environment:
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
   > **Note:** If you encounter an error like: `Error: pg_config executable not found`, you may need to install the `postgresql` package:
   >   ```commandline
   >   conda install postgresql
   >   ```

### Verify the Setup

Start a Python shell and test if Hail and gnomad_toolbox are working:
```python
import hail as hl
import gnomad_toolbox
hl.init()
print("Hail and gnomad_toolbox setup is complete!")
```

---

## Access gnomAD Data Locally with example notebooks

If you already have experience with gcloud and have no problem running these notebooks, you can skip this section.

### Install the Cloud Storage Connector

Hail uses the Google Cloud Storage Connector to read and write data from Google Cloud Storage. The easiest way to install the connector is to use the `install-gcs-connector` script provided by the Broad Institute:
```commandline
curl -sSL https://broad.io/install-gcs-connector | python3 - --auth-type UNAUTHENTICATED
```

### Use the Example Notebooks

The gnomAD tool-box package includes example notebooks to help you get started with loading and filtering gnomAD data.

1. Copy example notebooks to a new directory:
   ```commandline
   copy-gnomad-toolbox-notebooks /path/to/copy/notebooks/to
   ```
   > **Note:** If the specified directory already exists, you will need to provide a different path, or if you want to overwrite the existing directory, you will need to add the `--overwrite` flag:
   >   ```commandline
   >   copy-gnomad-toolbox-notebooks /path/to/copy/notebooks/to --overwrite
   >   ```

2. Start Jupyter with gnomad-toolbox specific configurations:
   - jupyter Notebook:
     ```commandline
     gnomad-toolbox-jupyter notebook
     ```
   - jupyter Lab:
     ```commandline
     gnomad-toolbox-jupyter lab
     ```

   > **Note:** These commands will start a Jupyter notebook/lab server and open a new tab in your default web browser. The notebook directory containing the example notebooks will be displayed.

3. Open the `explore_release_data.ipynb` notebook to learn how to load gnomAD release data:
   - Run all cells by clicking on the >> button in the toolbar (shown in the image below) or by selecting "Run All" from the "Cell" menu.
      ![jupyter notebook -- run all cells](images/jupyter_run_all.png)

4. Explore the other notebooks:
   - `intro_to_filtering_variant_data.ipynb`: Introduction to filtering variants.
   - `dive_into_secondary_analyses.ipynb`: Examples of some simple analyses using gnomAD data.

5. Try adding your own queries to the notebooks to explore the data further.
> **WARNING:** Avoid running queries on the full dataset as it may take a long time.

---

## Resources

### gnomAD:
- [gnomAD Toolbox Documentation](https://broadinstitute.github.io/gnomad-toolbox/)
- [gnomAD Browser](https://gnomad.broadinstitute.org/)
- [gnomAD Download Page](https://gnomad.broadinstitute.org/downloads)
- [gnomAD Forum](https://discuss.gnomad.broadinstitute.org)

### Hail:
- [Hail Documentation](https://hail.is/docs/0.2/index.html)
- [Hail Discussion Forum](https://discuss.hail.is/)

---

## Contributing

We welcome contributions! Please submit issues and pull requests on our [GitHub repository](https://github.com/broadinstitute/gnomad-toolbox).

---

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.
