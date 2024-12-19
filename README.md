# gnomad-toolbox
This repository provides a set of Python functions to simplify working with gnomAD Hail Tables. It includes tools for data access, filtering, and analysis.

## Repository structure
```
ggnomad_toolbox/
│
├── load_data.py         # Functions to load gnomAD release Hail Tables.
│
├── filtering/
│   ├── __init__.py
│   ├── constraint.py    # Functions to filter constraint metrics (e.g., observed/expected ratios).
│   ├── coverage.py      # Functions to filter variants or regions based on coverage thresholds.
│   ├── frequency.py     # Functions to filter variants by allele frequency thresholds.
│   ├── pext.py          # Functions to filter variants using predicted expression (pext) scores.
|   ├── variant.py       # Functions to filter to a specific variant or set of variants.
│   ├── vep.py           # Functions to filter variants based on VEP (Variant Effect Predictor) annotations.
│
├── analysis/
│   ├── __init__.py
│   ├── general.py       # General analysis functions, such as summarizing variant statistics.
│
├── notebooks/
│   ├── intro_to_release_data.ipynb    # Jupyter notebook introducing the loading of gnomAD release data.
```

## Setting Up Your Environment for Hail and gnomAD Toolbox

This guide provides step-by-step instructions to set up a working environment for
using Hail and the gnomad_toolbox.

Prerequisites

Before starting, ensure you have the following:
* Administrator access to your system to install software.
* Internet connection for downloading dependencies.


## Part 1: Setting Up Your Environment

### Install Miniconda
Miniconda is a lightweight distribution of Conda that includes just Conda and its dependencies.
1. Download Miniconda for your system from the [official website](https://docs.anaconda.com/miniconda/install/).
2. Follow the installation instructions for your operating system:

	•	Linux/macOS: Run the installer script in your terminal:
    ```
    bash Miniconda3-latest-Linux-x86_64.sh
    ```
    •	Windows: Run the installer executable and follow the installation wizard.
        (# TODO: Will we encourage users to use Windows?)
3. Confirm the installation by running:
   ```
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
3. Clone the gnomad-toolbox repository and install the dependencies:
   ```commandline
   cd /path/to/your/directory
   git clone https://github.com/broadinstitute/gnomad-toolbox.git
   cd gnomad-toolbox
   pip install -r requirements.txt
    ```
   You might encounter errors when installing the dependencies, such as `pg_config
   executable not found`. If so, you may need to install additional system packages.
   For example, on Ubuntu, you can install the `libpq-dev` package:
   ```commandline
    sudo apt-get install libpq-dev
    ```
   or on macOS, you can install the `postgresql` package:
   ```commandline
    brew install postgresql
    ```

### Verify the Setup
1.	Start a Python shell and test if Hail and gnomad_toolbox are working:
```commandline
import hail as hl
from gnomad_toolbox.analysis.general import get_variant_count_by_freq_bin
hl.init()
print("Hail and gnomad_toolbox setup is complete!")
```
Or open the notebooks:
```commandline
jupyter lab
```

## Part2: Accessing gnomAD Data Locally with example notebooks
If you already have experience with gcloud and have no problem running these notebooks,
you can skip this section.

### Install Google Cloud SDK (gcloud)

The Google Cloud SDK is required to interact with Google Cloud services and access gnomAD public data locally.
1. Follow the official Google Cloud SDK installation [guide](https://cloud.google.
   com/sdk/docs/install) for your operating system.
2. After installation, initialize gcloud to log in and set up your default project:
   ```
   gcloud init
   ```
3. You can check your gcloud config by:
   ```
   gcloud config list
   ```
   or set the default project:
   ```
   gcloud config set project {YOUR_PROJECT_ID}
   ```

## Configure a Service Account
You will need to create a service account in gcloud console IAM & Admin or using
gcloud CLI. Then you can create a key for service account and set the GOOGLE_APPLICATION_CREDENTIALS
variable to the path of the key file.
   ```commandline
    gcloud iam service-accounts keys create hail-local-sa-key.json --iam-account {YOUR_SERVICE_ACCOUNT}

    export GOOGLE_APPLICATION_CREDENTIALS=./hail-local-sa-key.json
   ```
Now, you can access gnomAD data locally using the gnomad_toolbox functions, however,
avoid running queries on the full dataset as it may take a long time and consume a
lot of resources, and most importantly, it may incur costs.
