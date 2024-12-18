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

# TODO: Add fully detailed info about how to install and open the notebooks.
## Getting started
### Install
pip install -r requirements.txt

### Opening the notebooks
jupyter lab

# Setting Up Your Environment for Hail and gnomAD Toolbox

This guide provides step-by-step instructions to set up a working environment for using Hail and the gnomad_toolbox. The steps include installing Miniconda, creating a Conda environment, installing the necessary dependencies, and configuring the service account.

Prerequisites

Before starting, ensure you have the following:
	•	Administrator access to your system to install software.
	•	Internet connection for downloading dependencies.

## Step 1: Install Miniconda
Miniconda is a lightweight distribution of Conda that includes just Conda and its dependencies.
1. Download Miniconda for your system from the [official website](https://docs.anaconda.com/miniconda/install/).
2. Follow the installation instructions for your operating system:
	•	Linux/macOS: Run the installer script in your terminal:
```
bash Miniconda3-latest-Linux-x86_64.sh
```
    •	Windows: Run the installer executable and follow the installation wizard.
3. Confirm the insstallation by running:
```
conda --version
```

## Step 2: Install Google Cloud SDK (gcloud)

The Google Cloud SDK is required to interact with Google Cloud services, including setting up authentication for Hail.
1. Follow the official Google Cloud SDK installation guide for your operating system.
2. After installation, initialize gcloud:
```
gcloud init
```
3. Authenticate with your Google account:
```
gcloud auth login
```
4. Set the default project:
```
gcloud config set project broad-mpg-gnomad
```

## Step 43: Configure the Service Account
```commandline
gcloud iam service-accounts keys create hail-local-sa-key.json --iam-account hail-local-sa@broad-mpg-gnomad.iam.gserviceaccount.com
export GOOGLE_APPLICATION_CREDENTIALS=path-to-your-key/hail-local-sa-key.json
```

## Step 4: Create a Conda Environment
1. Create a new Conda environment with a specific version of Hail (here we use Hail
   0.2.132 and Python 3.11):
```commandline
conda create -n gnomad-toolbox hail=0.2.132,python=3.11
```
2. Activate the Conda environment:
```commandline
conda activate gnomad-toolbox
```
3. Clone the gnomad-toolbox repository:
```commandline
cd /path/to/your/directory
git clone https://github.com/broadinstitute/gnomad-toolbox.git
cd gnomad-toolbox
pip install -r requirements.txt
pip install git+https://github.com/broadinstitute/gnomad-toolbox.git
```

## Step 5: Verify the Setup
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
