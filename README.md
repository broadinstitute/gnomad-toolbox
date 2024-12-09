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

## Getting started
### Install
pip install -r requirements.txt

### Opening the notebooks
jupyter lab
