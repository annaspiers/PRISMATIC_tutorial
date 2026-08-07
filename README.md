<img src="docs/img/logo.png" align="right" width="25%"/>

[![Python Package using Conda](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/RS-PRISMATIC/preprocessing/actions/workflows/python-package-conda.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/RS-PRISMATIC/preprocessing/HEAD)

# PRISMATIC workflow tutorial

<a href="https://de.cyverse.org/apps/de/29496b2c-86f4-11f1-a465-008cfa5ae3e1/launch?saved-launch-id=e5434bd9-15bc-426a-92ef-009b1d4851a7" target="_blank" rel="noopener noreferrer"><img src="https://de.cyverse.org/Powered-By-CyVerse-blue.svg"></a>

This is a tutorial version of the [PRISMATIC](https://github.com/annaspiers/PRISMATIC) workflow, demonstrating the end-to-end pipeline on a single 1 km × 1 km remote sensing tile. The full PRISMATIC pipeline generates initial conditions for [FATES](https://fates-users-guide.readthedocs.io/) (an open source, dynamic vegetation model) using NEON (National Ecological Observatory Network) remote sensing and field data, but this tutorial is scoped to one tile at one site rather than multi-site, multi-year runs.

The pipeline combines LiDAR point cloud processing, hyperspectral imagery, and forest inventory data to estimate forest structure (e.g. canopy layers) and classify plant functional types (PFTs) to produce FATES cohort/patch files.

Workshop presented at the annual Ecological Society of America Meeting in July 2026 in Salt Lake City, UT by A Spiers.

## Container access
1. Sign up for CyVerse account: https://user.cyverse.org/signup
1. Enroll in workshop: https://user.cyverse.org/workshops/214
1. Wait for approval (subscriptions approved automatically each hour)
1. Go to Discovery Environment: https://de.cyverse.org/dashboard
1. Search "PRISMATIC Tutorial"
1. Launch app

## Running the Tutorial

The tutorial is presented as a Jupyter notebook (`prismatic_workshop.ipynb`) that walks through each step sequentially for a single tile. Individual pipeline steps can also be run via the Python scripts in `initialize/`.

The environment (Python 3.10, conda-forge) includes PDAL, GDAL, geopandas, rasterio, scikit-learn, rpy2, hydra-core, and R packages (randomForest, caret, terra, sf). R helper scripts are called from Python via rpy2 using helper scripts (*_helper.R) in initialize/.


### Processing Modules (initialize/)

Each module contains one or more pipeline steps:
- `inventory.py` key functions: `download_trait_table()`, `download_veg_structure_data()`, `prep_veg_structure()`
- `plots.py` key functions: `download_polygons()`, `prep_polygons()` 
- `lidar.py` key functions: `download_lidar()`, `download_aop_bbox()`, `normalize_laz()`, `clip_lidar_by_plots()`
- `biomass.py` key functions: `prep_biomass()` 
- `lad.py` key functions: `prep_lad()` generates leaf area density profiles by size class 
- `hyperspectral.py` key functions: `download_hyperspectral()`, `prep_aop_imagery()`, `extract_spectra_from_polygon()`, `train_pft_classifier()`
- `generate_initial_conditions.py` key functions: `generate_initial_conditions()`

### Utilities (utils/)

- `download_functions.py` — NEON AOP file download helpers
- `allometry.py` — biomass allometry equations
- `neon_aop_hyperspectral.py` — hyperspectral data loading/manipulation
- `apply_brdf_corrections.py` — topographic/BRDF corrections for hyperspectral imagery
- `plot_partition.py` — spatial plot partitioning pipeline)
