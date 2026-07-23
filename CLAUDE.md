# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Does

PRISMATIC generates initial conditions for [FATES](https://fates-users-guide.readthedocs.io/) (a dynamic vegetation model) from NEON (National Ecological Observatory Network) remote sensing and field data. It combines LiDAR point cloud processing, hyperspectral imagery, and forest inventory data to classify plant functional types (PFTs) via a Random Forest ML model and produce FATES cohort/patch files.

Supported initialization types (set via `ic_type` in `conf/others/others.yaml`):
- `field_inv_plots` — from NEON forest inventory plots
- `rs_inv_plots` — remote sensing over inventory plots
- `rs_tower_ftpt` — remote sensing over eddy covariance tower footprint
- `rs_random_plots` — remote sensing over random plots
- `rs_wall2wall` — remote sensing over gridded wall-to-wall plots (default)

## Running the Pipeline

```bash
# Run all configured sites
python main.py

# Run a specific site
python main.py sites.global_run_params.run=SJER

# Run multiple specific sites
python main.py "sites.global_run_params.run=[SJER, SOAP]"

# Force rerun of a specific step for all sites
python main.py sites.global_run_params.force_rerun.prep_biomass=True

# Force rerun of a step for a specific site/year
python main.py sites.SJER.2019.force_rerun.prep_biomass=True

# Submit to Perlmutter HPC via SLURM
sbatch submit_main.sh
```

The pipeline uses Hydra for configuration, so any config key can be overridden on the command line.

## Linting

```bash
# Error-level checks (syntax errors, undefined names)
flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics

# Style checks (warnings only, non-blocking)
flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
```

There are no active tests (pytest is present but disabled in CI).

## Environment Setup

```bash
conda env create -f environment.yml
```

The environment (Python 3.10, conda-forge) includes PDAL, GDAL, geopandas, rasterio, scikit-learn, rpy2, hydra-core, and R packages (randomForest, caret, terra, sf). R helper scripts are called from Python via rpy2.

Data paths are configured in `conf/paths/paths.yaml` pointing to `/pscratch/sd/a/aspiers/data/` on Perlmutter.

## Architecture

### Configuration (Hydra)

All configuration flows through `conf/config.yaml`, which composes:
- `conf/paths/paths.yaml` — raw/intermediate/final data paths
- `conf/sites/sites.yaml` — per-site, per-year parameters and `force_rerun` flags
- `conf/others/others.yaml` — global processing parameters (`ic_type`, `hs_type`, `n_plots`, `plot_length`, `ntree`, `coords_bbox`, etc.)

### Processing Modules (`initialize/`)

Each module contains one or more pipeline steps. `main.py` calls them in order:

| Module | Key functions |
|---|---|
| `inventory.py` | `download_trait_table`, `download_veg_structure_data`, `prep_veg_structure` |
| `plots.py` | `download_polygons`, `prep_polygons` |
| `lidar.py` | `download_lidar`, `download_aop_bbox`, `normalize_laz`, `clip_lidar_by_plots` |
| `biomass.py` | `prep_biomass` |
| `lad.py` | `prep_lad` — leaf area density profiles by size class |
| `hyperspectral.py` | `download_hyperspectral`, `prep_aop_imagery`, `extract_spectra_from_polygon`, `train_pft_classifier` |
| `generate_initial_conditions.py` | `generate_initial_conditions` |

Many steps call R via rpy2 using helper scripts (`*_helper.R`) co-located in `initialize/`.

### Caching and Rerun Logic (`utils/utils.py`)

`build_cache_site()` checks whether intermediate outputs already exist and skips recomputation. Steps are forced to rerun either via `force_rerun` flags in `conf/sites/sites.yaml` or via CLI overrides. The `force_rerun()` decorator wraps individual functions.

### Utilities (`utils/`)

- `download_functions.py` — NEON AOP file download helpers
- `allometry.py` — biomass allometry equations
- `neon_aop_hyperspectral.py` — hyperspectral data loading/manipulation
- `apply_brdf_corrections.py` — topographic/BRDF corrections for hyperspectral imagery
- `plot_partition.py` — spatial plot partitioning
- `logger.py` — logging setup

### Current State of `main.py`

Most processing steps are currently commented out. Only data download steps are active: `download_lidar` / `download_aop_bbox`, `download_hyperspectral`, and `download_veg_structure_data`. Intermediate processing, ML training, and IC generation are commented out pending further development.
