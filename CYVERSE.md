# Running the PRISMATIC tutorial on CyVerse VICE

The files, `Dockerfile` `Justfile`, `jupyter_notebook_config.json`, and `.dockerignore` are
used to create a Docker image for [CyVerse VICE](https://cyverse.org/vice) (Visual Interactive
Computing Environment), so workshop participants launch a ready-to-run JupyterLab instance.

## Files

| File | Purpose |
|------|---------|
| `Dockerfile` | Builds the image on top of the CyVerse datascience base and bakes in the env + kernels. |
| `environment.yml` | The tutorial's conda environment (unchanged; the image builds from it directly). |
| `jupyter_notebook_config.json` | Disables JupyterLab's token/password auth. |
| `.dockerignore` | Keeps the build context small (only `environment.yml` and the Jupyter config are sent). |
| `justfile` | Build / push / publish recipes. |

## What's in the image

- **Base:** `harbor.cyverse.org/vice/jupyter/datascience:4.6.0` — a pinned CyVerse image with
  JupyterLab, GitHub CLI, GoCommands, and AI CLIs. Pinned (not `:latest`) so the environment
  can't drift between build day and workshop day.
- **The `prismatic_tutorial` conda env**, created from `environment.yml` at build time.
- **Two Jupyter kernels** registered from that env:
  - **PRISMATIC (Python)** — the Python 3.10 kernel with all the tutorial's packages.
  - **PRISMATIC (R)** — an R 4.4 kernel (via `r-irkernel`).
  - R is also available from Python through `rpy2`, which is how the tutorial's
    `initialize/*.R` helpers are driven.

### Notebook kernel

`prismatic_workshop.ipynb` is pinned to the `prismatic_tutorial` kernel in its metadata, so it
opens on the correct kernel automatically. If you author a new notebook, pick **PRISMATIC
(Python)** from the launcher or via *Kernel → Change Kernel*.

### Kernel install location

The Dockerfile installs the kernels with `--prefix /opt/conda` (into
`/opt/conda/share/jupyter/kernels`), **not** `--user`. In VICE a per-user persistent home is
mounted over `/home/jovyan` at runtime, which would hide anything installed under
`~/.local`. Installing into the image keeps the kernels visible for every user.

## Building and publishing

Prerequisites: [`docker`](https://docs.docker.com/get-docker/) with Buildx,
[`just`](https://github.com/casey/just), and push access to the Harbor project referenced in
the `justfile` (`harbor.cyverse.org/vice/prismatic-tutorial`).

```bash
just login          # authenticate to Harbor (once)
just build          # build the linux/amd64 image locally
# ...launch it and smoke-test (see below)...
just push           # push to Harbor

# or build + push in one step:
just publish

# override the tag for any recipe:
just tag=2026-07-workshop-v2 build
```

The image targets `linux/amd64` (what VICE runs).

### Smoke test before pushing

Run the image and open JupyterLab, then confirm:

- The launcher shows **PRISMATIC (Python)** and **PRISMATIC (R)** kernels
  (`jupyter kernelspec list` should show them under `/opt/conda/share/jupyter/kernels`).
- In the Python kernel: `import geopandas, rasterio, tensorflow, hydra` succeeds.
- In the R kernel: `library(sf); library(terra)` load and `R.version.string` reports 4.4.

## In VICE

VICE routes to JupyterLab on port 8888 with authentication handled at the ingress, so the image
ships with the JupyterLab token/password disabled (`jupyter_notebook_config.json`). Register the
published image as a VICE app pointing at
`harbor.cyverse.org/vice/prismatic-tutorial:<tag>`.

## Updating the environment

The image builds straight from `environment.yml`, so to change the environment, edit that file
and rebuild. Bump the image tag for a new workshop offering rather than overwriting an existing
tag, so past sessions stay reproducible.
