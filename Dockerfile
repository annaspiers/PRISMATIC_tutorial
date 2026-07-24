# PRISMATIC tutorial environment for CyVerse VICE.
# Base = CyVerse datascience image (JupyterLab + AI CLIs), pinned for reproducibility.
# Entrypoint, iRODS/gitconfig setup, and Jupyter config are inherited from the base image.
#
# Note: numbered tags (e.g. 4.6.0) do NOT bundle RStudio (only the moving :latest tag does).
# The workshop runs in JupyterLab, using R via the registered R kernel and rpy2, so RStudio
# is not needed here.
FROM harbor.cyverse.org/vice/jupyter/datascience:4.6.0

# ---- System libs occasionally needed to build geospatial / R packages (root) ----
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy entry script
COPY entry.sh /bin
RUN sudo chmod +x /bin/entry.sh

# ---- Conda env + kernels (jovyan owns /opt/conda in the Jupyter base image) ----
USER jovyan
WORKDIR /home/jovyan

# environment.yml sits beside this Dockerfile in the repo -> no vendoring needed.
COPY --chown=1000:100 environment.yml /home/jovyan/environment.yml
RUN mamba env create -f /home/jovyan/environment.yml && \
    mamba install -n prismatic_tutorial -y ipykernel r-irkernel && \
    mamba clean -afy

# Register the env's Python and R kernels with the already-installed JupyterLab.
# Run as jovyan so --user kernels land in jovyan's home where JupyterLab finds them.
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate prismatic_tutorial && \
    python -m ipykernel install --prefix /opt/conda --name prismatic_tutorial --display-name "PRISMATIC (Python)" \
        --env R_HOME /opt/conda/envs/prismatic_tutorial/lib/R \
        --env LD_LIBRARY_PATH /opt/conda/envs/prismatic_tutorial/lib/R/lib:/opt/conda/envs/prismatic_tutorial/lib \
        --env PATH /opt/conda/envs/prismatic_tutorial/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin && \
    R --quiet -e "IRkernel::installspec(name = 'prismatic_r', displayname = 'PRISMATIC (R)', prefix = '/opt/conda')"

# Default new terminals to the env.
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/jovyan/.bash_profile && \
    echo "conda activate prismatic_tutorial"   >> /home/jovyan/.bash_profile

# Copy the Jupyter configuration into the image.
COPY --chown=1000:100 jupyter_notebook_config.json /opt/conda/etc/jupyter/jupyter_notebook_config.json

# Copy the contents of this repository into the image.
COPY --chown=1000:100 . /home/jovyan/PRISMATIC_tutorial

# Restore the data-store working dir used by VICE.
WORKDIR /home/jovyan/data-store
