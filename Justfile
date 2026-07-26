# Build and publish the PRISMATIC tutorial image for CyVerse VICE.
#
#   just login                          # authenticate to Harbor (once)
#   just build                          # build locally, then smoke-test before pushing
#   just push                           # push the built image
#   just publish                        # build + push in one step
#   just test-notebook                  # run the notebook against the built image (no rebuild)
#   just clean                          # remove local test artifacts (data/, outputs/)
#   just tag=2026-07-workshop-v2 build  # override the tag for any recipe
#
# VICE runs linux/amd64, so the platform is pinned. On Apple Silicon the amd64
# build runs under emulation (correct for VICE, just slower).

image    := "harbor.cyverse.org/vice/prismatic-tutorial"
tag      := "2026-07-workshop"
platform := "linux/amd64"

# Build the image locally (loads it into your local Docker for smoke-testing).
build:
    docker buildx build --platform {{platform}} -t {{image}}:{{tag}} --load .

# Push a previously built image.
push:
    docker push {{image}}:{{tag}}

# Build for amd64 and push in one step (skips the local smoke-test).
publish:
    docker buildx build --platform {{platform}} -t {{image}}:{{tag}} --push .

# Authenticate to Harbor (needed before pushing).
login:
    docker login harbor.cyverse.org

# Run the workshop notebook headlessly against the built image, with your local
# working copy mounted in — so code edits are testable in seconds with no rebuild.
# Downloaded NEON data lands in ./data (gitignored) and persists between runs, so
# re-runs skip already-completed steps via the notebook's force_rerun cache.
#   just test-notebook            # run every cell
#   just test-notebook 2          # run only cells 0..2 (e.g. the downloads)
#   just test-notebook 10 120     # cap each cell at 120s to confirm heavy cells start
# NEON now requires an API token; export it in your shell and it is forwarded in:
#   NEON_API_TOKEN=... just test-notebook 1
test-notebook to="1000" timeout="0":
    docker run --rm \
      -v {{justfile_directory()}}:/home/jovyan/data-store/PRISMATIC_tutorial \
      -e R_HOME=/opt/conda/envs/prismatic_tutorial/lib/R \
      -e LD_LIBRARY_PATH=/opt/conda/envs/prismatic_tutorial/lib/R/lib:/opt/conda/envs/prismatic_tutorial/lib \
      -e NEON_API_TOKEN \
      --entrypoint bash {{image}}:{{tag}} \
      -lc 'cd "$TUTORIAL_HOME" && python tools/run_cells.py prismatic_workshop.ipynb --to {{to}} --cell-timeout {{timeout}}'

# Remove local test artifacts (downloaded NEON data and pipeline outputs).
clean:
    rm -rf {{justfile_directory()}}/data {{justfile_directory()}}/outputs
