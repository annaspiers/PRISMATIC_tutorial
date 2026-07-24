# Build and publish the PRISMATIC tutorial image for CyVerse VICE.
#
#   just login                          # authenticate to Harbor (once)
#   just build                          # build locally, then smoke-test before pushing
#   just push                           # push the built image
#   just publish                        # build + push in one step
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
