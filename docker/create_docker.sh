#!/usr/bin/env bash
set -euo pipefail

# Build the polnet Docker image.
# IMOD 4.11.25 is downloaded automatically during the build.
# Run from the project root.

echo "Building polnet_docker image..."
docker build -f docker/Dockerfile -t polnet_docker .
echo "Done. Run with: ./docker/run_docker.sh --config <yaml> --out_dir <dir>"