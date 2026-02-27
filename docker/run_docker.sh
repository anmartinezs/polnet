#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 --config <yaml> --out_dir <dir> [--data_dir <dir>] [--gpus] [-- <polnet args>]

Docker-level flags:
  --config  <yaml>   Required. Path to the YAML configuration file.
  --out_dir <dir>    Required. Host directory where output will be written.
  --data_dir <dir>   Optional. Host directory with input model files.
  --gpus             Optional. Pass --gpus all to docker for GPU support.

Anything after '--' is forwarded to the polnet CLI inside the container.
Examples:  -- -v          (INFO verbosity)
           -- -vv         (DEBUG verbosity)
           -- -s 12345    (override seed)
           -- -n 5        (generate 5 tomograms)
EOF
    exit 1
}

config_file=""
out_dir=""
data_dir=""
gpu_flag=""
polnet_args=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)   config_file="$2"; shift 2 ;;
        --out_dir)  out_dir="$2";     shift 2 ;;
        --data_dir) data_dir="$2";    shift 2 ;;
        --gpus)     gpu_flag="--gpus all"; shift ;;
        --)         shift; polnet_args=("$@"); break ;;
        *)          usage ;;
    esac
done

[[ -z "${config_file:-}" ]] && usage
[[ -z "${out_dir:-}" ]]     && usage

# Prevent Docker from creating the host output directory as root
mkdir -p "$out_dir"

mounts=(
    -v "$(realpath "$config_file")":/app/generation_config.yaml:ro
    -v "$(realpath "$out_dir")":/app/results
)
[[ -n "${data_dir:-}" ]] && mounts+=(-v "$(realpath "$data_dir")":/app/data:ro)

docker run --rm \
    -e PYTHONUNBUFFERED=1 \
    --user "$(id -u):$(id -g)" \
    ${gpu_flag} \
    "${mounts[@]}" \
    polnet_docker /app/generation_config.yaml "${polnet_args[@]}"