# PolNet

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://github.com/anmartinezs/polnet/blob/main/LICENSE)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/anmartinezs)

Python package for generating synthetic datasets of the cellular context for Cryo-Electron Tomography.

## Table of Contents

- [PolNet](#polnet)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Using conda (recommended)](#using-conda-recommended)
    - [Using venv](#using-venv)
    - [For developers](#for-developers)
      - [Available extras](#available-extras)
      - [Code style](#code-style)
    - [Verifying the installation](#verifying-the-installation)
  - [Usage](#usage)
    - [CLI reference](#cli-reference)
    - [Configuration](#configuration)
    - [Examples](#examples)
  - [Output](#output)
  - [Docker](#docker)
    - [Prerequisites](#prerequisites)
    - [Building the image](#building-the-image)
    - [Running a container](#running-a-container)
    - [Using the built-in example](#using-the-built-in-example)
    - [Cross-platform support](#cross-platform-support)
    - [Notes](#notes)
  - [Package description](#package-description)
  - [Third-party code](#third-party-code)
    - [Curvatubes performance notes](#curvatubes-performance-notes)
  - [Documentation](#documentation)
  - [Citations](#citations)

## Overview

PolNet simulates the cellular context found in cryo-electron tomograms. The pipeline is split into two main stages:

1. **Sample simulation** — generates a synthetic 3-D volume of interest (VOI) populated with membranes, filaments, cytosolic proteins, and membrane-bound proteins.
2. **Microscope simulation & tomogram reconstruction** — applies a TEM acquisition model (tilt-series projection, detector noise, misalignment) and reconstructs the final tomogram using IMOD tools.

All simulation parameters are controlled through a single YAML configuration file. The `polnet` CLI is the main entry point.

## Requirements

- **Python ≥ 3.10**
- **IMOD** must be installed and available on `$PATH`. PolNet calls several IMOD commands (`xyzproj`, `tilt`, `alterheader`) during the TEM simulation stage.
  Installation guide: https://bio3d.colorado.edu/imod/doc/guide.html
- **Git** for cloning the repository.

## Installation

### Using conda (recommended)

```console
git clone https://github.com/anmartinezs/polnet.git
cd polnet

conda create -n polnet python=3.11 pip
conda activate polnet

pip install .
```

> **Note:** PyTorch and nibabel (required by the curvatubes membrane generator) are
> core dependencies and installed automatically. If you have a CUDA-capable GPU,
> make sure the PyTorch version matches your CUDA driver. See
> [pytorch.org/get-started](https://pytorch.org/get-started/locally/) for details.

### Using venv

```console
git clone https://github.com/anmartinezs/polnet.git
cd polnet

python3 -m venv .venv
source .venv/bin/activate    # Linux / macOS
# .venv\Scripts\activate     # Windows

pip install .
```

### For developers

Install the package in editable mode with the development extras:

```console
pip install -e ".[dev]"
```

This installs the core package plus **black**, **pylint**, and **pre-commit**.

To also install documentation and GUI dependencies:

```console
pip install -e ".[all]"        # dev + docs + gui
```

#### Available extras

| Extra | Contents |
|-------|----------|
| `dev` | black, pylint, pre-commit |
| `docs` | sphinx, sphinx-rtd-theme, sphinx-autodoc-typehints |
| `gui` | jupyter, ipyfilechooser, ipywidgets, tqdm, wget |
| `all` | dev + docs + gui |

#### Code style

- **Formatter:** black (line-length 79)
- **Linter:** pylint (default rules)
- **Docstring style:** Google style

Configuration for both tools is in `pyproject.toml`.

### Verifying the installation

After installation, the `polnet` command should be available:

```console
polnet --version
```

## Usage

### CLI reference

```
polnet <config.yaml> [options]
```

| Flag | Description |
|------|-------------|
| `<config.yaml>` | Path to the YAML configuration file (required). |
| `-v` / `-vv` | Increase console verbosity: `-v` for INFO, `-vv` for DEBUG. Default is WARNING. |
| `-o`, `--log-dir <dir>` | Directory for log files. Defaults to the output folder from the config. |
| `-s`, `--seed <int>` | Override the random seed defined in the config. |
| `-n`, `--ntomos <int>` | Override the number of tomograms to generate. |
| `--version` | Print version and exit. |

### Configuration

The pipeline is driven by a YAML file. An example is provided at `config/all_features.yaml`. The top-level sections are:

```yaml
metadata:       # Version, author, date (informational).
folders:        # root, input, and output paths.
global:         # ntomos, seed.
sample:         # VOI shape, voxel size, membranes, filaments, pns, pms.
tem:            # TEM configuration file path.
```

All structural model files (`.mbs`, `.flms`, `.pns`, `.pms`, `.tem`) referenced in the config are resolved relative to the `folders.input` directory. Default models are included in `data/`.

### Examples

```console
# Generate one tomogram with default settings (WARNING on console)
polnet config/all_features.yaml

# Verbose output (INFO level)
polnet config/all_features.yaml -v

# Debug output, custom seed, 5 tomograms
polnet config/all_features.yaml -vv -s 12345 -n 5

# Custom log directory
polnet config/all_features.yaml --log-dir /tmp/polnet_logs
```

## Output

After running `polnet config.yaml`, the output directory (set by `folders.output` in the config) will contain one sub-folder per tomogram plus a shared labels table:

```
results/                              # folders.output
├── labels_table.csv                  # Global mapping: model file → integer label
└── Tomo001/                          # Per-tomogram directory (Tomo001, Tomo002, …)
    ├── tomo_001_den.mrc              # Ground-truth density volume
    ├── tomo_001_lbl.mrc              # Integer label volume (segmentation ground truth)
    ├── tomo_001_poly_den.vtp         # VTK PolyData surfaces of all placed objects
    ├── tomo_001_poly_skel.vtp        # VTK PolyData skeleton / centrelines
    ├── tomo_001_motif_list.csv       # Per-monomer ground truth (CSV)
    ├── tomo_001_mics_snr<X>.mrc      # Simulated tilt-series micrographs (if TEM enabled)
    ├── tomo_001_rec_snr<X>.mrc       # Reconstructed tomogram (if TEM enabled)
    └── tem/                          # TEM working directory (intermediate files)
```

| File | Format | Description |
|------|--------|-------------|
| `labels_table.csv` | TSV | Maps every input model file to a unique integer label, shared across all tomograms. |
| `*_den.mrc` | MRC | Noise-free 3-D density volume at the configured voxel size. |
| `*_lbl.mrc` | MRC | Integer label volume — each voxel carries the label of the object that occupies it (0 = background). |
| `*_poly_den.vtp` | VTK | Triangle-mesh surfaces of all placed macromolecules and membranes, with associated scalar fields. |
| `*_poly_skel.vtp` | VTK | Skeleton / centreline polydata (filament backbones, polymer chains). |
| `*_motif_list.csv` | TSV | One row per placed monomer: type, label, PDB code, polymer index, position (X, Y, Z), orientation quaternion (Q1–Q4). |
| `*_mics_snr<X>.mrc` | MRC | Simulated tilt-series stack at the configured SNR. Only written when TEM simulation is enabled. |
| `*_rec_snr<X>.mrc` | MRC | Weighted back-projection reconstruction from the tilt series. Only written when TEM simulation is enabled. |

## Docker

The Docker image packages polnet and IMOD into a self-contained container so you can generate tomograms without installing either on your host.

### Prerequisites

- **Docker** must be installed and running.
- Internet access during build (IMOD 4.11.25 is downloaded automatically).

### Building the image

From the project root:

```console
./docker/create_docker.sh
```

This builds the `polnet_docker` image with IMOD 4.11.25 and all bundled data.

### Running a container

Use the helper script:

```console
./docker/run_docker.sh \
    --config config/all_features.yaml \
    --out_dir ./results
```

| Flag | Description |
|------|-------------|
| `--config <yaml>` | **Required.** Path to your YAML configuration file. |
| `--out_dir <dir>` | **Required.** Host directory where output will be written. |
| `--data_dir <dir>` | Optional. Host directory with input model files (`.mbs`, `.pns`, etc.). Mounted read-only at `/app/data` inside the container, overriding the built-in data. |
| `--gpus` | Optional. Pass `--gpus all` to Docker for GPU support (curvatubes membranes). |
| `-- <args>` | Optional. Everything after `--` is forwarded to the `polnet` CLI inside the container. |

The container runs as your current user (`uid:gid`) so output files are owned by you.

**Forwarding polnet CLI flags:**

```console
# Verbose output (INFO)
./docker/run_docker.sh --config config/all_features.yaml --out_dir ./results -- -v

# Debug output + custom seed + 5 tomograms
./docker/run_docker.sh --config config/all_features.yaml --out_dir ./results -- -vv -s 12345 -n 5

# GPU-accelerated curvatubes membranes
./docker/run_docker.sh --config config/all_features.yaml --out_dir ./results --gpus -- -v
```

### Using the built-in example

The image ships with `config/all_features.yaml` and all required input data. To run it directly:

```console
mkdir -p results
docker run --rm \
    -v "$(pwd)/results":/app/outdir \
    --user "$(id -u):$(id -g)" \
    polnet_docker
```

### Cross-platform support

The Docker **image** is Linux-based and works on all platforms via [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Linux, macOS, Windows with WSL 2).

The helper **scripts** (`create_docker.sh`, `run_docker.sh`) require a Bash shell. On Windows you can run them from:
- **WSL 2** (recommended) — the scripts work as-is.
- **Git Bash** — ships with [Git for Windows](https://gitforwindows.org/).

Alternatively, run the `docker build` and `docker run` commands shown above directly from PowerShell or Command Prompt.

### Notes

- The IMOD self-extracting `.sh` installer derives its version from its own filename (via a `sed` on `imod_[0-9.]*`). The Dockerfile preserves the original filename during download; renaming it would break the version extraction and cause a tar path mismatch at install time.
- The build verifies that the three required IMOD binaries (`xyzproj`, `tilt`, `alterheader`) are present and executable before finalising the image.
- PyTorch and nibabel (required by curvatubes membranes) are installed by default. To use GPU-accelerated membrane generation, pass `--gpus` to `run_docker.sh` (requires [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) on the host).

## Package description

```
polnet/                     Project root
├── src/polnet/             Core Python package
│   ├── cli.py              CLI entry point (polnet command)
│   ├── pipeline.py         Top-level orchestrator (gen_tomos)
│   ├── tomogram.py         SynthTomo: single tomogram lifecycle
│   ├── logging_conf.py     Logger setup (console + rotating file)
│   ├── sample/             Sample simulation
│   │   ├── sample.py       SyntheticSample: VOI management
│   │   ├── membranes/      Membrane generators (sphere, ellipsoid, toroid, curvatubes)
│   │   ├── filaments/      Filament generators (microtubules, actin)
│   │   ├── polymers/       Base polymer classes (Monomer, Polymer, Network)
│   │   ├── pns/            Cytosolic protein generators
│   │   └── pms/            Membrane-bound protein generators
│   ├── tem/                TEM simulation and reconstruction
│   │   ├── tem.py          TEM acquisition model (projection, noise, reconstruction)
│   │   └── tem_file.py     TEM parameter file parser
│   └── utils/              Shared utilities
│       ├── affine.py       Affine transformations and rotation utilities
│       ├── lio.py          File I/O (MRC, VTK, CSV)
│       ├── poly.py         Polygon and surface operations
│       └── utils.py        General-purpose helpers
├── src/external/           Vendored third-party packages (see Third-party code)
│   └── cvtub/              curvatubes library (MIT, Anna Song)
├── config/                 Example YAML configuration files
├── data/                   Input structural models and templates
│   ├── in_10A/             Macromolecular models at 10 Å voxel size
│   ├── in_5A/              Macromolecular models at 5 Å voxel size
│   ├── in_8A/              Macromolecular models at 8 Å voxel size
│   ├── in_flms/            Filament configuration files (.flms)
│   ├── in_mbs/             Membrane configuration files (.mbs)
│   ├── in_tem/             TEM configuration files (.tem)
│   └── templates/          Atomic models (PDB) and density maps (MRC)
├── docker/                 Docker configuration
└── docs/                   Documentation
    ├── source/             Sphinx documentation sources
    ├── build/              Sphinx-generated output (not committed)
    └── molecules_table.md  Descriptions of PDB models used
```

## Third-party code

PolNet vendors certain third-party libraries under `src/external/` to keep them bundled with the package while clearly separated from the PolNet source code. Each vendored package retains its original license.

| Package | Location | License | Copyright | Upstream |
|---------|----------|---------|-----------|----------|
| **cvtub** (curvatubes) | `src/external/cvtub/` | MIT | © 2021 Anna Song | [github.com/annasongmaths/curvatubes](https://github.com/annasongmaths/curvatubes) |

The full license text for each vendored package is included in its directory (e.g. `src/external/cvtub/LICENSE`). Every source file carries a header comment identifying its origin and license. Modifications made to vendored code are:

1. **Relative imports** — absolute `from cvtub.xxx` imports converted to relative `from .xxx` so the package works within the PolNet `src/` layout. Each changed line is marked with a `# polnet: relative import` comment.
2. **GPU VRAM cleanup** (`generator.py` only) — a comprehensive 7-step cleanup block was added at the end of `_generate_shape()` to prevent VRAM leaks when generating multiple surfaces in sequence. The block is marked with `# polnet: comprehensive GPU VRAM cleanup`.
3. **Code formatting** — all `.py` files reformatted with `black` (line-length 79) for consistency with the rest of the polnet codebase. No functional changes.

### Curvatubes performance notes

The curvatubes membrane generator runs a GPU-accelerated phase-field optimisation that can be **very slow** for large tomogram shapes. Generation times for a single membrane on a `300×300×300` grid are typically in the range of **30 minutes to several hours**, depending on the GPU, `CT_MAXEVAL`, and the curvature coefficients. A detailed analysis of computation time vs. tomogram shape is provided in [Seghiri, Gallego Nicolás et al. (2026)](https://doi.org/10.64898/2026.01.15.699326).

**GPU VRAM** is another limiting factor. The optimisation allocates multiple full-volume tensors (phase field, gradient, Hessian, FFT filters, optimizer state) on the GPU simultaneously. Empirical measurements on $S \times S \times 250$ volumes follow the quadratic relationship:

$$\text{VRAM (MB)} \approx 0.05664 \cdot S^{2}$$

| Tomogram shape ($S \times S \times 250$) | Predicted VRAM |
|------------------------------------------|----------------|
| $S = 128$  | ~0.9 GB  |
| $S = 200$  | ~2.2 GB  |
| $S = 300$  | ~5.0 GB  |
| $S = 400$  | ~8.8 GB  |
| $S = 500$  | ~13.8 GB |
| $S = 600$  | ~19.9 GB |

If VRAM is insufficient, PyTorch will raise an `OutOfMemoryError`. You can control the trade-off between quality and resource consumption with:

- **`CT_MAXEVAL`** — maximum number of energy evaluations (recommended: 6000).
- **`CT_INFORM_EVERY`** — log progress every N evaluations (recommended: 500; set to 0 to disable).

Both are set in the `.mbs` configuration file (see `data/in_mbs/curvatubes.mbs` for an annotated example).

## Documentation

Full API documentation is generated with [Sphinx](https://www.sphinx-doc.org/) from Google-style docstrings. To build it locally:

```console
pip install -e ".[docs]"
cd docs
make html          # Linux / macOS
# make.bat html     # Windows
```

The generated site will be at `docs/build/html/index.html`.

Additional resources:

- **Molecules table:** `docs/molecules_table.md` — detailed descriptions of the PDB models used to create macromolecular structural units provided in `data/`.

## Citations

Original article:

 - Martinez-Sanchez A., Lamm L., Jasnin M. & Phelippeau H. (2024) "Simulating the cellular context in synthetic datasets for cryo-electron tomography" *IEEE Transactions on Medical Imaging* 43(11), 3742-3754. [10.1109/TMI.2024.3398401](https://doi.org/10.1109/TMI.2024.3398401)

Realistic membranes extension:

 - Seghiri, R., Gallego Nicolás, J. D., Brandt, R., Meroño, M. A., Lefevre, P., Hajarolasvadi, N., Baum, D., Heebner, J., Phelippeau, H., Doux, P., & Martinez-Sanchez, A. (2026). TomoSegNet: Augmented membrane segmentation for cryo-electron tomography by simulating the cellular context. bioRxiv. [10.64898/2026.01.15.699326](https://doi.org/10.64898/2026.01.15.699326)

Curvatubes article (third party code):

 - Song, A. (2022). Generation of Tubular and Membranous Shape Textures with Curvature Functionals. *Journal of Mathematical Imaging and Vision* 64, 17–40. [10.1007/s10851-021-01049-9](https://doi.org/10.1007/s10851-021-01049-9)