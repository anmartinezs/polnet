# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-02-25

### Added

- **Curvatubes membrane generator** — GPU-accelerated phase-field membrane generation via the vendored `cvtub` library (Anna Song, MIT). Registered as the `"curvatubes"` membrane type. Supports Helfrich + polynomial curvature parameterisation with conserved and non-conserved gradient flows.
- **YAML configuration system** — replaced Python-script configs with declarative YAML files (`config/all_features.yaml`, `config/curvatubes.yaml`).
- **CLI entry point** (`polnet`) — new command-line interface with `-v`/`-vv` verbosity, `--seed`, `--ntomos`, and `--log-dir` options.
- **Centralised logging** — rotating file handler + console handler via `logging_conf.py`; replaces all ad-hoc `print()` statements across the codebase.
- **Per-monomer ground-truth export** — `motif_list.csv` with type, label, code, polymer index, position (x/y/z), and orientation quaternion (q1–q4) per monomer.
- **Global labels table** — `labels_table.csv` shared across all tomograms for consistent multi-tomogram labelling.
- **Membrane factory pattern** — decorator-based `MbFactory` registry for pluggable membrane generators.
- **Sphinx API documentation** — full autodoc RST coverage for all modules with Google-style docstring support via Napoleon.
- **Vendored `cvtub` package** under `src/external/cvtub/` with MIT license preserved. Header comments in every file identify the origin and license.
- **NOTICE file** for Apache-2.0 third-party attribution compliance.
- **Docker workflow** — rewritten `create_docker.sh` and `run_docker.sh` with volume-mount support for config, data, and output directories. IMOD 4.11.25 is downloaded automatically at build time.

### Changed

- **Complete package restructure** — flat module layout replaced with hierarchical sub-packages under `src/polnet/`: `sample/membranes/`, `sample/filaments/`, `sample/pns/`, `sample/pms/`, `sample/polymers/`, `tem/`, `utils/`.
- **Membrane generators refactored** — `MbGen` abstract base class with `_build()` + `generate_set()` separation; concrete generators (sphere, ellipsoid, toroid, curvatubes) inherit shared occupancy-loop logic.
- **`Mb` data class** — standalone membrane data holder with validation, properties, and `insert_density_svol()` method.
- **`SynthTomo` lifecycle** — clean keyword-only constructor with type validation; explicit `generate_sample()` → `simulate_tem()` → `save()` pipeline.
- **`SyntheticSample`** now tracks entity IDs, label registry, motifs, and per-type structure/voxel counts.
- **Config file parsers** (`MbFile`, `FlmsFile`, `PnFile`, `PmFile`, `TEMFile`) use `ast.literal_eval` for safe value parsing.
- **Getter methods → `@property`** throughout the codebase.
- **Docker image** rebuilt on `python:3.11-bookworm` (Debian 12) with IMOD 4.11.25 (downloaded at build time). Uses entrypoint wrapper to source IMOD environment. The IMOD installer filename is preserved during download (it derives its version from its own name).

### Fixed

- All `print()` calls replaced with structured `logging` calls.
- All bare `assert` statements replaced with proper exceptions (`ValueError`, `TypeError`, `RuntimeError`, etc.).
- VTK `SetValue` `TypeError` when writing integer arrays.
- **Comprehensive GPU VRAM cleanup** in curvatubes integration — 9 leak sources identified and fixed (nn.Modules, optimizer state, global tensors, autograd graph, closure reference cycles).
- MRC files opened in read-only mode where appropriate.
- Import organisation standardised (stdlib → third-party → local).

### Not Yet Migrated

- **Legacy scripts** (`scripts/`) — the v1.0.0 Python scripts have not been updated to use the new v1.1.0 API and YAML configuration system. They will be adapted in a future release.
- **Jupyter notebooks** (`gui/`) — the interactive notebooks for dataset generation, membrane alignment, model creation, and atomic-to-density conversion still depend on the v1.0.0 flat module layout and will not work with the current package structure. Updated notebooks will be provided in a future release.

## [1.0.0] - 2025-07-29

### Added

- Initial release of Polnet, a comprehensive tool for generating synthetic cryo-electron tomograms
- Core simulation capabilities for biological structures:
  - **Membranes**: Simulation of cellular membranes with spherical, ellipsoidal, and toroidal geometries
  - **Filaments**: Simulation of cytoskeletal structures including actin filaments and microtubules with helicoidal geometry
  - **Protein Complexes**: Simulation of globular protein clusters using Self-Avoiding Walk on Lattice (SAWLC) networks
  - **Membrane Proteins**: Simulation of membrane-bound protein complexes
- **Data Generation Pipeline**:
  - Configurable tomogram dimensions and voxel sizes
  - Adjustable occupancy and density parameters
  - Support for multiple biological structures in single tomograms
- **Output Formats**:
  - Ground truth density maps (`.mrc` format)
  - Segmentation label maps (`.mrc` format)  
  - 3D polydata models (`.vtp` format) for visualization
  - Simulated micrograph tilt-series for realistic cryo-ET data
- **TEM Simulation**:
  - Configurable tilt angles and detector parameters
  - Signal-to-noise ratio controls
  - Misalignment simulation for realistic data

### Technical Details

- Python-based implementation with scientific computing stack
- Support for Docker containerization
