"""Tomogram generation pipeline.

Orchestrates the creation of synthetic Cryo-ET tomograms by parsing
a YAML configuration, constructing SynthTomo instances, and driving
the sample generation → TEM simulation → output sequence.

Called by :func:`polnet.cli.app`. Can also be used programmatically::

    from polnet.pipeline import gen_tomos
    gen_tomos(config_dict)

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import csv
import random
from pathlib import Path

import numpy as np

from .logging_conf import _LOGGER as logger
from .tomogram import SynthTomo


def _save_labels_table(
    out_dir: Path,
    membranes: list[str],
    filaments: list[str],
    proteins: list[str],
    mb_proteins: list[str],
) -> None:
    """Write a global labels_table.csv mapping every model file to a static label.

    The label counter starts at 1 and increments by 1 for each model entry,
    in order: membranes → filaments → proteins → membrane proteins.
    This file is shared across all tomograms to guarantee consistent labelling.

    Args:
        out_dir: Root output directory.
        membranes: Relative paths to .mbs config files.
        filaments: Relative paths to .flms config files.
        proteins: Relative paths to .pns config files.
        mb_proteins: Relative paths to .pms config files.
    """
    out_file = out_dir / "labels_table.csv"
    label = 1
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["MODEL", "LABEL"], delimiter="\t"
        )
        writer.writeheader()
        for model_list in (membranes, filaments, proteins, mb_proteins):
            for model_name in model_list:
                writer.writerow({"MODEL": model_name, "LABEL": label})
                label += 1
    logger.info("Labels table (%d entries) saved to %s", label - 1, out_file)


def gen_tomos(config: dict) -> None:
    """Generate one or more synthetic tomograms from a configuration dict.

    Args:
        config: Parsed YAML configuration with required keys:
            - folders.root (str | None): Project root; None → inferred from package location.
            - folders.input (str): Relative path to input data directory.
            - folders.output (str): Relative path to output directory.
            - global.seed (int | None): RNG seed for reproducibility.
            - global.ntomos (int): Number of tomograms to generate.
            - sample.voi_shape (list[int]): Volume shape [X, Y, Z] in voxels.
            - sample.voi_offset (list[int]): VOI offset in voxels.
            - sample.v_size (float): Voxel size in Angstroms.
            - sample.membranes (list[str]): Relative paths to .mbs config files.
            - sample.filaments (list[str]): Relative paths to .flms config files.
            - sample.pns (list[str]): Relative paths to .pns config files.
            - sample.pms (list[str]): Relative paths to .pms config files.
            - tem.config (str | None): Relative path to .tem config file.

    Raises
        KeyError: If required config keys are missing.
        FileNotFoundError: If the input data directory does not exist.
        FileNotFoundError: If the TEM config file does not exist.
    """

    logger.info("Starting tomogram generation pipeline.")

    logger.info("Parsing configuration and validating input paths.")

    if "folders" not in config:
        raise KeyError("Config missing 'folders' section.")
    root_dir = Path(config["folders"].get("root") or Path.cwd())

    data_rpath = config["folders"].get("input", None)
    if data_rpath is None:
        raise KeyError("Config missing 'input' path in 'folders' section.")
    data_apath = Path(root_dir / data_rpath)
    if not data_apath.exists():
        raise FileNotFoundError(f"Data path {data_apath} does not exist.")

    out_rpath = config["folders"].get("output", None)
    if out_rpath is None:
        raise KeyError("Config missing 'output' path in 'folders' section.")
    out_apath = Path(root_dir / out_rpath)
    if not out_apath.exists():
        out_apath.mkdir(parents=True, exist_ok=True)

    seed = config["global"].get("seed", None)
    n_tomos = config["global"].get("ntomos", None)
    if n_tomos is None:
        raise KeyError("Config missing 'ntomos' in 'global' section.")

    voi_shape = config["sample"].get("voi_shape", None)
    if voi_shape is None:
        raise KeyError("Config missing 'voi_shape' in 'sample' section.")

    voi_offs = config["sample"].get("voi_offset", None)
    if voi_offs is None:
        raise KeyError("Config missing 'voi_offset' in 'sample' section.")

    v_size = config["sample"].get("v_size", None)
    if v_size is None:
        raise KeyError("Config missing 'v_size' in 'sample' section.")

    membranes = config["sample"].get("membranes", [])
    filaments = config["sample"].get("filaments", [])
    proteins = config["sample"].get("pns", [])
    mb_proteins = config["sample"].get("pms", [])

    tem_file_rpath = config.get("tem", {}).get("config", None)
    tem_apath = data_apath / tem_file_rpath if tem_file_rpath else None
    if tem_apath is not None and not tem_apath.exists():
        raise FileNotFoundError(f"TEM config file {tem_apath} does not exist.")

    logger.info("Setting random seed to %s for reproducibility.", seed)
    random.seed(seed)
    np.random.seed(seed)

    logger.info(
        "Saving labels table at %s and generating %d tomogram(s) "
        "with VOI shape %s, offset %s, and voxel size %s \u00c5.",
        out_apath,
        n_tomos,
        voi_shape,
        voi_offs,
        v_size,
    )
    _save_labels_table(out_apath, membranes, filaments, proteins, mb_proteins)

    for tomo_id in range(n_tomos):
        logger.info("Generating tomogram %d/%d.", tomo_id + 1, n_tomos)
        synth_tomo = SynthTomo(
            index=tomo_id + 1,
            mbs_file_list=membranes,
            flms_file_list=filaments,
            pns_file_list=proteins,
            pms_file_list=mb_proteins,
            tem_file_path=tem_apath,
            data_path=data_apath,
            tomo_dir=out_apath / f"Tomo{tomo_id + 1:03d}",
            shape=voi_shape,
            v_size=v_size,
            offset=voi_offs,
        )

        logger.info(
            "Generating sample for tomogram %d/%d.", tomo_id + 1, n_tomos
        )
        synth_tomo.generate_sample()

        logger.info("Simulating TEM for tomogram %d/%d.", tomo_id + 1, n_tomos)
        synth_tomo.simulate_tem()

        logger.info("Saving outputs for tomogram %d/%d.", tomo_id + 1, n_tomos)
        synth_tomo.save()
        synth_tomo.print_summary()

        logger.info(
            "Tomogram %d/%d generation complete.", tomo_id + 1, n_tomos
        )
