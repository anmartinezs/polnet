"""Synthetic tomogram lifecycle orchestrator.

Encapsulates the full single-tomogram pipeline: sample geometry
generation, TEM simulation, and structured output writing.  Used
as the primary execution unit in
:func:`polnet.pipeline.gen_tomos`.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import shutil
import time
from pathlib import Path

from .logging_conf import _LOGGER as logger
from .sample import (
    FlmsFile,
    MbFile,
    PmFile,
    PnFile,
    SyntheticSample,
)
from .tem import (
    TEM,
    TEMFile,
)
from .utils import lio


class SynthTomo:
    """A single synthetic tomogram.

    Owns the full lifecycle: sample generation -> TEM simulation -> output.

    Args:
        index: Tomogram index (1-based), used in output filenames.
        mbs_file_list: Relative paths to membrane .mbs config files.
        flms_file_list: Relative paths to helix .flms config files.
        pns_file_list: Relative paths to cytosolic protein .pns config files.
        pms_file_list: Relative paths to membrane protein .pms config files.
        tem_file_path: Absolute path to the .tem config file, or None to skip TEM.
        data_path: Absolute path to the input data directory.
        tomo_dir: Absolute path to this tomogram's output directory.
        shape: VOI shape (X, Y, Z) in voxels.
        v_size: Voxel size in Angstroms.
        offset: VOI offset in voxels.

    Raises:
        TypeError: If any of the arguments are of the wrong type or have invalid values.
    """

    def __init__(
        self,
        *,
        index: int,
        mbs_file_list: list,
        flms_file_list: list,
        pns_file_list: list,
        pms_file_list: list,
        tem_file_path: Path | None,
        data_path: Path,
        tomo_dir: Path,
        shape: tuple,
        v_size: float,
        offset: tuple,
    ):
        if tem_file_path is not None and not isinstance(tem_file_path, Path):
            raise TypeError("tem_file_path must be a Path object or None")
        if not isinstance(data_path, Path):
            raise TypeError("data_path must be a Path object")
        if not isinstance(tomo_dir, Path):
            raise TypeError("tomo_dir must be a Path object")
        if (
            not isinstance(shape, list)
            or not all(isinstance(s, int) for s in shape)
            or len(shape) != 3
        ):
            raise TypeError("shape must be a list of 3 integers.")
        if not isinstance(v_size, float) or not v_size > 0:
            raise TypeError("v_size must be a positive float")
        if not isinstance(offset, list) or len(offset) != 3:
            raise TypeError("offset must be a list of 3 integers")

        self._id = index
        self._mbs_files = mbs_file_list
        self._flms_files = flms_file_list
        self._pns_files = pns_file_list
        self._pms_files = pms_file_list
        self._tem_file_path = tem_file_path
        self._sample = None
        self._temic = None
        self._snr = None
        self._data_path = data_path
        self._tomo_dir = tomo_dir
        self._shape = shape
        self._v_size = v_size
        self._offset = offset

    def reset(self) -> None:
        """Reset the synthetic tomogram by clearing the sample and TEM simulation results.

        Returns:
            None
        """
        self._sample = None
        self._temic = None
        self._snr = None
        if self._tomo_dir.exists() and self._tomo_dir.is_dir():
            shutil.rmtree(self._tomo_dir)
            logger.info("Output directory %s cleared.", self._tomo_dir)

    def generate_sample(self) -> None:
        """Build the synthetic sample.

        Raises:
            RuntimeError: If the sample has already been generated.
        """

        if self._sample is not None:
            raise RuntimeError("Sample has already been generated.")

        logger.info("Generating synthetic sample.")

        tic = time.time()

        self._sample = SyntheticSample(
            shape=self._shape, v_size=self._v_size, offset=self._offset
        )

        logger.info("Adding membranes to the sample.")

        for mb_file_rpath in self._mbs_files:
            mb_file_apath = self._data_path / mb_file_rpath
            mb_file = MbFile()
            mb_params = mb_file.load(mb_file_apath)

            self._sample.add_set_membranes(params=mb_params, max_mbtries=10)

        logger.info(
            "Membranes added. Adding helicoidal networks to the sample."
        )

        for flms_file_rpath in self._flms_files:
            flms_file_apath = self._data_path / flms_file_rpath
            flms_file = FlmsFile()
            flms_params = flms_file.load(flms_file_apath)

            self._sample.add_helicoidal_network(params=flms_params)

        logger.info(
            "Helicoidal networks added. Adding cytosolic proteins to the sample."
        )

        for pn_file_rpath in self._pns_files:
            pn_file_apath = self._data_path / pn_file_rpath
            pn_file = PnFile()
            pn_params = pn_file.load(pn_file_apath)

            self._sample.add_set_cproteins(
                params=pn_params,
                data_path=self._data_path,
                surf_dec=0.9,
                mmer_tries=20,
                pmer_tries=100,
            )

        logger.info(
            "Cytosolic proteins added. Adding membrane-bound proteins to the sample."
        )

        for pm_file_rpath in self._pms_files:
            pm_file_apath = self._data_path / pm_file_rpath
            pm_file = PmFile()
            pm_params = pm_file.load(pm_file_apath)

            self._sample.add_set_mb_proteins(
                params=pm_params,
                data_path=self._data_path,
                surf_dec=0.9,
                mmer_tries=20,
                pmer_tries=100,
            )

        toc = time.time()
        logger.info("Synthetic sample generated in %.2f seconds.", toc - tic)

    def simulate_tem(self) -> None:
        """Simulate TEM imaging and reconstruct the tomogram.

        Requires that the sample has already been generated.
        Skips TEM simulation if no TEM config file was provided.
        """
        if self._sample is None:
            raise RuntimeError("Sample has not been generated yet.")
        if self._tem_file_path is None:
            logger.warning(
                "No TEM config file provided. Skipping TEM simulation."
            )
            return

        self._tomo_dir.mkdir(parents=True, exist_ok=True)

        self._temic = TEM(self._tomo_dir / "tem")
        tem_file = TEMFile()
        tem_params = tem_file.load(self._tem_file_path)

        self._snr = self._temic.simulate(
            density=self._sample.density,
            tem_params=tem_params,
            v_size=self._sample.v_size,
        )

        logger.info("TEM simulation complete (SNR=%s).", self._snr)

    def save(self) -> None:
        """Write all output files to tomo_dir.

        Outputs: density MRC, labels MRC, poly VTP, skel VTP, motif_list.csv.
        If TEM was run: micrographs MRC and reconstruction MRC (with SNR tag).
        """

        output_folder = self._tomo_dir
        if output_folder is None or not isinstance(output_folder, Path):
            raise TypeError("output_folder must be a Path object.")
        output_folder.mkdir(parents=True, exist_ok=True)

        self._save_motif_list()

        den_path = output_folder / f"tomo_{self._id:03d}_den.mrc"
        lio.write_mrc(
            self._sample.density,
            den_path,
            v_size=self._sample.v_size,
        )

        lbl_path = output_folder / f"tomo_{self._id:03d}_lbl.mrc"
        lio.write_mrc(
            self._sample.labels,
            lbl_path,
            v_size=self._sample.v_size,
        )

        if self._sample.poly_vtp is not None:
            poly_den_path = output_folder / f"tomo_{self._id:03d}_poly_den.vtp"
            lio.save_vtp(
                self._sample.poly_vtp,
                poly_den_path,
            )
        else:
            logger.warning("No poly_vtp data to save.")

        if self._sample.skel_vtp is not None:
            poly_skel_path = (
                output_folder / f"tomo_{self._id:03d}_poly_skel.vtp"
            )
            lio.save_vtp(
                self._sample.skel_vtp,
                poly_skel_path,
            )
        else:
            logger.warning("No skel_vtp data to save.")

        if self._temic is not None:
            snr_tag = f"_snr{self._snr}" if self._snr is not None else ""
            tem_dir = self._tomo_dir / "tem"

            mics_src = tem_dir / "out_micrographs.mrc"
            if mics_src.exists():
                mics_dst = (
                    output_folder / f"tomo_{self._id:03d}_mics{snr_tag}.mrc"
                )
                shutil.copyfile(mics_src, mics_dst)
            else:
                logger.warning(
                    "Micrographs file not found in TEM working directory."
                )

            rec_src = tem_dir / "out_rec3d.mrc"
            if rec_src.exists():
                rec_dst = (
                    output_folder / f"tomo_{self._id:03d}_rec{snr_tag}.mrc"
                )
                shutil.copyfile(rec_src, rec_dst)
            else:
                logger.warning(
                    "Reconstruction file not found in TEM working directory."
                )

    def _save_motif_list(self) -> None:
        """Save per-monomer ground truth to a tab-separated CSV file.

        Columns: Type, Label, Code, Polymer, X, Y, Z, Q1, Q2, Q3, Q4
        """
        motifs = self._sample.motifs
        if len(motifs) == 0:
            logger.warning("No motifs to save.")
            return

        out_file = self._tomo_dir / f"tomo_{self._id:03d}_motif_list.csv"
        header = "\t".join(motifs.dtype.names)

        with open(out_file, "w", newline="", encoding="utf-8") as f:
            f.write(header + "\n")
            for row in motifs:
                values = [
                    row["type"],
                    str(row["label"]),
                    row["code"],
                    str(row["polymer"]),
                    f"{row['x']:.8g}",
                    f"{row['y']:.8g}",
                    f"{row['z']:.8g}",
                    f"{row['q1']:.8g}",
                    f"{row['q2']:.8g}",
                    f"{row['q3']:.8g}",
                    f"{row['q4']:.8g}",
                ]
                f.write("\t".join(values) + "\n")

        logger.info(
            "Motif list (%d entries) saved to %s", len(motifs), out_file
        )

    def print_summary(self) -> None:
        """Show a summary of the synthetic sample.

        Returns:
            None
        """
        if self._sample is None:
            logger.error("No sample generated yet.")
            return

        self._sample.print_summary()
