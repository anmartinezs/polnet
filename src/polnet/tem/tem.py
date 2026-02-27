"""Transmission Electron Microscopy (TEM) simulation models.

Wraps IMOD command-line binaries (``xyzproj``, ``tilt``,
``alterheader``) to project and reconstruct tomographic volumes.

.. warning::
    IMOD must be installed and accessible on the system ``PATH``.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import math
import random
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import shift as sp_shift

from ..logging_conf import _LOGGER as logger
from ..utils import lio

## IMOD commands
IMOD_CMD_XYZPROJ = "xyzproj"
IMOD_CMD_TILT = "tilt"
IMOD_CMD_AHEADER = "alterheader"


class TEM:
    """Wrapper around IMOD command-line tools for TEM simulation.

    Manages a working directory for intermediate MRC files and
    tilt-angle lists, and exposes methods to generate a tilt
    series, reconstruct a 3-D volume, add noise, and adjust
    headers.

    Note:
        IMOD must be installed and on the system PATH.
    """

    def __init__(self, work_dir):
        """Initialise a TEM simulator.

        Args:
            work_dir (str | Path): Directory for intermediate
                files (created if it does not exist).
        """
        self.__work_dir = Path(work_dir)
        self.__log_file = self.__work_dir / "TEM.log"
        self.__create_work_dir()
        self.__vol_file = self.__work_dir / "in_vol.mrc"
        self.__micgraphs_file = self.__work_dir / "out_micrographs.mrc"
        self.__tangs_file = self.__work_dir / "out_tangs.tlt"
        self.__rec3d_file = self.__work_dir / "out_rec3d.mrc"

    def __create_work_dir(self):
        """Create the working directory if it does not already exist."""
        self.__work_dir.mkdir(parents=True, exist_ok=True)

    def __save_tangs_file(self, angs):
        """Write tilt angles to an IMOD-format .tlt file.

        Args:
            angs (iterable): Sequence of tilt angles in degrees.
        """
        with open(self.__tangs_file, "w", encoding="utf-8") as file:
            for ang in angs:
                file.write(str(ang) + "\n")

    def __load_tangs_file(self):
        """Read the tilt angles .tlt file into a numpy array.

        Returns:
            numpy.ndarray: Array of tilt angles in degrees.
        """
        angs = list()
        with open(self.__tangs_file, "r", encoding="utf-8") as file:
            for line in file:
                angs.append(float(line.strip()))
        return np.asarray(angs)

    def gen_tilt_series_imod(self, vol, angs, ax="X", mode="real"):
        """Generate a 2-D tilt series from a 3-D volume via IMOD xyzproj.

        Writes the volume as MRC, calls ``xyzproj``, and saves the
        tilt angles to the working directory.

        Args:
            vol (numpy.ndarray): Input 3-D tomogram.
            angs (iterable | range): Tilt angles in degrees.
            ax (str): Tilt axis: 'X' (default), 'Y', or 'Z'.
            mode (str): Output MRC mode: 'byte', 'int', or 'real'
                (default).

        Raises:
            TypeError: If vol is not a 3-D numpy array.
            ValueError: If angs is empty, or ax/mode is invalid.
            subprocess.CalledProcessError: If xyzproj fails.
            IOError: If the log file cannot be written.
        """

        if not isinstance(vol, np.ndarray) or vol.ndim != 3:
            raise TypeError("vol must be a 3-D numpy array.")
        if not hasattr(angs, "__len__") or len(angs) == 0:
            raise ValueError("angs must be a non-empty iterable.")
        if ax not in ("X", "Y", "Z"):
            raise ValueError("ax must be 'X', 'Y' or 'Z'.")
        if mode not in ("byte", "int", "real"):
            raise ValueError("mode must be 'byte', 'int' or 'real'.")

        # Call to IMOD binary (xyzproj)

        xyzproj_cmd = [IMOD_CMD_XYZPROJ]
        in_vol_path = self.__vol_file
        lio.write_mrc(vol, in_vol_path)
        xyzproj_cmd += ["-inp", str(in_vol_path)]
        out_vol_path = self.__micgraphs_file
        xyzproj_cmd += ["-o", str(out_vol_path)]
        xyzproj_cmd += ["-ax", ax]
        if isinstance(angs, range):
            xyzproj_cmd += ["-an"]
            tangles = (
                str(angs.start) + "," + str(angs.stop) + "," + str(angs.step)
            )
        else:
            xyzproj_cmd += ["-ta"]
            tangles = str(angs[0])
            for i in range(1, len(angs)):
                tangles += "," + str(angs[i])
        xyzproj_cmd += [tangles]
        if mode == "byte":
            xyzproj_cmd += ["-m", "0"]
        elif mode == "int":
            xyzproj_cmd += ["-m", "1"]
        elif mode == "real":
            xyzproj_cmd += ["-m", "2"]

        try:
            with open(self.__log_file, "a", encoding="utf-8") as file_log:
                file_log.write(
                    "\n["
                    + time.strftime("%c")
                    + "]RUNNING COMMAND:-> "
                    + " ".join(xyzproj_cmd)
                    + "\n"
                )
                subprocess.check_call(
                    xyzproj_cmd, stdout=file_log, stderr=file_log
                )
            self.__save_tangs_file(angs)
        except subprocess.CalledProcessError:
            logger.error("Error calling the command: %s", xyzproj_cmd)
            raise
        except IOError:
            logger.error("Log file could not be written: %s", self.__log_file)
            raise

    def add_detector_noise(self, snr):
        """Add Gaussian readout noise to the simulated micrographs.

        Dark-current noise is one order of magnitude lower than
        readout noise and is neglected.  The noise level is
        determined by the desired SNR relative to the mean
        foreground intensity.

        Args:
            snr (float): Target signal-to-noise ratio (linear
                scale); must be > 0.

        Raises:
            FileNotFoundError: If the micrograph MRC file is
                missing.
            ValueError: If snr <= 0.
        """

        if not self.__micgraphs_file.exists():
            raise FileNotFoundError(
                f"Micrographs file not found: {self.__micgraphs_file}"
            )
        if snr <= 0:
            raise ValueError("snr must be greater than 0.")

        mics = lio.load_mrc(self.__micgraphs_file)

        rng = np.random.default_rng()
        # Adding readout noise with Gaussian distribution to every micrograph
        for i in range(mics.shape[2]):
            mic = mics[:, :, i]
            mask = mic > 0
            mn = mic[mask].mean()
            sg_fg = mn / snr
            mics[:, :, i] = mic + rng.normal(mn, sg_fg, mic.shape)

        lio.write_mrc(mics, self.__micgraphs_file)

    def recon3D_imod(self, thick=None):
        """Reconstruct a 3-D tomogram from the tilt series via IMOD tilt.

        Reads the micrograph stack and tilt-angle file from the
        working directory, calls ``tilt``, and saves the result
        after swapping and flipping the axes to match the input
        coordinate system.

        Args:
            thick (int, optional): Output tomogram thickness along
                Z in voxels; None (default) uses the input volume
                X dimension.

        Raises:
            ValueError: If thick <= 0.
            subprocess.CalledProcessError: If tilt fails.
            IOError: If the log file cannot be written.
        """

        tilt_cmd = [IMOD_CMD_TILT]
        vol = lio.load_mrc(self.__vol_file, mmap=True, no_saxes=False)
        tilt_cmd += ["-inp", str(self.__micgraphs_file)]
        tilt_cmd += ["-output", str(self.__rec3d_file)]
        tilt_cmd += ["-TILTFILE", str(self.__tangs_file)]
        if thick is None:
            tilt_cmd += ["-THICKNESS", str(vol.shape[0])]
        else:
            if thick <= 0:
                raise ValueError("thick must be greater than 0.")
            tilt_cmd += ["-THICKNESS", str(thick)]

        try:
            with open(self.__log_file, "a", encoding="utf-8") as file_log:
                file_log.write(
                    "\n["
                    + time.strftime("%c")
                    + "]RUNNING COMMAND:-> "
                    + " ".join(tilt_cmd)
                    + "\n"
                )
                subprocess.check_call(
                    tilt_cmd, stdout=file_log, stderr=file_log
                )
        except subprocess.CalledProcessError:
            logger.error("Error calling the command: %s", tilt_cmd)
            raise
        except IOError:
            logger.error("Log file could not be written: %s", self.__log_file)
            raise

        # Swap Y-Z axes from the output given by IMOD
        hold_rec = np.swapaxes(lio.load_mrc(self.__rec3d_file), 1, 2)
        # Flip Z-axis
        lio.write_mrc(np.flip(hold_rec, axis=2), self.__rec3d_file)

    def set_header(self, data="mics", p_size=None, origin=None):
        """Set pixel size and/or origin in an MRC header via IMOD alterheader.

        Args:
            data (str): Target file: 'mics' (default) for the
                micrograph stack, or 'rec3d' for the reconstruction.
            p_size (array-like, optional): Voxel size as (dx, dy,
                dz) in Angstroms.
            origin (array-like, optional): Origin as (ox, oy, oz)
                in Angstroms.

        Raises:
            ValueError: If data is invalid, or p_size/origin have
                != 3 elements.
            subprocess.CalledProcessError: If alterheader fails.
            IOError: If the log file cannot be written.
        """
        if data not in ("mics", "rec3d"):
            raise ValueError("data must be 'mics' or 'rec3d'.")
        if p_size is not None:
            if not hasattr(p_size, "__len__") or len(p_size) != 3:
                raise ValueError("p_size must have exactly 3 elements.")
        if origin is not None:
            if not hasattr(origin, "__len__") or len(origin) != 3:
                raise ValueError("origin must have exactly 3 elements.")

        aheader_cmd = [IMOD_CMD_AHEADER]
        if data == "mics":
            aheader_cmd += [
                str(self.__micgraphs_file),
            ]
        else:
            aheader_cmd += [
                str(self.__rec3d_file),
            ]
        if p_size is not None:
            aheader_cmd += [
                "-del",
                str(p_size[0]) + "," + str(p_size[1]) + "," + str(p_size[2]),
            ]
        if origin is not None:
            aheader_cmd += [
                "-org",
                str(origin[0]) + "," + str(origin[1]) + "," + str(origin[2]),
            ]

        try:
            with open(self.__log_file, "a", encoding="utf-8") as file_log:
                file_log.write(
                    "\n["
                    + time.strftime("%c")
                    + "]RUNNING COMMAND:-> "
                    + " ".join(aheader_cmd)
                    + "\n"
                )
                subprocess.check_call(
                    aheader_cmd, stdout=file_log, stderr=file_log
                )
        except subprocess.CalledProcessError:
            logger.error("Error calling the command: %s", aheader_cmd)
            raise
        except IOError:
            logger.error("Log file could not be written: %s", self.__log_file)
            raise

    def invert_mics_den(self):
        """Negate the density values of all micrographs in the stack."""
        mics = lio.load_mrc(self.__micgraphs_file)
        lio.write_mrc(-1 * mics, self.__micgraphs_file)

    def add_mics_misalignment(self, mn, mx, n_sigma=0):
        """Add sinusoidal X/Y misalignment to each micrograph.

        The misalignment magnitude follows the model:
        ``f(θ) = mn + mx * sin(θ) / sin(θ_max)``
        with optional Gaussian noise on top.

        Args:
            mn (float): Misalignment at zero tilt.
            mx (float): Misalignment at maximum tilt;
                must be >= mn.
            n_sigma (float): Standard deviation of additional
                Gaussian noise on the shifts (default 0).
        """

        if mx < mn:
            raise ValueError("mx must be >= mn.")

        rng = np.random.default_rng()

        mics = lio.load_mrc(self.__micgraphs_file)
        angs = np.abs(np.radians(self.__load_tangs_file()))
        n_angs = len(angs)
        shifts = (
            mn
            + np.sin(angs) / np.sin(angs.max())
            + rng.normal(0, n_sigma, n_angs)
        )
        split_fs = rng.uniform(0, 1, n_angs)
        for i, shift, split_f in zip(range(n_angs), shifts, split_fs):
            shift_x = shift / math.sqrt(split_f + 1)
            shift_y = split_f * shift_x
            mics[:, :, i] = sp_shift(
                mics[:, :, i],
                (shift_x, shift_y),
                output=None,
                order=3,
                mode="constant",
                cval=0.0,
                prefilter=True,
            )

        lio.write_mrc(mics, self.__micgraphs_file)

    def _resolve_snr(self, detector_snr):
        """Resolve SNR value from config (scalar, range, or None).

        Args:
            detector_snr (float | list | tuple | None): SNR
                specification — a scalar, a two-element
                ``[lo, hi]`` range, or None to skip noise.

        Returns:
            float | None: Resolved SNR value, or None.
        """
        if detector_snr is None:
            return None
        if isinstance(detector_snr, (list, tuple)) and len(detector_snr) >= 2:
            return round(
                (detector_snr[1] - detector_snr[0]) * random.random()
                + detector_snr[0],
                2,
            )
        elif isinstance(detector_snr, (list, tuple)):
            return detector_snr[0]
        return detector_snr

    def simulate(self, density, tem_params, v_size):
        """Run the full TEM simulation pipeline.

        Args:
            density: 3D density volume.
            tem_params: Dict from .tem config file.
            v_size: Voxel size in Angstroms.

        Returns:
            snr (float or None): The SNR used for detector noise.
        """
        tilt_rg = tem_params["TILT_ANGS_RG"]
        tilt_step = tem_params["TILT_ANGS_STEP"]
        tilt_angs = np.arange(tilt_rg[0], tilt_rg[1], tilt_step)

        logger.debug("Generating tilt series.")
        self.gen_tilt_series_imod(density, tilt_angs, ax="Y")

        self.add_mics_misalignment(
            tem_params["MALIGN_MIN"],
            tem_params["MALIGN_MAX"],
            tem_params["MALIGN_SG"],
        )

        snr = self._resolve_snr(tem_params.get("DETECTOR_SNR"))
        if snr is not None:
            logger.debug("Adding detector noise with SNR=%s.", snr)
            self.add_detector_noise(snr)

        self.invert_mics_den()
        self.set_header(data="mics", p_size=(v_size, v_size, v_size))

        logger.debug("Reconstructing 3D tomogram.")
        self.recon3D_imod()
        self.set_header(
            data="rec3d", p_size=(v_size, v_size, v_size), origin=(0, 0, 0)
        )

        return snr
