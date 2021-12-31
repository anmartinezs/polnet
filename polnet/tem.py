"""
Models for Transmision Electron Microscopy
    * VERY IMPORTANT: this package makes usege of some IMOD binaries, so they must be installed in the system

"""

__author__ = 'Antonio Martinez-Sanchez'

import os
import time
import subprocess
import numpy as np
from polnet import lio

## IMOD commands
IMOD_CMD_XYZPROJ = 'xyzproj'
IMOD_CMD_TILT = 'tilt'


class TEM:
    """
    Class for modeling a Transmission Electron Microscope
    """

    def __init__(self, work_dir):
        """
        :param work_dir: path to working directory where intermediate files are stored
        """
        self.__work_dir = work_dir
        self.__log_file = self.__work_dir + '/TEM.log'
        self.__create_work_dir()
        self.__vol_file = self.__work_dir + '/in_vol.mrc'
        self.__micgraphs_file = self.__work_dir + '/out_micrographs.mrc'
        self.__tangs_file = self.__work_dir + '/out_tangs.tlt'
        self.__rec3d_file = self.__work_dir + '/out_rec3d.mrc'

    def __create_work_dir(self):
        """
        Create the working directory
        """
        if not os.path.exists(self.__work_dir):
            os.mkdir(self.__work_dir)

    def __save_tangs_file(self, angs, out_file):
        """
        Sores an tilt angles file according IMOD format
        :param angs: non-empty iterable with the tilt angles
        :param out_file: output path
        """
        with open(out_file, 'w') as file:
            for ang in angs:
                file.write(str(ang) + '\n')

    def gen_tilt_series_imod(self, vol, angs, ax='X', mode='real'):
        """
        Generates the 2D projection series from a 3D volume using 'xyzproj' IMOD binary
        :param vol: input 3D volume
        :param angs: non-empty iterable with the tilt angles
        :param ax: tilt axis, either 'X', 'Y' or 'Z' (default 'X')
        :param mode: mode of output file, valid: 'byte', 'int' or 'real' (default)
        """

        # Input parsing
        assert isinstance(vol, np.ndarray) and (len(vol.shape) == 3)
        assert hasattr(angs, '__len__') and (len(angs) > 0)
        assert (ax == 'X') or (ax == 'Y') or (ax == 'Z')
        assert (mode == 'byte') or (mode == 'int') or (mode == 'real')

        # Call to IMOD binary (xyzproj)

        # Building the command
        xyzproj_cmd = [IMOD_CMD_XYZPROJ]
        in_vol_path = self.__vol_file
        lio.write_mrc(vol, in_vol_path)
        xyzproj_cmd += ['-inp', in_vol_path]
        out_vol_path = self.__micgraphs_file
        xyzproj_cmd += ['-o', out_vol_path]
        xyzproj_cmd += ['-ax', ax]
        xyzproj_cmd += ['-ta']
        tangles = str(angs[0])
        for i in range(1, len(angs)):
            tangles += ',' + str(angs[i])
        xyzproj_cmd += [tangles]
        if mode == 'byte':
            xyzproj_cmd += ['-m', '0']
        elif mode == 'int':
            xyzproj_cmd += ['-m', '1']
        elif mode == 'real':
            xyzproj_cmd += ['-m', '2']

        # Command calling
        try:
            with open(self.__log_file, 'a') as file_log:
                file_log.write('\n[' + time.strftime("%c") + ']RUNNING COMMAND:-> ' + ' '.join(xyzproj_cmd) + '\n')
                subprocess.call(xyzproj_cmd, stdout=file_log, stderr=file_log)
            self.__save_tangs_file(angs, self.__tangs_file)
        except subprocess.CalledProcessError:
            print('ERROR: Error calling the command:', xyzproj_cmd)
            raise subprocess.CalledProcessError
        except IOError:
            print('ERROR: Log file could not be written:', self.__log_file)
            raise IOError

    def add_detector_noise(self, snr):
        """
        Add detector noise to micrographs. Readout noise has Gaussian distribution and dark current is typically
        one order magnitude lower, so the last one is neglected.
        :param snr: target snr to determine the level of noise to be added, linear scale so greater than zero
        """

        # Input parsing
        assert os.path.exists(self.__micgraphs_file)
        assert snr > 0

        # Load micrographs
        mics = lio.load_mrc(self.__micgraphs_file)

        # Adding readout noise with Gaussian distribution to every micrograph
        for i in range(mics.shape[2]):
            mic = mics[:, :, i]
            mask = mic > 0
            mn = mic[mask].mean()
            sg_fg = mn / snr
            mics[:, :, i] = mic + np.random.normal(mn, sg_fg, mic.shape)

        # Update micrographs file
        lio.write_mrc(mics, self.__micgraphs_file)

    def recon3D_imod(self):
        """
        Performs a 3D reconstruction from the tilted series micrograph using 'tilt' IMOD binary
        """

        # Call to IMOD binary (tilt)

        # Building the command
        tilt_cmd = [IMOD_CMD_TILT]
        vol = lio.load_mrc(self.__vol_file, mmap=True, no_saxes=False)
        tilt_cmd += ['-inp', self.__micgraphs_file]
        tilt_cmd += ['-output', self.__rec3d_file]
        tilt_cmd += ['-TILTFILE', self.__tangs_file]
        tilt_cmd += ['-THICKNESS', str(vol.shape[0])]

        # Command calling
        try:
            with open(self.__log_file, 'a') as file_log:
                file_log.write('\n[' + time.strftime("%c") + ']RUNNING COMMAND:-> ' + ' '.join(tilt_cmd) + '\n')
                subprocess.call(tilt_cmd, stdout=file_log, stderr=file_log)
        except subprocess.CalledProcessError:
            print('ERROR: Error calling the command:', tilt_cmd)
            raise subprocess.CalledProcessError
        except IOError:
            print('ERROR: Log file could not be written:', self.__log_file)
            raise IOError

        # Swap Y-Z axes from the output given by IMOD
        lio.write_mrc(np.swapaxes(lio.load_mrc(self.__rec3d_file), 1, 2), self.__rec3d_file)


    def invert_mics_den(self):
        """
        Invert micrographs densities
        """
        mics = lio.load_mrc(self.__micgraphs_file)
        lio.write_mrc(-1 * mics, self.__micgraphs_file)