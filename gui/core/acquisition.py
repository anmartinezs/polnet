import time
import random
import os

import pandas as pd
import shutil

from polnet import lio
from polnet import tem
from polnet.utils import clean_dir


def generate_acquisition(
    out_dir,
    voi_vsize,
    tilt_angs,
    detector_snr,
    malign_mn,
    malign_mx,
    malign_sg,
    acquisition_name,
    verbose,
):
    """
    Generate a set of synthetic tomograms from a list of motifs
    """

    tem_dir = out_dir + f"/tem_{acquisition_name}"
    tomos_dir = out_dir + f"/tomos_{acquisition_name}"

    os.makedirs(tomos_dir, exist_ok=True)
    os.makedirs(tem_dir, exist_ok=True)

    if verbose:
        print("tem_dir", tem_dir)
        print("tomos_dir", tomos_dir)

    clean_dir(tem_dir)
    clean_dir(tomos_dir)

    # Load motif
    input_motif_path = out_dir + "/digital_sample_motif_list.csv"
    motif_csv = pd.read_csv(input_motif_path, delimiter="\t")

    density_names = motif_csv["Density"].unique()

    for density_name in density_names:
        tomod_id = density_name.split("_")[-1].replace(".mrc", "")

        tomo_den_out = density_name

        # TEM for 3D reconstructions
        if verbose:
            print("start TEM")
        start_tem_time = time.time()
        temic = tem.TEM(tem_dir)
        vol = lio.load_mrc(tomo_den_out)
        start_time_tilt = time.time()
        temic.gen_tilt_series_imod(vol, tilt_angs, ax="Y")
        temic.add_mics_misalignment(malign_mn, malign_mx, malign_sg)
        end_time_tilt = time.time()
        if verbose:
            print(
                f"tilt series took {end_time_tilt - start_time_tilt} seconds"
            )
        if detector_snr is not None:
            if hasattr(detector_snr, "__len__"):
                if len(detector_snr) >= 2:
                    snr = round(
                        (detector_snr[1] - detector_snr[0]) * random.random()
                        + detector_snr[0],
                        2,
                    )
                else:
                    snr = detector_snr[0]
            else:
                snr = detector_snr
            temic.add_detector_noise(snr)
        temic.invert_mics_den()
        temic.set_header(data="mics", p_size=(voi_vsize, voi_vsize, voi_vsize))
        temic.recon3D_imod()
        temic.set_header(
            data="rec3d",
            p_size=(voi_vsize, voi_vsize, voi_vsize),
            origin=(0, 0, 0),
        )
        if detector_snr is not None:
            out_mics, out_tomo_rec = (
                tomos_dir + "/tomo_mics_" + str(tomod_id)
                # + "_snr"
                # + str(snr)
                + ".mrc",
                tomos_dir + "/tomo_rec_" + str(tomod_id)
                # + "_snr"
                # + str(snr)
                + ".mrc",
            )
        else:
            out_mics, out_tomo_rec = (
                tomos_dir + "/tomo_mics_" + str(tomod_id) + ".mrc",
                tomos_dir + "/tomo_rec_" + str(tomod_id) + ".mrc",
            )
        shutil.copyfile(tem_dir + "/out_micrographs.mrc", out_mics)
        shutil.copyfile(tem_dir + "/out_rec3d.mrc", out_tomo_rec)

        motif_csv.loc[
            motif_csv["Density"] == density_name,
            "Micrographs",
        ] = out_mics

        motif_csv.loc[
            motif_csv["Density"] == density_name,
            "Tomo3D",
        ] = out_tomo_rec

        end_tem_time = time.time()
        if verbose:
            print(
                f"end tomo gen, it took {end_tem_time - end_time_tilt} seconds"
            )
            print(f"end TEM, it took {end_tem_time - start_tem_time} seconds")

    # Save csv to output_path
    motif_csv.to_csv(
        out_dir + f"/tomos_{acquisition_name}_motif_list.csv",
        index=False,
        sep="\t",
    )

    # Remove dir with temp files
    clean_dir(tem_dir)
    os.rmdir(tem_dir)
