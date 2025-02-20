import sys
import os
import multiprocessing
from tqdm import tqdm
import time
import yaml

sys.path.append("/app")

from gui.core.digital_sample import generate_digital_sample
from gui.core.acquisition import generate_acquisition


# Load parameters from yaml config file
path = "/app/generation_config.yaml"

with open(path, "r") as stream:
    config = yaml.safe_load(stream)


# Common options
step = config["step"]
out_dir = os.path.join("/app/outdir", config["out_dir"])
root_path = os.path.join("/app", config["root_path"])
verbose = config["verbose"]
multi_process = config["multi_process"]
ntomos = config["ntomos"]
voi_vsize = config["voi_vsize"]


if step == "sample":
    nb_tomos_per_process_max = config["nb_tomos_per_process_max"]
    voi_shape = config["voi_shape"]
    mmer_tries = config["mmer_tries"]
    pmer_tries = config["pmer_tries"]
    membranes_list = config["membranes_list"]
    helix_list = config["helix_list"]
    proteins_list = config["proteins_list"]
    mb_proteins_list = config["mb_proteins_list"]
    surf_dec = config["surf_dec"]

    # Let's complete args
    # if list are none then list empty

    membranes_list = (
        [os.path.join(root_path, file) for file in membranes_list]
        if membranes_list
        else []
    )
    helix_list = (
        [os.path.join(root_path, file) for file in helix_list]
        if helix_list
        else []
    )
    proteins_list = (
        [os.path.join(root_path, file) for file in proteins_list]
        if proteins_list
        else []
    )
    mb_proteins_list = (
        [os.path.join(root_path, file) for file in mb_proteins_list]
        if mb_proteins_list
        else []
    )
    voi_offs = [
        [4, voi_shape[0] - 4],
        [4, voi_shape[1] - 4],
        [4, voi_shape[2] - 4],
    ]

    generate_func = generate_digital_sample


elif step == "acquisition":
    tilt_angs = config["tilt_angs"]
    detector_snr = config["detector_snr"]
    malign_mn = config["malign_mn"]
    malign_mx = config["malign_mx"]
    malign_sg = config["malign_sg"]
    acquisition_name = config["acquisition_name"]

    # Let's complete args
    tilt_angs = range(*tilt_angs)  # [-60, 60, 3]

    generate_func = generate_acquisition


start_time = time.time()


if multi_process:
    nb_cpu = multiprocessing.cpu_count()
    nb_processes_simult = nb_cpu // 2

    print(f"Number of CPUs: {nb_cpu},used: {nb_processes_simult}")

    if step == "sample":

        if ntomos <= nb_processes_simult:
            ntomos_per_process_list = [1 for _ in range(ntomos)]
        else:
            ntomos_per_process_round = min(
                ntomos // nb_processes_simult, nb_tomos_per_process_max
            )

            if ntomos_per_process_round < nb_tomos_per_process_max:
                ntomos_per_process_list = [
                    ntomos_per_process_round
                ] * nb_processes_simult

                for i in range(ntomos % nb_processes_simult):
                    ntomos_per_process_list[i] += 1

            else:
                ntomos_per_process_list = [ntomos_per_process_round] * (
                    ntomos // nb_tomos_per_process_max
                )

                if ntomos % nb_tomos_per_process_max != 0:
                    ntomos_per_process_list.append(
                        ntomos % nb_tomos_per_process_max
                    )

        nb_process_total = len(ntomos_per_process_list)
        nb_processes_simult = min(nb_processes_simult, nb_process_total)
        ntomos_per_process_avg = ntomos / nb_process_total

        print(
            f"Number of tomograms: {ntomos,sum(ntomos_per_process_list)} (both should be equal)"
        )
        print(
            f"Average Number of tomograms per process: {ntomos_per_process_avg}"
        )
        print(f"Number of processes: {nb_process_total}")
        print(f"Number of simultaneous processes: {nb_processes_simult}")

        args = [
            (
                cur_ntomos_per_process,
                voi_shape,
                os.path.join(
                    out_dir,
                    f"generation_{id_process:0{len(str(nb_process_total))}}",
                ),
                voi_offs,
                voi_vsize,
                mmer_tries,
                pmer_tries,
                membranes_list,
                helix_list,
                proteins_list,
                mb_proteins_list,
                surf_dec,
                verbose,
            )
            for id_process, cur_ntomos_per_process in zip(
                range(nb_process_total), ntomos_per_process_list
            )
        ]
    else:
        # Find all directories in out_dir that start with "generation_"
        dirs = [
            os.path.join(out_dir, d)
            for d in os.listdir(out_dir)
            if os.path.isdir(os.path.join(out_dir, d))
            and d.startswith("generation_")
        ]

        nb_process_total = len(dirs)
        nb_processes_simult = min(nb_processes_simult, nb_process_total)
        ntomos_per_process_avg = ntomos / nb_process_total

        args = [
            (
                dir_path,
                voi_vsize,
                tilt_angs,
                detector_snr,
                malign_mn,
                malign_mx,
                malign_sg,
                acquisition_name,
                verbose,
            )
            for dir_path in dirs
        ]

    with multiprocessing.Pool(processes=nb_processes_simult) as pool:
        result = pool.starmap_async(generate_func, args)
        with tqdm(
            total=ntomos, desc=f"Tomos completed (estimation): "
        ) as pbar:
            while not result.ready():
                time.sleep(
                    1
                )  # Sleep for a short duration to avoid excessive CPU usage
                remaining_tasks = (
                    result._number_left
                )  # Get the number of remaining tasks
                pbar.n = (
                    nb_process_total - remaining_tasks
                ) * ntomos_per_process_avg
                # Set the progress bar to the correct value
                pbar.refresh()  # Refresh the progress bar

        result.get()

end_time = time.time()
print(
    f"{step} done in {end_time - start_time:.2f} seconds, {(end_time - start_time)/ntomos:.2f} seconds per tomo."
)

# Save yaml file in outdir
if step == "sample":
    out_config_path = os.path.join(out_dir, "sample_config.yaml")
elif step == "acquisition":
    out_config_path = os.path.join(
        out_dir, f"acquisition_{acquisition_name}_config.yaml"
    )

with open(out_config_path, "w") as f:
    yaml.dump(config, f)
