import sys
import csv
import time
import random
from tqdm import tqdm

from polnet.utils import *
from polnet import lio
from polnet import tem
from polnet import poly as pp
from polnet.network import (
    NetSAWLC,
    NetSAWLCInter,
    NetHelixFiber,
    NetHelixFiberB,
)
from polnet.polymer import FiberUnitSDimer, MTUnit, MB_DOMAIN_FIELD_STR
from polnet.stomo import (
    MmerFile,
    MbFile,
    SynthTomo,
    SetTomos,
    HelixFile,
    MTFile,
    ActinFile,
    MmerMbFile,
)
from polnet.lrandom import (
    EllipGen,
    SphGen,
    TorGen,
    PGenHelixFiberB,
    PGenHelixFiber,
    SGenUniform,
    SGenProp,
    OccGen,
)
from polnet.membrane import SetMembranes


def generate_digital_sample(
    ntomos,
    voi_shape,
    out_dir,
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
):
    # Common tomogram settings
    root_path = os.path.realpath(os.path.join(os.getcwd(), "data"))

    if verbose:
        print("root_path", root_path)

    # OUTPUT FILES
    # out_dir = os.path.realpath(root_path + '/../data_generated/polnet_test') # '/out_all_tomos_9-10' # '/only_actin' # '/out_rotations'
    os.makedirs(out_dir, exist_ok=True)

    if verbose:
        print("out_dir", out_dir)

    sample_dir = out_dir + "/sample"

    os.makedirs(sample_dir, exist_ok=True)

    if verbose:
        print("sample_dir", sample_dir)

    # OUTPUT LABELS
    lbl_mb = 1
    lbl_ac = 2
    lbl_mt = 3
    lbl_cp = 4
    lbl_mp = 5
    # lbl_br = 6

    ##### Main procedure

    set_stomos = SetTomos()
    vx_um3 = (voi_vsize * 1e-4) ** 3

    # Preparing intermediate directories
    clean_dir(sample_dir)

    if verbose:
        print("dir clean")

    # Save labels table
    unit_lbl = 1
    header_lbl_tab = ["MODEL", "LABEL"]
    with open(out_dir + "/labels_table.csv", "w") as file_csv:
        writer_csv = csv.DictWriter(
            file_csv, fieldnames=header_lbl_tab, delimiter="\t"
        )
        writer_csv.writeheader()
        for i in range(len(membranes_list)):
            writer_csv.writerow(
                {
                    header_lbl_tab[0]: membranes_list[i],
                    header_lbl_tab[1]: unit_lbl,
                }
            )
            unit_lbl += 1
        for i in range(len(helix_list)):
            writer_csv.writerow(
                {header_lbl_tab[0]: helix_list[i], header_lbl_tab[1]: unit_lbl}
            )
            unit_lbl += 1
        for i in range(len(proteins_list)):
            writer_csv.writerow(
                {
                    header_lbl_tab[0]: proteins_list[i],
                    header_lbl_tab[1]: unit_lbl,
                }
            )
            unit_lbl += 1
        for i in range(len(mb_proteins_list)):
            writer_csv.writerow(
                {
                    header_lbl_tab[0]: mb_proteins_list[i],
                    header_lbl_tab[1]: unit_lbl,
                }
            )
            unit_lbl += 1

    if verbose:
        tomos_ids = tqdm(range(ntomos), desc="Generating tomograms")
    else:
        tomos_ids = range(ntomos)

    # Loop for tomograms
    if verbose:
        print("Tomo loop!")
    for tomod_id in tomos_ids:
        if verbose:
            print("Generating tomogram number:", tomod_id)
        hold_time = time.time()

        # Generate the VOI and tomogram density
        if isinstance(voi_shape, str):
            voi = lio.load_mrc(voi_shape) > 0
            voi_off = np.zeros(shape=voi.shape, dtype=bool)
            voi_off[
                voi_offs[0][0] : voi_offs[0][1],
                voi_offs[1][0] : voi_offs[1][1],
                voi_offs[2][0] : voi_offs[2][1],
            ] = True
            voi = np.logical_and(voi, voi_off)
            del voi_off
        else:
            voi = np.zeros(shape=voi_shape, dtype=bool)
            voi[
                voi_offs[0][0] : voi_offs[0][1],
                voi_offs[1][0] : voi_offs[1][1],
                voi_offs[2][0] : voi_offs[2][1],
            ] = True
            voi_inital_invert = np.invert(voi)
        bg_voi = voi.copy()
        voi_voxels = voi.sum()
        tomo_lbls = np.zeros(shape=voi_shape, dtype=np.float32)
        tomo_den = np.zeros(shape=voi.shape, dtype=np.float32)
        synth_tomo = SynthTomo()
        poly_vtp, mbs_vtp, skel_vtp = None, None, None
        entity_id = 1
        mb_voxels, ac_voxels, mt_voxels, cp_voxels, mp_voxels = 0, 0, 0, 0, 0
        set_mbs = None

        # Membranes loop
        count_mbs, hold_den = 0, None
        for p_id, p_file in enumerate(membranes_list):
            if verbose:
                print("\tProcessing file:", p_file)

            # Loading the membrane file
            memb = MbFile()
            memb.load_mb_file(p_file)

            # Generating the occupancy
            hold_occ = memb.get_occ()
            if hasattr(hold_occ, "__len__"):
                hold_occ = OccGen(hold_occ).gen_occupancy()

            # Membrane random generation by type
            param_rg = (
                memb.get_min_rad(),
                math.sqrt(3) * max(voi_shape) * voi_vsize,
                memb.get_max_ecc(),
            )

            if memb.get_type() == "sphere":
                mb_sph_generator = SphGen(radius_rg=(param_rg[0], param_rg[1]))
                set_mbs = SetMembranes(
                    voi,
                    voi_vsize,
                    mb_sph_generator,
                    param_rg,
                    memb.get_thick_rg(),
                    memb.get_layer_s_rg(),
                    hold_occ,
                    memb.get_over_tol(),
                    bg_voi=bg_voi,
                )
                set_mbs.build_set(verbosity=verbose)
                hold_den = set_mbs.get_tomo()
                if memb.get_den_cf_rg() is not None:
                    hold_den *= mb_sph_generator.gen_den_cf(
                        memb.get_den_cf_rg()[0], memb.get_den_cf_rg()[1]
                    )
            elif memb.get_type() == "ellipse":
                mb_ellip_generator = EllipGen(
                    radius_rg=param_rg[:2], max_ecc=param_rg[2]
                )
                set_mbs = SetMembranes(
                    voi,
                    voi_vsize,
                    mb_ellip_generator,
                    param_rg,
                    memb.get_thick_rg(),
                    memb.get_layer_s_rg(),
                    hold_occ,
                    memb.get_over_tol(),
                    bg_voi=bg_voi,
                )
                set_mbs.build_set(verbosity=verbose)
                hold_den = set_mbs.get_tomo()
                if memb.get_den_cf_rg() is not None:
                    hold_den *= mb_ellip_generator.gen_den_cf(
                        memb.get_den_cf_rg()[0], memb.get_den_cf_rg()[1]
                    )
            elif memb.get_type() == "toroid":
                mb_tor_generator = TorGen(radius_rg=(param_rg[0], param_rg[1]))
                set_mbs = SetMembranes(
                    voi,
                    voi_vsize,
                    mb_tor_generator,
                    param_rg,
                    memb.get_thick_rg(),
                    memb.get_layer_s_rg(),
                    hold_occ,
                    memb.get_over_tol(),
                    bg_voi=bg_voi,
                )
                set_mbs.build_set(verbosity=verbose)
                hold_den = set_mbs.get_tomo()
                if memb.get_den_cf_rg() is not None:
                    hold_den *= mb_tor_generator.gen_den_cf(
                        memb.get_den_cf_rg()[0], memb.get_den_cf_rg()[1]
                    )
            else:
                if verbose:
                    print(
                        "ERROR: Membrane type",
                        memb.get_type(),
                        "not recognized!",
                    )
                sys.exit()

            # Density tomogram updating
            voi = set_mbs.get_voi()
            mb_mask = set_mbs.get_tomo() > 0
            mb_mask[voi_inital_invert] = False
            tomo_lbls[mb_mask] = entity_id
            count_mbs += set_mbs.get_num_mbs()
            mb_voxels += (tomo_lbls == entity_id).sum()
            tomo_den = np.maximum(tomo_den, hold_den)
            hold_vtp = set_mbs.get_vtp()
            pp.add_label_to_poly(hold_vtp, entity_id, "Entity", mode="both")
            pp.add_label_to_poly(hold_vtp, lbl_mb, "Type", mode="both")
            if poly_vtp is None:
                poly_vtp = hold_vtp
                skel_vtp = hold_vtp
            else:
                poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
                skel_vtp = pp.merge_polys(skel_vtp, hold_vtp)
            synth_tomo.add_set_mbs(
                set_mbs, "Membrane", entity_id, memb.get_type()
            )
            entity_id += 1

        # Get membranes poly
        if set_mbs is not None:
            mbs_vtp = vtk.vtkPolyData()
            mbs_vtp.DeepCopy(poly_vtp)

        # Loop for Helicoidal structures
        count_actins, count_mts = 0, 0
        for p_id, p_file in enumerate(helix_list):
            if verbose:
                print("\tProcessing file:", p_file)

            # Loading the helix file
            helix = HelixFile()
            helix.load_hx_file(p_file)

            # print("helix", helix.get_type())

            # Generating the occupancy
            hold_occ = helix.get_occ()
            if hasattr(hold_occ, "__len__"):
                hold_occ = OccGen(hold_occ).gen_occupancy()

            # print("hold_occ", hold_occ)

            # Helicoida random generation by type
            if helix.get_type() == "mt":
                helix = MTFile()
                helix.load_mt_file(p_file)
                # Fiber unit generation
                funit = MTUnit(
                    helix.get_mmer_rad(),
                    helix.get_rad(),
                    helix.get_nunits(),
                    voi_vsize,
                )
                model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
                # Helix Fiber parameters model
                pol_generator = PGenHelixFiber()
                # Network generation
                net_helix = NetHelixFiber(
                    voi,
                    voi_vsize,
                    helix.get_l() * helix.get_mmer_rad() * 2,
                    model_surf,
                    pol_generator,
                    hold_occ,
                    helix.get_min_p_len(),
                    helix.get_hp_len(),
                    helix.get_mz_len(),
                    helix.get_mz_len_f(),
                    helix.get_over_tol(),
                    (helix.get_rad() + 0.5 * helix.get_mmer_rad()) * 2.4,
                )
                if helix.get_min_nmmer() is not None:
                    net_helix.set_min_nmmer(helix.get_min_nmmer())
                net_helix.build_network()
            elif helix.get_type() == "actin":
                # print("go actin")
                helix = ActinFile()
                helix.load_ac_file(p_file)
                # print("loaded")
                # Fiber unit generation
                funit = FiberUnitSDimer(helix.get_mmer_rad(), voi_vsize)
                model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
                # print("generated")
                # Helix Fiber parameters model
                pol_generator = PGenHelixFiberB()
                # Network generation
                net_helix = NetHelixFiberB(
                    voi,
                    voi_vsize,
                    helix.get_l() * helix.get_mmer_rad() * 2,
                    model_surf,
                    pol_generator,
                    hold_occ,
                    helix.get_min_p_len(),
                    helix.get_hp_len(),
                    helix.get_mz_len(),
                    helix.get_mz_len_f(),
                    helix.get_bprop(),
                    helix.get_p_branch(),
                    helix.get_over_tol(),
                )

                if helix.get_min_nmmer() is not None:
                    net_helix.set_min_nmmer(helix.get_min_nmmer())
                net_helix.build_network()

                # Geting branches poly
                br_vtp = pp.points_to_poly_spheres(
                    points=[
                        [0, 0, 0],
                    ],
                    rad=helix.get_mmer_rad(),
                )

                lio.save_vtp(
                    net_helix.get_branches_vtp(shape_vtp=br_vtp),
                    sample_dir + "/poly_br_" + str(tomod_id) + ".vtp",
                )
            else:
                if verbose:
                    print(
                        "ERROR: Helicoidal type",
                        helix.get_type(),
                        "not recognized!",
                    )
                sys.exit()

            # Density tomogram updating
            # voi = net_helix.get_voi()
            # tomo_den = np.maximum(tomo_den, net_helix.get_tomo())
            model_mask = model_svol < 0.05
            # off = .5 * np.asarray(model_svol.shape) - center
            net_helix.insert_density_svol(
                model_mask, voi, voi_vsize, merge="min", off_svol=None
            )
            if helix.get_den_cf_rg() is None:
                cte_val = 1
            else:
                cte_val = pol_generator.gen_den_cf(
                    helix.get_den_cf_rg()[0], helix.get_den_cf_rg()[1]
                )
            net_helix.insert_density_svol(
                model_svol * cte_val, tomo_den, voi_vsize, merge="max"
            )
            hold_lbls = np.zeros(shape=tomo_lbls.shape, dtype=np.float32)
            net_helix.insert_density_svol(
                np.invert(model_mask), hold_lbls, voi_vsize, merge="max"
            )
            tomo_lbls[hold_lbls > 0] = entity_id
            # lio.write_mrc(hold_lbls.astype(np.float32), '/fs/pool/pool-lucic2/antonio/polnet/riboprot/hold.mrc')
            hold_vtp = net_helix.get_vtp()
            hold_skel_vtp = net_helix.get_skel()
            pp.add_label_to_poly(hold_vtp, entity_id, "Entity", mode="both")
            pp.add_label_to_poly(
                hold_skel_vtp, entity_id, "Entity", mode="both"
            )
            if helix.get_type() == "mt":
                pp.add_label_to_poly(hold_vtp, lbl_mt, "Type", mode="both")
                pp.add_label_to_poly(
                    hold_skel_vtp, lbl_mt, "Type", mode="both"
                )
                count_mts += net_helix.get_num_pmers()
                mt_voxels += (tomo_lbls == entity_id).sum()
            elif helix.get_type() == "actin":
                pp.add_label_to_poly(hold_vtp, lbl_ac, "Type", mode="both")
                pp.add_label_to_poly(
                    hold_skel_vtp, lbl_ac, "Type", mode="both"
                )
                count_actins += net_helix.get_num_pmers()
                ac_voxels += (tomo_lbls == entity_id).sum()
            if poly_vtp is None:
                poly_vtp = hold_vtp
                skel_vtp = hold_skel_vtp
            else:
                poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
                skel_vtp = pp.merge_polys(skel_vtp, hold_skel_vtp)
            synth_tomo.add_network(
                net_helix, "Helix", entity_id, code=helix.get_type()
            )
            entity_id += 1

        # Loop for the list of input proteins loop
        count_prots = 0
        model_surfs, models, model_masks, model_codes = (
            list(),
            list(),
            list(),
            list(),
        )
        for p_id, p_file in enumerate(proteins_list):
            if verbose:
                print("\tProcessing file:", p_file)

            # Loading the protein
            protein = MmerFile(p_file)

            # Generating the occupancy
            hold_occ = protein.get_pmer_occ()
            if hasattr(hold_occ, "__len__"):
                hold_occ = OccGen(hold_occ).gen_occupancy()

            # Genrate the SAWLC network associated to the input protein
            # Polymer parameters
            # To read macromolecular models first we try to find the absolute path and secondly the relative to root_path
            try:
                model = lio.load_mrc(protein.get_mmer_svol())
            except FileNotFoundError:
                model = lio.load_mrc(root_path + protein.get_mmer_svol())

            # print("model", model.shape)

            # model = lio.load_mrc(root_path + '/' + protein.get_mmer_svol())
            model = lin_map(model, lb=0, ub=1)
            model = vol_cube(model)
            model_mask = model < protein.get_iso()
            model[model_mask] = 0
            model_surf = pp.iso_surface(
                model, protein.get_iso(), closed=False, normals=None
            )
            if surf_dec is not None:
                model_surf = pp.poly_decimate(model_surf, surf_dec)

            # print("model", model.shape)

            center = 0.5 * np.asarray(model.shape, dtype=float)
            # Monomer centering
            model_surf = pp.poly_translate(model_surf, -center)
            # Voxel resolution scaling
            model_surf = pp.poly_scale(model_surf, voi_vsize)
            model_surfs.append(model_surf)
            surf_diam = pp.poly_diam(model_surf) * protein.get_pmer_l()
            models.append(model)
            # print("model_mask", model_mask.shape, model_mask.sum())
            model_masks.append(model_mask)
            model_codes.append(protein.get_mmer_id())

            # Network generation
            # print("surfdiam", surf_diam)
            pol_l_generator = PGenHelixFiber()

            net_sawlc = NetSAWLC(
                voi,
                voi_vsize,
                protein.get_pmer_l() * surf_diam,
                model_surf,
                protein.get_pmer_l_max(),
                pol_l_generator,
                hold_occ,
                protein.get_pmer_over_tol(),
                poly=None,
                svol=model < protein.get_iso(),
                tries_mmer=mmer_tries,
                tries_pmer=pmer_tries,
            )
            # net_sawlc = NetSAWLCInter(voi, voi_vsize, surf_diams, model_surfs, protein.get_pmer_l_max(),
            #                           pol_l_generator, pol_s_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(),
            #                           poly=None, svols=model_masks, codes=model_codes, compaq=5.5)
            net_sawlc.build_network()

            # Density tomogram updating
            net_sawlc.insert_density_svol(
                model_mask, voi, voi_vsize, merge="min"
            )
            net_sawlc.insert_density_svol(
                model, tomo_den, voi_vsize, merge="max"
            )
            hold_lbls = np.zeros(shape=tomo_lbls.shape, dtype=np.float32)
            net_sawlc.insert_density_svol(
                np.invert(model_mask), hold_lbls, voi_vsize, merge="max"
            )
            tomo_lbls[hold_lbls > 0] = entity_id
            count_prots += net_sawlc.get_num_mmers()
            if verbose:
                print(f"{net_sawlc.get_num_mmers()} proteins added.")
            cp_voxels += (tomo_lbls == entity_id).sum()
            hold_vtp = net_sawlc.get_vtp()
            hold_skel_vtp = net_sawlc.get_skel()
            pp.add_label_to_poly(hold_vtp, entity_id, "Entity", mode="both")
            pp.add_label_to_poly(
                hold_skel_vtp, entity_id, "Entity", mode="both"
            )
            pp.add_label_to_poly(hold_vtp, lbl_cp, "Type", mode="both")
            pp.add_label_to_poly(hold_skel_vtp, lbl_cp, "Type", mode="both")
            if poly_vtp is None:
                poly_vtp = hold_vtp
                skel_vtp = hold_skel_vtp
            else:
                poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
                skel_vtp = pp.merge_polys(skel_vtp, hold_skel_vtp)
            synth_tomo.add_network(
                net_sawlc, "SAWLC", entity_id, code=protein.get_mmer_id()
            )
            entity_id += 1

        # Loop for the list of input proteins loop
        count_mb_prots = 0
        if mbs_vtp is None:
            if len(mb_proteins_list) > 0:
                if verbose:
                    print(
                        "WARNING: membrane proteins can not inserted because there is no membrane surfaces!"
                    )
        else:
            model_surfs, surf_diams, models, model_masks, model_codes = (
                list(),
                list(),
                list(),
                list(),
                list(),
            )
            for p_id, p_file in enumerate(mb_proteins_list):
                if verbose:
                    print("\tProcessing file:", p_file)

                # Loading the membrane protein
                protein = MmerMbFile(p_file)

                # Generating the occupancy
                hold_occ = protein.get_pmer_occ()
                if hasattr(hold_occ, "__len__"):
                    hold_occ = OccGen(hold_occ).gen_occupancy()

                # Insert membrane bound densities in a Polymer
                # Polymer parameters
                # To read macromolecular models first we try to find the absolute path and secondly the relative to root_path
                try:
                    model = lio.load_mrc(protein.get_mmer_svol())
                except FileNotFoundError:
                    model = lio.load_mrc(root_path + protein.get_mmer_svol())
                model = lin_map(model, lb=0, ub=1)
                model_mask = model < protein.get_iso()
                model[model_mask] = 0
                model_surf = iso_surface(
                    model, protein.get_iso(), closed=False, normals=None
                )
                center = (
                    protein.get_mmer_center()
                )  # .5 * np.asarray(model.shape, dtype=float)
                if center is None:
                    center = 0.5 * (np.asarray(model.shape, dtype=float) - 1)
                    off = np.asarray((0.0, 0.0, 0.0))
                else:
                    center = np.asarray(center)
                    off = 0.5 * np.asarray(model.shape) - center
                # Adding membrane domain to monomer surface
                mb_domain_mask = np.ones(shape=model.shape, dtype=bool)
                hold_mb_z_height = protein.get_mb_z_height()
                if hold_mb_z_height is None:
                    hold_mb_z_height = int(round(center[2] + 2.5 / voi_vsize))
                for z in range(hold_mb_z_height + 1, model.shape[2]):
                    mb_domain_mask[:, :, z] = 0
                pp.add_sfield_to_poly(
                    model_surf,
                    mb_domain_mask,
                    MB_DOMAIN_FIELD_STR,
                    dtype="float",
                    interp="NN",
                    mode="points",
                )
                # Monomer centering
                model_surf = pp.poly_translate(model_surf, -center)
                # Voxel resolution scaling
                model_surf = pp.poly_scale(model_surf, voi_vsize)
                surf_diam = pp.poly_diam(model_surf)
                pol_l_generator = PGenHelixFiber()
                # Network generation
                if protein.get_pmer_reverse_normals():
                    mbs_vtp = pp.poly_reverse_normals(mbs_vtp)
                net_sawlc = NetSAWLC(
                    voi,
                    voi_vsize,
                    protein.get_pmer_l() * surf_diam,
                    model_surf,
                    protein.get_pmer_l_max(),
                    pol_l_generator,
                    hold_occ,
                    protein.get_pmer_over_tol(),
                    poly=mbs_vtp,
                    svol=model < protein.get_iso(),
                    tries_mmer=mmer_tries,
                    tries_pmer=pmer_tries,
                )
                # net_sawlc = NetSAWLCInter(voi, voi_vsize, protein.get_pmer_l() * surf_diam, model_surf, protein.get_pmer_l_max(),
                #                      pol_l_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(), poly=mb_poly,
                #                      svol=model < protein.get_iso())
                net_sawlc.build_network()
                # voi = net_sawlc.get_voi()

                # lio.write_mrc(voi.astype(np.float32), root_path + '/hold_voi.mrc')
                # lio.write_mrc(set_mbs.get_tomo().astype(np.float32), root_path + '/hold_den.mrc')

                # Density tomogram updating
                net_sawlc.insert_density_svol(
                    model_mask, voi, voi_vsize, merge="min"
                )
                net_sawlc.insert_density_svol(
                    model, tomo_den, voi_vsize, merge="max"
                )
                hold_lbls = np.zeros(shape=tomo_lbls.shape, dtype=np.float32)
                net_sawlc.insert_density_svol(
                    np.invert(model_mask), hold_lbls, voi_vsize, merge="max"
                )
                tomo_lbls[hold_lbls > 0] = entity_id
                count_mb_prots += net_sawlc.get_num_mmers()
                if verbose:
                    print(
                        f"{net_sawlc.get_num_mmers()} membrane bound proteins added."
                    )
                mp_voxels += (tomo_lbls == entity_id).sum()
                hold_vtp = net_sawlc.get_vtp()
                hold_skel_vtp = net_sawlc.get_skel()
                pp.add_label_to_poly(
                    hold_vtp, entity_id, "Entity", mode="both"
                )
                pp.add_label_to_poly(
                    hold_skel_vtp, entity_id, "Entity", mode="both"
                )
                pp.add_label_to_poly(hold_vtp, lbl_mp, "Type", mode="both")
                pp.add_label_to_poly(
                    hold_skel_vtp, lbl_mp, "Type", mode="both"
                )
                if poly_vtp is None:
                    poly_vtp = hold_vtp
                    skel_vtp = hold_skel_vtp
                else:
                    poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
                    skel_vtp = pp.merge_polys(skel_vtp, hold_skel_vtp)
                synth_tomo.add_network(
                    net_sawlc,
                    "Mb-SAWLC",
                    entity_id,
                    code=protein.get_mmer_id(),
                )
                entity_id += 1

        # Tomogram statistics
        if verbose:
            print("\t\t-Tomogram", str(tomod_id), "Density Statistics:")
            print(
                "\t\t\t+Membranes:",
                count_mbs,
                "#, ",
                mb_voxels * vx_um3,
                "um**3, ",
                100.0 * (mb_voxels / voi_voxels),
                "%",
            )
            print(
                "\t\t\t+Actin:",
                count_actins,
                "#, ",
                ac_voxels * vx_um3,
                "um**3, ",
                100.0 * (ac_voxels / voi_voxels),
                "%",
            )
            print(
                "\t\t\t+Microtublues:",
                count_mts,
                "#, ",
                mt_voxels * vx_um3,
                "um**3, ",
                100.0 * (mt_voxels / voi_voxels),
                "%",
            )
            print(
                "\t\t\t+Proteins:",
                count_prots,
                "#, ",
                cp_voxels * vx_um3,
                "um**3, ",
                100.0 * (cp_voxels / voi_voxels),
                "%",
            )
            print(
                "\t\t\t+Membrane proteins:",
                count_mb_prots,
                "#, ",
                mp_voxels * vx_um3,
                "um**3, ",
                100.0 * (mp_voxels / voi_voxels),
                "%",
            )
        counts_total = (
            count_mbs + count_actins + count_mts + count_prots + count_mb_prots
        )
        total_voxels = (
            mb_voxels + ac_voxels + mt_voxels + cp_voxels + mp_voxels
        )

        if verbose:
            print(
                "\t\t\t+Total:",
                counts_total,
                "#, ",
                total_voxels * vx_um3,
                "um**3, ",
                100.0 * (total_voxels / voi_voxels),
                "%",
            )
        if verbose:
            print(
                "\t\t\t+Time for generation: ",
                (time.time() - hold_time) / 60,
                "mins",
            )

        # Storing simulated density results
        tomo_den_out = sample_dir + "/tomo_den_" + str(tomod_id) + ".mrc"
        lio.write_mrc(tomo_den, tomo_den_out, v_size=voi_vsize)
        synth_tomo.set_den(tomo_den_out)
        tomo_lbls_out = sample_dir + "/tomo_lbls_" + str(tomod_id) + ".mrc"
        lio.write_mrc(tomo_lbls, tomo_lbls_out, v_size=voi_vsize)
        poly_den_out = sample_dir + "/poly_den_" + str(tomod_id) + ".vtp"
        lio.save_vtp(poly_vtp, poly_den_out)
        synth_tomo.set_poly(poly_den_out)
        poly_skel_out = sample_dir + "/poly_skel_" + str(tomod_id) + ".vtp"
        lio.save_vtp(skel_vtp, poly_skel_out)

        # We don't compute acquisition here
        synth_tomo.set_mics("no_acquisition.mrc")
        synth_tomo.set_tomo("no_acquisition.mrc")

        # Update the set
        set_stomos.add_tomos(synth_tomo)

    # Storing tomograms CSV file
    set_stomos.save_csv(out_dir + "/digital_sample_motif_list.csv")
    if verbose:
        print(out_dir + "/digital_sample_motif_list.csv")
        print("Successfully terminated. (" + time.strftime("%c") + ")")

    time.sleep(1)
