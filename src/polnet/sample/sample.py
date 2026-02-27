"""Synthetic sample representation for cryo-ET simulation.

Defines :class:`SyntheticSample`, which models the full 3-D volume
of interest (VOI) containing membranes, filaments, cytosolic
proteins, and membrane-bound proteins, together with their
ground-truth geometric annotations.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import random
from pathlib import Path

import numpy as np
import vtk

from .filaments import FlmsFactory
from .membranes import MbFactory
from .pms import (
    PmGen,
    PmSAWLCPolyNet,
)
from .pns import (
    PnGen,
    PnSAWLCNet,
)
from ..logging_conf import _LOGGER as logger
from ..utils.poly import (
    add_label_to_poly,
    merge_polys,
)

MOTIF_DTYPE = np.dtype(
    [
        ("type", "U12"),  # "membrane", "filament", "cprotein", "mbprotein"
        ("label", np.int32),  # entity_id
        ("code", "U32"),  # "mt", "actin", "pdb_1bxn", "sphere", etc.
        (
            "polymer",
            np.int32,
        ),  # polymer index (networks) or point index (membranes)
        ("x", np.float64),
        ("y", np.float64),
        ("z", np.float64),
        ("q1", np.float64),  # quaternion w (networks) or normal nx (membranes)
        ("q2", np.float64),
        ("q3", np.float64),
        ("q4", np.float64),  # quaternion z (networks) or 0 (membranes)
    ]
)


class SyntheticSample:
    """A model for a synthetic Cryo-ET sample."""

    OUTPUT_LABELS = {
        "membrane": 1,
        "actin": 2,
        "microtubule": 3,
        "cprotein": 4,
        "mbprotein": 5,
    }

    def __init__(self, shape: tuple, v_size: float, offset=(4, 4, 4)):
        """Constructor.

        Args:
            shape (tuple): The shape of the sample.
            v_size (float): The voxel size of the sample in angstroms.
            offset (tuple, optional): The offset for the sample. Defaults to (4, 4, 4).
        """
        self._shape = shape
        self._v_size = v_size
        self._offset = offset
        self.reset()

    def reset(self):
        """Reset the sample to its initial state. Call this method to clear all data."""
        self._voi = np.zeros(shape=self._shape, dtype=bool)
        self._voi[
            self._offset[0] : self._shape[0] - self._offset[0],
            self._offset[1] : self._shape[1] - self._offset[1],
            self._offset[2] : self._shape[2] - self._offset[2],
        ] = True
        self._voi_voxels = self._voi.sum()
        self._bg_voi = self._voi.copy()
        self._labels = np.zeros(
            shape=self._shape, dtype=np.uint8
        )  # Up to 255 unique labels, which should be sufficient for most samples. Adjust dtype if more are needed.
        self._density = np.zeros(shape=self._shape, dtype=np.float32)
        self._poly_vtp = None
        self._skel_vtp = None
        self._mbs_vtp = None
        self._structure_counts = {}
        self._voxel_counts = {}
        self._entity_id_counter = 1
        self._label_registry = []
        self._motifs = np.empty(0, dtype=MOTIF_DTYPE)

    @property
    def density(self) -> np.ndarray:
        """Deep copy of the current density tomogram."""
        return self._density.copy()

    @property
    def labels(self) -> np.ndarray:
        """Deep copy of the per-voxel entity-ID label volume."""
        return self._labels.copy()

    @property
    def poly_vtp(self) -> vtk.vtkPolyData | None:
        """Merged polygon surface of all inserted components, or None."""
        return self._poly_vtp

    @property
    def skel_vtp(self) -> vtk.vtkPolyData | None:
        """Merged skeleton polydata of all inserted components, or None."""
        return self._skel_vtp

    @property
    def mbs_vtp(self) -> vtk.vtkPolyData | None:
        """Merged membrane-only polygon surface, or None."""
        return self._mbs_vtp

    @property
    def v_size(self) -> float:
        """Voxel size in Angstroms."""
        return self._v_size

    @property
    def motifs(self) -> np.ndarray:
        """Per-monomer ground truth as structured numpy array.

        Fields:
            type (membrane, actin, microtubule, cprotein, mbprotein)
            label (entity_id)
            code (model identifier)
            polymer (polymer index for networks or point index for membranes)
            x/y/z (coordinates of motif center)
            q1/q2/q3/q4 (quaternion w/x/y/z for networks or normal nx/ny/nz/0 for membranes).
        """
        return self._motifs

    @property
    def label_registry(self) -> list:
        """List of (model_name, entity_id) tuples for all components in the sample."""
        return self._label_registry.copy()

    def structure_count(self, struct_type: str) -> int:
        """Get the structure count for a given type.

        Args:
            struct_type (str): The type of structure to retrieve (e.g. 'membrane', 'actin').

        Returns:
            int: The structure count for the specified type.
        """
        if struct_type not in self._structure_counts:
            logger.warning("%s structure count not found.", struct_type)
        return self._structure_counts.get(struct_type, 0)

    def _commit_component(
        self,
        type_key: str,
        model_name: str,
        count: int,
        label_mask: np.ndarray,
        hold_vtp: vtk.vtkPolyData,
        hold_skel_vtp=None,
        is_membrane: bool = False,
    ) -> None:
        """Shared bookkeeping after inserting any component.

        Args:
            type_key (str): The key for the component type (e.g. 'actin', 'microtubule', 'membrane').
            model_name (str): The model identifier for the label registry (e.g. 'sphere', 'mt', 'pdb_1bxn').
            count (int): The count of the component (e.g. number of filaments, number of membranes).
            label_mask (np.ndarray): A boolean mask indicating where to place the labels for this component in the sample.
            hold_vtp: The polydata for the component.
            hold_skel_vtp: The skeletal polydata for the component, if applicable.
            is_membrane (bool, optional): Whether the component is a membrane (used for special handling of membrane VTP). Defaults to False.
        """
        self._label_registry.append((model_name, self._entity_id_counter))

        self._labels[label_mask > 0] = self._entity_id_counter

        if type_key not in self._structure_counts:
            self._structure_counts[type_key] = 0
            self._voxel_counts[type_key] = 0
        self._structure_counts[type_key] += count
        self._voxel_counts[type_key] += (
            self._labels == self._entity_id_counter
        ).sum()

        type_label = self.OUTPUT_LABELS.get(type_key, self._entity_id_counter)
        add_label_to_poly(
            hold_vtp, self._entity_id_counter, "Entity", mode="both"
        )
        add_label_to_poly(hold_vtp, type_label, "Type", mode="both")
        if hold_skel_vtp is None:
            hold_skel_vtp = hold_vtp
        else:
            add_label_to_poly(
                hold_skel_vtp, self._entity_id_counter, "Entity", mode="both"
            )
            add_label_to_poly(hold_skel_vtp, type_label, "Type", mode="both")

        if self._poly_vtp is None:
            self._poly_vtp = hold_vtp
            self._skel_vtp = hold_skel_vtp
        else:
            self._poly_vtp = merge_polys(self._poly_vtp, hold_vtp)
            self._skel_vtp = merge_polys(self._skel_vtp, hold_skel_vtp)

        if is_membrane:
            if self._mbs_vtp is None:
                self._mbs_vtp = hold_vtp
            else:
                self._mbs_vtp = merge_polys(self._mbs_vtp, hold_vtp)

        self._entity_id_counter += 1

    def _collect_network_motifs(self, net, m_type: str, code: str) -> None:
        """Extract per-monomer motifs from a Network and append to the motif array.

        Args:
            net: A Network object with get_pmers_list().
            m_type: Motif type string ("filament", "cprotein", "mbprotein").
            code: Identifier string (e.g. "mt", "actin", "pdb_1bxn").
        """
        rows = []
        label = self._entity_id_counter
        for pmer_id, pmer in enumerate(net.pmers_list):
            for m_id in range(pmer.num_monomers):
                c = pmer.get_mmer_center(m_id)
                q = pmer.get_mmer_rotation(m_id)
                rows.append(
                    (
                        m_type,
                        label,
                        code,
                        pmer_id,
                        c[0],
                        c[1],
                        c[2],
                        q[0],
                        q[1],
                        q[2],
                        q[3],
                    )
                )
        if rows:
            chunk = np.array(rows, dtype=MOTIF_DTYPE)
            self._motifs = np.concatenate((self._motifs, chunk))

    def _collect_membrane_motifs(self, vtp, mb_type: str) -> None:
        """Extract per-point motifs from a membrane VTP surface.

        Surface normals are stored as a pseudo-quaternion (nx, ny, nz, 0).

        Args:
            vtp: A vtkPolyData with point normals.
            mb_type: Membrane type string (e.g. "sphere", "ellipsoid").
        """
        n_points = vtp.GetNumberOfPoints()
        normals = vtp.GetPointData().GetNormals()
        if normals is None:
            logger.warning(
                "Membrane VTP has no normals; skipping membrane motifs."
            )
            return
        label = self._entity_id_counter
        rows = []
        for i in range(n_points):
            x, y, z = vtp.GetPoint(i)
            nx, ny, nz = normals.GetTuple(i)
            rows.append(
                ("membrane", label, mb_type, i, x, y, z, nx, ny, nz, 0.0)
            )
        if rows:
            chunk = np.array(rows, dtype=MOTIF_DTYPE)
            self._motifs = np.concatenate((self._motifs, chunk))

    def _stamp_network(self, net, mask, model) -> np.ndarray:
        """Insert a network into the VOI and density, and return the full-volume label mask.

        Args:
            net: A Network object with insert_density_svol().
            mask: Boolean mask of the monomer subvolume (True = empty/background, False = occupied).
            model: Density subvolume to stamp into the sample density.

        Returns:
            np.ndarray: The full-volume label mask for this component.
        """
        net.insert_density_svol(mask, self._voi, self._v_size, merge="min")
        net.insert_density_svol(
            model, self._density, self._v_size, merge="max"
        )
        hold_lbls = np.zeros(shape=self._shape, dtype=np.float32)
        net.insert_density_svol(
            np.invert(mask), hold_lbls, self._v_size, merge="max"
        )
        return hold_lbls

    def add_set_membranes(
        self, params: dict, max_mbtries: int = 10, grow: int = 0
    ) -> None:
        """Generate and add a set of membranes to the sample.

        Args:
            params (dict): Parameters for the membrane generator class. Should include 'MB_TYPE' key.
            max_mbtries (int, optional): Maximum number of tries to add the membrane. Defaults to 10.
            grow (int, optional): Number of voxels to grow the membrane mask in the VOI. Defaults to 0.

        Raises:
            KeyError: if 'MB_TYPE' key is not in params.
        """

        if "MB_TYPE" not in params:
            raise KeyError(
                "params must include 'MB_TYPE' key specifying the membrane type."
            )

        mb_type = params["MB_TYPE"]
        mb_generator = MbFactory.create(mb_type, params)
        result = mb_generator.generate_set(
            voi=self._voi,
            bg_voi=self._bg_voi,
            v_size=self._v_size,
            max_mbtries=max_mbtries,
            grow=grow,
        )
        if result.num_mbs == 0:
            logger.info("No %s membranes generated.", mb_type)
            return

        logger.info(
            "Inserted %d membranes of type '%s' " "with occupancy %.4f %%.",
            result.num_mbs,
            mb_type,
            100.0 * result.mb_occupancy,
        )

        self._voi = result.voi
        self._density = np.maximum(self._density, result.density)
        self._collect_membrane_motifs(result.vtp, mb_type)
        self._commit_component(
            type_key="membrane",
            model_name=mb_type,
            count=result.num_mbs,
            label_mask=result.mask,
            hold_vtp=result.vtp,
            hold_skel_vtp=result.vtp,
            is_membrane=True,
        )

    def add_helicoidal_network(self, params: dict) -> None:
        """Generate and add a helicoidal fiber network to the sample."""
        hx_type = params["FLMS_TYPE"]

        occ = params["HX_PMER_OCC"]
        if isinstance(occ, (list, tuple)):
            occ = random.uniform(occ[0], occ[1])

        fiber_unit, param_gen, NetworkCls, net_kwargs = FlmsFactory.create(
            hx_type, params, self._v_size
        )
        model_svol = fiber_unit.tomo
        model_surf = fiber_unit.vtp
        model_mask = model_svol < 0.05

        net = NetworkCls(
            voi=self._voi,
            v_size=self._v_size,
            l_length=params["HX_PMER_L"] * params["HX_MMER_RAD"] * 2,
            m_surf=model_surf,
            gen_hfib_params=param_gen,
            occ=occ,
            min_p_len=params["HX_MIN_P_LEN"],
            hp_len=params["HX_HP_LEN"],
            mz_len=params["HX_MZ_LEN"],
            mz_len_f=params["HX_MZ_LEN_F"],
            over_tolerance=params.get("HX_OVER_TOL", 0),
            **net_kwargs,
        )
        if params.get("HX_MIN_NMMER") is not None:
            net.set_min_nmmer(params["HX_MIN_NMMER"])
        net.build_network()

        if net.num_pmers == 0:
            logger.info(
                "No %s fibers generated "
                "(occupancy target may be unreachable).",
                hx_type,
            )
            return

        den_cf_rg = params.get("HX_DEN_CF_RG")
        den_cf = (
            param_gen.gen_den_cf(den_cf_rg[0], den_cf_rg[1])
            if den_cf_rg
            else 1.0
        )
        hold_lbls = self._stamp_network(net, model_mask, model_svol * den_cf)

        type_key = "microtubule" if hx_type == "mt" else "actin"

        # Must run before _commit_component (reads _entity_id_counter before it is incremented)
        self._collect_network_motifs(net, "filament", hx_type)
        self._commit_component(
            type_key=type_key,
            model_name=hx_type,
            count=net.num_pmers,
            label_mask=hold_lbls,
            hold_vtp=net.vtp,
            hold_skel_vtp=net.get_skel(),
        )

    def add_set_cproteins(
        self,
        params: dict,
        data_path: Path,
        surf_dec: float = 0.9,
        mmer_tries: int = 20,
        pmer_tries: int = 100,
    ) -> None:
        """Generate and add a set of cytosolic proteins to the sample. Parameters for the protein generator class should be provided via the params dict.

        Args:
            params (dict): Parameters for the protein generator class. Should include 'type' key.
            data_path (Path): Path to the data directory containing the model files.
            surf_dec (float, optional): Surface decimation factor. Defaults to 0.9.
        """

        pn_generator = PnGen.from_params(
            params, data_path=data_path, surf_dec=surf_dec
        )
        pn_generator.set_scale(self._v_size)

        set_pns = PnSAWLCNet(
            voi=self._voi,
            v_size=self._v_size,
            l_length=pn_generator.pmer_l * pn_generator.surf_diam,
            m_surf=pn_generator.surf,
            max_p_length=pn_generator.pmer_l_max,
            occ=pn_generator.rnd_occ(),
            over_tolerance=pn_generator.over_tolerance,
            poly=None,
            svol=pn_generator.svol,
            tries_mmer=mmer_tries,
            tries_pmer=pmer_tries,
        )

        set_pns.build_network()

        if set_pns.num_mmers == 0:
            logger.info(
                "No %s proteins generated.",
                params.get("MMER_ID", "unknown"),
            )
            return

        hold_lbls = self._stamp_network(
            set_pns, pn_generator.mask, pn_generator.model
        )

        # Must run before _commit_component (reads __entity_id_counter before it is incremented)
        self._collect_network_motifs(
            set_pns, "cprotein", params.get("MMER_ID", "unknown")
        )

        self._commit_component(
            type_key="cprotein",
            model_name=params["MMER_ID"],
            count=set_pns.num_mmers,
            label_mask=hold_lbls,
            hold_vtp=set_pns.vtp,
            hold_skel_vtp=set_pns.get_skel(),
        )

    def add_set_mb_proteins(
        self,
        params: dict,
        data_path: Path,
        surf_dec: float = 0.9,
        mmer_tries: int = 20,
        pmer_tries: int = 100,
    ) -> None:
        """Generate and add a set of membrane-bound proteins to the sample.

        Requires membranes to have been added first. Proteins are placed on
        membrane surfaces using a constrained SAWLC network.

        Args:
            params (dict): Parameters from the .pms config file. Must include 'MMER_ID'.
            data_path (Path): Path to the data directory containing model files.
            surf_dec (float, optional): Surface decimation factor. Defaults to 0.9.
            mmer_tries (int, optional): Monomer placement attempts per step. Defaults to 20.
            pmer_tries (int, optional): Polymer placement attempts. Defaults to 100.
        """
        if self._mbs_vtp is None:
            logger.warning("No membrane surfaces; skipping membrane protein.")
            return

        pm_gen = PmGen.from_params(
            params, data_path=data_path, surf_dec=surf_dec, v_size=self._v_size
        )
        pm_gen.set_scale(self._v_size)

        net = PmSAWLCPolyNet(
            voi=self._voi,
            v_size=self._v_size,
            l_length=pm_gen.pmer_l * pm_gen.surf_diam,
            m_surf=pm_gen.surf,
            max_p_length=pm_gen.pmer_l_max,
            occ=pm_gen.rnd_occ(),
            poly=self._mbs_vtp,
            reverse_normals=pm_gen.reverse_normals,
            over_tolerance=pm_gen.over_tolerance,
            svol=pm_gen.svol,
            tries_mmer=mmer_tries,
            tries_pmer=pmer_tries,
        )
        net.build_network()

        if net.num_mmers == 0:
            logger.info(
                "No %s membrane proteins generated.",
                params.get("MMER_ID", "unknown"),
            )
            return

        hold_lbls = self._stamp_network(net, pm_gen.mask, pm_gen.model)

        # Must run before _commit_component (reads _entity_id_counter before it is incremented)
        self._collect_network_motifs(
            net, "mbprotein", params.get("MMER_ID", "unknown")
        )

        self._commit_component(
            type_key="mbprotein",
            model_name=params["MMER_ID"],
            count=net.num_mmers,
            label_mask=hold_lbls,
            hold_vtp=net.vtp,
            hold_skel_vtp=net.get_skel(),
        )

    def print_summary(self) -> None:
        """Prints a summary of the sample contents.

        Returns:
            None
        """
        message = "Synthetic Sample Summary:\n"
        message += f"  Shape: {self._shape}\n"
        message += f"  Voxel Size: {self._v_size} Å\n"
        message += f"  VOI Voxels: {self._voi_voxels}\n"
        message += "  Structure Counts:\n"
        for struct_type, count in self._structure_counts.items():
            message += f"    {struct_type}:\n"
            message += f"      Count: {count}\n"
            message += f"      Occupancy: {100.0 * self._voxel_counts.get(struct_type, 0) / self._voi_voxels:.4f} %\n"
        logger.info(message)
