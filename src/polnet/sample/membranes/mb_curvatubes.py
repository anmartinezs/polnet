"""Phase-field membrane generator using the curvatubes library.

Integrates the curvatubes library (Anna Song, MIT License) for
generating tubular and membranous shape textures via phase-field
gradient descent on a curvature functional.

The curvature functional is:

.. math::

    F(S) = \\int \\{ a_{2,0} \\kappa_1^2 + a_{1,1} \\kappa_1 \\kappa_2
           + a_{0,2} \\kappa_2^2 + b_{1,0} \\kappa_1 + b_{0,1} \\kappa_2
           + c \\} \\, dH^2

Supports two parameterisations:

- **polynomial**: Direct specification of a20, a11, a02, b10, b01, c.
- **helfrich**: Specification via bending rigidity (kappa_b), Gaussian
  curvature modulus (kappa_G), and spontaneous curvature (H0).

References:
    Song, A. (2021). "Generation of Tubular and Membranous Shape
    Textures with Curvature Functionals." Journal of Mathematical
    Imaging and Vision, 64(1), 17--40.
    https://doi.org/10.1007/s10851-021-01049-9

    Source: https://github.com/annasongmaths/curvatubes
    License: MIT -- Copyright (c) 2021 Anna SONG

:maintainer: Juan Diego Gallego Nicolas
"""

import contextlib
import gc
import io
import os
import random
import shutil
import tempfile
import threading
from pathlib import Path

import numpy as np
import vtk
from scipy.ndimage import (
    binary_dilation,
    distance_transform_edt,
    gaussian_filter,
)

from .mb import (
    Mb,
    MbError,
    MbGen,
    MbSetResult,
)
from .mb_factory import MbFactory
from ...logging_conf import _LOGGER as logger
from ...utils.affine import poly_scale
from ...utils.poly import add_sfield_to_poly
from ...utils.utils import (
    density_norm,
    iso_surface,
    lin_map,
    poly_threshold,
)


def _helfrich_to_polynomial(
    kappa_b: float,
    kappa_G: float,
    H0: float,
) -> tuple[float, float, float, float, float, float]:
    """Convert Helfrich parameters to polynomial coefficients.

    The Helfrich functional is:

    .. math::

        F_H = \\int \\{ \\kappa_b / 2 \\cdot (H - H_0)^2
              + \\kappa_G \\cdot K \\} \\, dA

    where :math:`H = (\\kappa_1 + \\kappa_2) / 2` (mean curvature)
    and :math:`K = \\kappa_1 \\kappa_2` (Gaussian curvature).

    In the curvatubes convention the sum curvature
    :math:`H_{cv} = \\kappa_1 + \\kappa_2 = 2H` is used, giving
    the mapping:

    .. math::

        a_{2,0} = a_{0,2} = \\kappa_b / 8 \\\\
        a_{1,1} = \\kappa_b / 4 + \\kappa_G \\\\
        b_{1,0} = b_{0,1} = -\\kappa_b \\cdot H_0 / 2 \\\\
        c = \\kappa_b \\cdot H_0^2 / 2

    Args:
        kappa_b: Bending rigidity.
        kappa_G: Gaussian curvature modulus.
        H0: Spontaneous curvature (standard convention).

    Returns:
        Tuple ``(a20, a11, a02, b10, b01, c)``.
    """
    a20 = kappa_b / 8.0
    a02 = kappa_b / 8.0
    a11 = kappa_b / 4.0 + kappa_G
    b10 = -kappa_b * H0 / 2.0
    b01 = b10
    c = kappa_b * H0**2 / 2.0
    return a20, a11, a02, b10, b01, c


@contextlib.contextmanager
def _suppress_curvatubes_output():
    """Suppress all display, print, and file output from curvatubes.

    Creates a temporary directory for snapshot files and patches
    ``matplotlib.pyplot.show`` to a no-op.

    Yields:
        str: Path to a temporary snapshot directory (trailing separator).
    """
    import matplotlib.pyplot as plt

    orig_show = plt.show
    plt.show = lambda *_a, **_kw: None

    tmpdir = Path(tempfile.mkdtemp(prefix="polnet_cvtub_"))
    (tmpdir / "Curves").mkdir(parents=True, exist_ok=True)

    with contextlib.redirect_stdout(io.StringIO()):
        try:
            yield str(tmpdir) + os.sep
        finally:
            plt.show = orig_show
            plt.close("all")
            shutil.rmtree(tmpdir, ignore_errors=True)


try:
    from external.cvtub import generator as _cvtub_gen
    from external.cvtub import utils as _cvtub_utils
except ImportError as exc:
    raise ImportError(
        "The vendored curvatubes library was not found at "
        "'external.cvtub'.  Ensure the 'src/external/cvtub/' "
        "package exists and that 'src/' is on the Python path."
    ) from exc


@MbFactory.register("curvatubes")
class CurvatubesGen(MbGen):
    """Phase-field membrane generator using curvatubes.

    Generates a membrane that fills the entire tomogram volume via
    gradient descent on a curvature functional.  The zero level set
    of the optimised phase field defines the membrane midplane.

    Unlike the other generators, this one overrides
    :meth:`generate_set` because the curvatubes optimisation is a
    single expensive GPU pass that produces the complete membrane
    rather than individual localisable units.

    References:
        Song, A. (2021). JMIV, 64(1), 17--40.
        MIT License -- Copyright (c) 2021 Anna SONG.
    """

    # pylint: disable=too-many-instance-attributes

    def __init__(
        self,
        thick_rg: tuple[float, float],
        layer_s_rg: tuple[float, float],
        occ_rg: tuple[float, float],
        over_tol: float,
        mb_den_cf_rg: tuple[float, float],
        *,
        # curvatubes phase-field
        eps: float,
        delta_x: float,
        xi: float,
        M0: float,
        flow_type: str,
        mode: str,
        # optimiser
        optim_method: str,
        maxeval: int,
        lr: float,
        sigma_blur: float,
        # polynomial coefficients
        a20: float,
        a11: float,
        a02: float,
        b10: float,
        b01: float,
        c: float,
        # adam-specific
        betas: tuple[float, float] = (0.9, 0.999),
        eps_adam: float = 1e-8,
        weight_decay: float = 0.0,
        amsgrad: bool = False,
        # bfgs-specific
        bfgs_max_iter: int = 20,
        history_size: int = 100,
        line_search_fn: str | None = "strong_wolfe",
        # tracking
        display_it_nb: int = 500,
        fill_curve_nb: int = 10,
        inform_every: int = 500,
        # initialisation
        init_prop: float | None = None,
    ) -> None:
        """Constructor.

        Args:
            thick_rg: ``(min, max)`` bilayer thickness in angstroms.
            layer_s_rg: ``(min, max)`` Gaussian layer sigma in angstroms.
            occ_rg: ``(min, max)`` target occupancy fraction
                (kept for interface compatibility; not used by the
                occupancy loop which is bypassed).
            over_tol: Overlap tolerance fraction (not used).
            mb_den_cf_rg: ``(min, max)`` density contrast factor.
            eps: Phase-field diffuse interface width (epsilon).
            delta_x: Mathematical length per pixel.
            xi: Regularisation parameter (typically ``1e-6``).
            M0: Target average mass level (``-1`` to ``1``).
                Ignored when *flow_type* is ``'averm0'`` (mass is
                inferred from the initial condition).
            flow_type: Gradient flow — ``'L2'``, ``'averm0'``,
                or ``'cons'``.
            mode: Boundary conditions — ``'periodic'`` or
                ``'replicate'``.
            optim_method: Optimiser — ``'adam'`` or ``'bfgs'``.
            maxeval: Maximum number of evaluations.
            lr: Learning rate.
            sigma_blur: Gaussian blur sigma (pixels) applied to
                the phase field at every step.
            a20, a11, a02, b10, b01, c: Polynomial curvature
                coefficients.
            betas: Adam betas.
            eps_adam: Adam epsilon.
            weight_decay: Adam weight decay.
            amsgrad: Adam AMSGrad flag.
            bfgs_max_iter: L-BFGS max iterations per step.
            history_size: L-BFGS history size.
            line_search_fn: L-BFGS line search function.
            display_it_nb: Internal tracking interval (iterations).
            fill_curve_nb: Energy-curve sampling interval.
            inform_every: How often (in evaluations) to log
                optimisation progress to the polnet logger.
                Set to ``0`` to disable progress logging.
            init_prop: Fraction of ``+1`` pixels in the random
                initialisation.  Defaults to ``(1 + M0) / 2`` so
                the initial mean matches the target mass.
        """
        super().__init__(
            thick_rg=thick_rg,
            layer_s_rg=layer_s_rg,
            occ_rg=occ_rg,
            over_tol=over_tol,
            mb_den_cf_rg=mb_den_cf_rg,
        )
        self._eps = eps
        self._delta_x = delta_x
        self._xi = xi
        self._M0 = M0
        self._flow_type = flow_type
        self._mode = mode
        self._optim_method = optim_method
        self._maxeval = maxeval
        self._lr = lr
        self._sigma_blur = sigma_blur
        self._a20 = a20
        self._a11 = a11
        self._a02 = a02
        self._b10 = b10
        self._b01 = b01
        self._c = c
        self._betas = betas
        self._eps_adam = eps_adam
        self._weight_decay = weight_decay
        self._amsgrad = amsgrad
        self._bfgs_max_iter = bfgs_max_iter
        self._history_size = history_size
        self._line_search_fn = line_search_fn
        self._display_it_nb = display_it_nb
        self._fill_curve_nb = fill_curve_nb
        self._inform_every = inform_every
        self._init_prop = (
            init_prop if init_prop is not None else (1.0 + M0) / 2.0
        )

    @classmethod
    def from_params(cls, params: dict) -> "CurvatubesGen":
        """Create a :class:`CurvatubesGen` from a parameter dictionary.

        Supports both ``'polynomial'`` and ``'helfrich'``
        parameterisations.  When ``CT_PARAMETERIZATION`` is
        ``'helfrich'``, polynomial coefficients are derived from
        ``CT_KAPPA_B``, ``CT_KAPPA_G``, and ``CT_H0``.

        Args:
            params: Dictionary of membrane parameters (typically
                loaded from a ``.mbs`` file).

        Returns:
            CurvatubesGen instance.
        """
        parameterization = params.get("CT_PARAMETERIZATION", "polynomial")
        if parameterization == "helfrich":
            kappa_b = params["CT_KAPPA_B"]
            kappa_G = params.get("CT_KAPPA_G", 0.0)
            H0 = params.get("CT_H0", 0.0)
            a20, a11, a02, b10, b01, c_val = _helfrich_to_polynomial(
                kappa_b, kappa_G, H0
            )
        elif parameterization == "polynomial":
            a20 = params.get("CT_A20", 1.0)
            a11 = params.get("CT_A11", 0.0)
            a02 = params.get("CT_A02", 1.0)
            b10 = params.get("CT_B10", 0.0)
            b01 = params.get("CT_B01", 0.0)
            c_val = params.get("CT_C", 0.0)
        else:
            raise ValueError(
                f"Unknown CT_PARAMETERIZATION: {parameterization!r}. "
                "Must be 'polynomial' or 'helfrich'."
            )

        M0 = params.get("CT_M0", -0.5)

        return cls(
            thick_rg=params.get("MB_THICK_RG", (25.0, 35.0)),
            layer_s_rg=params.get("MB_LAYER_S_RG", (0.5, 2.0)),
            occ_rg=params.get("MB_OCC_RG", (0.001, 0.003)),
            over_tol=params.get("MB_OVER_TOL", 0.0),
            mb_den_cf_rg=params.get("MB_DEN_CF_RG", (0.3, 0.5)),
            eps=params.get("CT_EPS", 0.15),
            delta_x=params.get("CT_DELTA_X", 0.03),
            xi=params.get("CT_XI", 1e-6),
            M0=M0,
            flow_type=params.get("CT_FLOW_TYPE", "averm0"),
            mode=params.get("CT_MODE", "periodic"),
            optim_method=params.get("CT_OPTIM_METHOD", "adam"),
            maxeval=params.get("CT_MAXEVAL", 5000),
            lr=params.get("CT_LR", 0.05),
            sigma_blur=params.get("CT_SIGMA_BLUR", 1.0),
            a20=a20,
            a11=a11,
            a02=a02,
            b10=b10,
            b01=b01,
            c=c_val,
            betas=params.get("CT_BETAS", (0.9, 0.999)),
            eps_adam=params.get("CT_EPS_ADAM", 1e-8),
            weight_decay=params.get("CT_WEIGHT_DECAY", 0.0),
            amsgrad=params.get("CT_AMSGRAD", False),
            bfgs_max_iter=params.get("CT_BFGS_MAX_ITER", 20),
            history_size=params.get("CT_HISTORY_SIZE", 100),
            line_search_fn=params.get("CT_LINE_SEARCH_FN", "strong_wolfe"),
            display_it_nb=params.get("CT_DISPLAY_IT_NB", 500),
            fill_curve_nb=params.get("CT_FILL_CURVE_NB", 10),
            inform_every=params.get("CT_INFORM_EVERY", 500),
            init_prop=params.get("CT_INIT_PROP", None),
        )

    def _build_optim_props(self) -> dict:
        """Build the optimiser-properties dict expected by curvatubes."""
        props: dict = {
            "maxeval": self._maxeval,
            "sigma_blur": self._sigma_blur,
            "display_it_nb": self._display_it_nb,
            "fill_curve_nb": self._fill_curve_nb,
            "lr": self._lr,
        }
        if self._optim_method == "adam":
            props.update(
                {
                    "betas": self._betas,
                    "eps_adam": self._eps_adam,
                    "weight_decay": self._weight_decay,
                    "amsgrad": self._amsgrad,
                }
            )
        elif self._optim_method == "bfgs":
            props.update(
                {
                    "bfgs_max_iter": self._bfgs_max_iter,
                    "history_size": self._history_size,
                    "line_search_fn": self._line_search_fn,
                }
            )
        return props

    def _run_curvatubes(self, voi_shape: tuple[int, int, int]) -> np.ndarray:
        """Run the curvatubes phase-field optimisation.

        Args:
            voi_shape: Volume shape ``(Z, X, Y)`` in voxels.

        Returns:
            numpy.ndarray: Optimised phase field *u* (float32).

        Raises:
            MbError: If CUDA is unavailable or the optimisation
                diverges.
        """
        import torch

        if not torch.cuda.is_available():
            raise MbError(
                "Curvatubes requires a CUDA-capable GPU.  "
                "torch.cuda.is_available() returned False."
            )

        params = (
            self._eps,
            self._a20,
            self._a11,
            self._a02,
            self._b10,
            self._b01,
            self._c,
        )
        optim_props = self._build_optim_props()

        # Random initialisation
        if self._flow_type == "cons":
            v0 = torch.stack(
                [
                    _cvtub_utils.random_init(voi_shape, prop=self._init_prop),
                    _cvtub_utils.random_init(voi_shape, prop=self._init_prop),
                    _cvtub_utils.random_init(voi_shape, prop=self._init_prop),
                ]
            )
        else:
            v0 = _cvtub_utils.random_init(voi_shape, prop=self._init_prop)

        # For 'averm0' the mass is set by the initial condition;
        # curvatubes requires M0=None in that case.
        M0_val = None if self._flow_type == "averm0" else self._M0

        logger.info(
            "Running curvatubes optimisation "
            "(%s, %s, %s) on %s grid with maxeval=%d ...",
            self._optim_method,
            self._flow_type,
            self._mode,
            voi_shape,
            self._maxeval,
        )

        # Progress
        stop_event = threading.Event()
        inform_every = self._inform_every
        maxeval = self._maxeval

        def _monitor() -> None:
            """Periodically log curvatubes progress."""
            last_reported = -1
            while not stop_event.is_set():
                n = getattr(_cvtub_gen, "n_evals", None)
                if n is not None and inform_every > 0:
                    step = (n // inform_every) * inform_every
                    if step > last_reported and n > 0:
                        last_reported = step
                        E_vals = getattr(_cvtub_gen, "E_curve", [])
                        E_str = (
                            f", E={E_vals[-1]:.4e}" if len(E_vals) > 0 else ""
                        )
                        logger.info(
                            "  curvatubes progress: " "eval %d / %d%s",
                            n,
                            maxeval,
                            E_str,
                        )
                stop_event.wait(2.0)  # poll every 2 s

        monitor = threading.Thread(target=_monitor, daemon=True)
        monitor.start()

        try:
            with _suppress_curvatubes_output() as snapshot_folder:
                u = _cvtub_gen._generate_shape(
                    v0=v0,
                    params=params,
                    delta_x=self._delta_x,
                    xi=self._xi,
                    optim_method=self._optim_method,
                    optim_props=optim_props,
                    flow_type=self._flow_type,
                    mode=self._mode,
                    M0=M0_val,
                    snapshot_folder=snapshot_folder,
                    exp_title="polnet_",
                    cond_take_snapshot=None,
                    display_all=False,
                    return_var=False,
                    return_energy=False,
                    check_viable=False,
                )
        finally:
            stop_event.set()
            monitor.join(timeout=5.0)

        u_np = u.detach().cpu().numpy().astype(np.float32)

        # Free GPU memory — both the returned tensor and any
        # module-level globals that curvatubes may have leaked.
        del u, v0

        # Null out curvatubes module-level globals that may still
        # hold GPU tensor references (see VRAM audit in generator.py).
        for _attr in ("u", "uu", "params2", "M02"):
            if hasattr(_cvtub_gen, _attr):
                setattr(_cvtub_gen, _attr, None)

        gc.collect()
        torch.cuda.empty_cache()

        if np.any(np.isnan(u_np)):
            raise MbError(
                "Curvatubes optimisation diverged (NaN values).  "
                "Try reducing the learning rate or adjusting "
                "the regularisation parameter (xi)."
            )

        logger.info("Curvatubes optimisation complete.")
        return u_np

    def _build(self, voi_shape: tuple[int, int, int], v_size: float) -> Mb:
        """Generate a phase-field membrane filling the entire volume.

        The optimised phase field *u* ranges approximately from
        ``-1`` to ``+1``.  The membrane midplane is the ``u = 0``
        level set.  A geometric bilayer profile is built from it:

        1. Binarise at ``u = 0`` and extract a 1-voxel boundary
           shell via binary dilation.
        2. Compute an EDT (Euclidean distance transform) from the
           shell.
        3. **mask** -- all voxels within half the physical bilayer
           thickness of the shell (``EDT < t_v``).
        4. **surface** -- marching-cubes iso-surface on the smooth
           phase field at ``u = 0`` (sub-voxel accuracy).
        5. **density** -- ``mask − shell``, giving two leaflet
           bands separated by a dark midplane, Gaussian-smoothed
           and normalised.

        This follows the algorithm described in Gallego-Nicolas
        et al. (2026), bioRxiv 2026.01.15.699326.

        Args:
            voi_shape: Volume shape in voxels.
            v_size: Voxel size in angstroms.

        Returns:
            Mb: Constructed membrane.

        Raises:
            MbError: If the optimisation fails or no valid surface
                is produced.
        """
        thick = random.uniform(self._thick_rg[0], self._thick_rg[1])
        layer_s = random.uniform(self._layer_s_rg[0], self._layer_s_rg[1])

        # Convert physical parameters to voxel units
        t_v = 0.5 * thick / v_size  # half-thickness in voxels
        s_v = layer_s / v_size  # Gaussian sigma in voxels

        if t_v < 1.5:
            logger.warning(
                "Bilayer half-thickness t_v=%.2f vx is below "
                "1.5 vx (thick=%.1f A, v_size=%.1f A).  "
                "Clamping to 1.5 vx for a valid bilayer profile.",
                t_v,
                thick,
                v_size,
            )
            t_v = 1.5

        u_np = self._run_curvatubes(voi_shape)

        # ---- 1-voxel boundary shell at the zero level set ----
        # binary_dilation with the default structuring element
        # (6-connectivity / face-neighbours) produces the thinnest
        # possible shell on the u < 0 side of the boundary.
        ct_bin = u_np > 0
        ct_shell = binary_dilation(ct_bin).astype(np.float32) - ct_bin.astype(
            np.float32
        )

        if ct_shell.sum() == 0:
            raise MbError(
                "Curvatubes produced no interface region "
                "(phase field is entirely in a single bulk "
                "phase; the zero level set is empty)."
            )

        # ---- distance transform from the shell ----
        ct_dist = distance_transform_edt(1.0 - ct_shell)

        # ---- mask (band of half-thickness t_v around shell) ----
        mask = ct_dist < t_v

        if mask.sum() == 0:
            raise MbError(
                "Membrane mask is empty after applying the "
                f"thickness band (thick={thick:.1f} A, "
                f"t_v={t_v:.2f} vx).  The bilayer thickness "
                "may be too small for the voxel size."
            )

        # ---- surface (zero level set on the smooth field) ----
        # We extract the iso-surface from the smooth phase field
        # (not the blocky distance transform) for sub-voxel
        # triangulation quality.
        try:
            surf = iso_surface(u_np, 0.0)
        except Exception as exc:
            raise MbError(
                "Failed to extract iso-surface from the phase " f"field: {exc}"
            ) from exc

        if surf.GetNumberOfPoints() == 0:
            raise MbError("Iso-surface extraction yielded no points.")

        add_sfield_to_poly(
            surf,
            mask,
            "mb_mask",
            dtype="int",
            interp="NN",
            mode="points",
        )
        surf = poly_threshold(surf, "mb_mask", mode="points", low_th=0.5)

        # ---- density (bilayer profile) ----
        # mask − shell gives +1 on both leaflet volumes (inner
        # and outer) and 0 on the 1-voxel midplane and outside.
        # After Gaussian smoothing, this produces the
        # characteristic bright-dark-bright profile of a
        # phospholipid bilayer in cryo-EM.
        G = mask.astype(np.float32) - ct_shell

        density = lin_map(
            density_norm(gaussian_filter(G, s_v), inv=True),
            ub=0,
            lb=1,
        )

        return Mb(
            voi_shape=voi_shape,
            v_size=v_size,
            thick=thick,
            layer_s=layer_s,
            density=density,
            mask=mask,
            surf=surf,
        )

    def generate_set(
        self,
        *,
        voi: np.ndarray,
        v_size: float,
        bg_voi: np.ndarray | None = None,
        max_mbtries: int = 10,
        grow: int = 0,
    ) -> MbSetResult:
        """Build a single full-volume curvatubes membrane.

        Overrides :meth:`MbGen.generate_set` because curvatubes
        membranes fill the entire tomogram in a single expensive
        GPU optimisation pass — the iterative occupancy loop of
        the base class is not applicable.

        Args:
            voi: 3-D bool array (volume of interest, copied).
            v_size: Voxel size in angstroms.
            bg_voi: Optional background mask for additional
                masking.
            max_mbtries: Unused (interface compatibility).
            grow: VOI dilation iterations for the inserted
                membrane.

        Returns:
            MbSetResult with the membrane data.
        """
        shape = voi.shape

        try:
            mb = self._build(shape, v_size)
        except MbError as exc:
            logger.error("Curvatubes membrane generation failed: %s", exc)
            return MbSetResult(
                voi=voi.copy(),
                density=np.zeros(shape, dtype=np.float32),
                mask=np.zeros(shape, dtype=bool),
                vtp=vtk.vtkPolyData(),
                num_mbs=0,
                mb_occupancy=0.0,
            )

        if bg_voi is not None:
            mb.masking(bg_voi)

        if mb.vol <= 0:
            logger.warning(
                "Curvatubes membrane has zero volume after masking."
            )
            return MbSetResult(
                voi=voi.copy(),
                density=np.zeros(shape, dtype=np.float32),
                mask=np.zeros(shape, dtype=bool),
                vtp=vtk.vtkPolyData(),
                num_mbs=0,
                mb_occupancy=0.0,
            )

        cf = self.rnd_cf()
        density = np.zeros(shape, dtype=np.float32)
        mask_acc = np.zeros(shape, dtype=bool)
        voi_out = voi.copy()

        mb.insert_density_svol(density, merge="max", mode="tomo")
        mb.insert_density_svol(mask_acc, merge="max", mode="mask")
        mb.insert_density_svol(voi_out, merge="min", mode="voi", grow=grow)

        density *= cf
        surfs = poly_scale(mb.vtp, v_size)

        voi_sum = voi.sum()
        occ = mask_acc.sum() / voi_sum if voi_sum > 0 else 0.0

        logger.info(
            "Curvatubes membrane inserted (occupancy %.4f %%).",
            100.0 * occ,
        )

        return MbSetResult(
            voi=voi_out,
            density=density,
            mask=mask_acc,
            vtp=surfs,
            num_mbs=1,
            mb_occupancy=occ,
        )
