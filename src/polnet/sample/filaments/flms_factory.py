"""Factory for helical filament networks.

:class:`FlmsFactory` instantiates the correct fiber unit,
parameter generator, and network class (microtubule ``"mt"`` or
actin ``"actin"``) from a parsed ``.flms`` configuration dict.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

from .fiber_unit import (
    FiberUnitSDimer,
    MTUnit,
)
from .flms_gen import (
    FlmsParamGen,
    HxParamGenBranched,
)
from .flms_network import (
    NetHelixFiber,
    NetHelixFiberB,
)


class FlmsFactory:
    """Factory for creating helical filament network components.

    Instantiates the correct fiber unit, parameter generator, and
    network class from a parsed ``.flms`` configuration dict.
    """

    @classmethod
    def create(cls, hx_type: str, params: dict, v_size: float):
        """Create filament components for the given helix type.

        Args:
            hx_type (str): Helix type identifier (``"mt"`` for
                microtubule, ``"actin"`` for actin).
            params (dict): Parsed ``.flms`` configuration dict.
            v_size (float): Voxel size in Angstroms.

        Returns:
            tuple: A 4-tuple ``(fiber_unit, param_gen, NetworkCls,
            net_kwargs)`` where *fiber_unit* is the structural
            repeat, *param_gen* the stochastic parameter
            generator, *NetworkCls* the network class to
            instantiate, and *net_kwargs* extra keyword arguments
            for its constructor.

        Raises:
            ValueError: If *hx_type* is not ``"mt"`` or
                ``"actin"``.
        """
        if hx_type == "mt":
            fiber_unit = MTUnit(
                sph_rad=params["HX_MMER_RAD"],
                mt_rad=params["MT_RAD"],
                n_units=int(params["MT_NUNITS"]),
                v_size=v_size,
            )
            param_gen = FlmsParamGen()
            net_cls = NetHelixFiber
            net_kwargs = dict(
                unit_diam=(params["MT_RAD"] + 0.5 * params["HX_MMER_RAD"])
                * 2.4,
            )
        elif hx_type == "actin":
            fiber_unit = FiberUnitSDimer(
                sph_rad=params["HX_MMER_RAD"],
                v_size=v_size,
            )
            param_gen = HxParamGenBranched()
            net_cls = NetHelixFiberB
            net_kwargs = dict(
                b_prop=params["A_BPROP"],
                max_p_branch=int(params["A_MAX_P_BRANCH"]),
            )
        else:
            raise ValueError(f"Unsupported helix type: {hx_type}")

        return fiber_unit, param_gen, net_cls, net_kwargs
