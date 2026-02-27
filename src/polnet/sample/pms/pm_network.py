"""SAWLC polymer network constrained to a membrane surface.

Defines :class:`NetSAWLCPoly`, the membrane-bound counterpart
of :class:`~polnet.sample.pns.pn_network.NetSAWLC`.  Monomers
are placed by sampling insertion points directly on the VTK
PolyData surface that represents the lipid bilayer.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

from ..pns.pn_network import NetSAWLC
from ...utils.poly import poly_reverse_normals


class NetSAWLCPoly(NetSAWLC):
    """SAWLC network where monomers are constrained to a polydata surface.

    This is the membrane-bound counterpart of :class:`NetSAWLC`.
    The ``poly`` parameter is required (cannot be None).
    """

    def __init__(
        self,
        voi,
        v_size,
        l_length,
        m_surf,
        max_p_length,
        occ,
        poly,
        reverse_normals=True,
        over_tolerance=0,
        svol=None,
        tries_mmer=50,
        tries_pmer=10,
    ):
        if poly is None:
            raise ValueError("NetSAWLCPoly requires a membrane poly surface.")
        if reverse_normals:
            poly = poly_reverse_normals(poly)
        super().__init__(
            voi=voi,
            v_size=v_size,
            l_length=l_length,
            m_surf=m_surf,
            max_p_length=max_p_length,
            occ=occ,
            over_tolerance=over_tolerance,
            poly=poly,
            svol=svol,
            tries_mmer=tries_mmer,
            tries_pmer=tries_pmer,
        )
