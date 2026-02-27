#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 9 11:37:02 2020

@author: Anna SONG

Computes the curvature diagrams of a surface {u = 0} defined by a phase-field u
with profile close to a tanh profile. Curvature diagrams are histograms of the
diffuse principal curvatures (kap_{1,eps}, kap_{2,eps}) on the surface, computed
with respect to a triangular mesh of the surface.

The statistics of these curvatures characterize the *texture* of a 3D shape.

Source: https://github.com/annasongmaths/curvatubes
License: MIT — Copyright (c) 2021 Anna SONG
See LICENSE in this directory.

NOTE (polnet integration): Modifications from the upstream source:

1. **Relative imports** — absolute ``from cvtub.xxx`` imports converted to
   relative ``from .xxx`` so the package works when vendored under
   ``src/external/cvtub/``.  Each changed line is marked with
   ``# polnet: relative import``.

2. **Code formatting** — reformatted with ``black`` (line-length 79) on
   2026-02-26 for consistency with the rest of the polnet codebase.
   No functional changes.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn

import skimage
from scipy.interpolate import RegularGridInterpolator as RGI

import torch

dtype = (
    torch.cuda.FloatTensor
)  # if torch.cuda.is_available() else torch.FloatTensor

from .utils import (
    double_W_prime,
    manual_softplus,
)  # polnet: relative import (was ``from cvtub.utils ...``)
from .energy import (
    auxiliary_function,
)  # polnet: relative import (was ``from cvtub.energy ...``)


def kap_eps(u, eps=0.02, delta_x=0.01, mode="periodic", xi=1e-6):

    if type(u) != torch.Tensor:
        u = torch.tensor(u).type(dtype)
    grad, Hess_diag, Hess_off = auxiliary_function(u, eps, delta_x, mode)

    dz_u = grad[..., 0]
    dx_u = grad[..., 1]
    dy_u = grad[..., 2]
    Hzz_u = Hess_diag[..., 0]
    Hxx_u = Hess_diag[..., 1]
    Hyy_u = Hess_diag[..., 2]
    Hzx_u = Hess_off[..., 0]
    Hzy_u = Hess_off[..., 1]
    Hxy_u = Hess_off[..., 2]

    norm_grad_sq = (grad**2).sum(-1)
    ngd_u = (norm_grad_sq + xi**2).sqrt()
    Wprim = double_W_prime(u)

    # nHn is n^T H n where H is the hessian of u and n approximates the direction
    # of the gradient (not well defined when grad u ~ 0) up to a small term xi
    nHn = (
        Hzz_u * dz_u**2
        + Hxx_u * dx_u**2
        + Hyy_u * dy_u**2
        + 2 * (Hzx_u * dz_u * dx_u + Hzy_u * dz_u * dy_u + Hxy_u * dx_u * dy_u)
    )
    nHn /= norm_grad_sq + xi**2

    Tra = -eps * (Hzz_u + Hxx_u + Hyy_u) + Wprim / eps

    Nor_pow2 = (
        (eps**2) * ((Hess_diag**2).sum(-1) + 2 * (Hess_off**2).sum(-1))
        + Wprim**2 / (eps**2)
        - 2 * Wprim * nHn
    )

    positive_part = manual_softplus(2 * Nor_pow2 - Tra**2)
    sqrt_part = torch.sqrt(positive_part)

    # Heps = Tra / (eps * ngd_u)
    # Keps = (Tra**2 - Nor_pow2) / (2 * eps**2 * (norm_grad_sq + xi**2) )
    kap1_eps = (Tra + sqrt_part) / (2 * eps * ngd_u)
    kap2_eps = (Tra - sqrt_part) / (2 * eps * ngd_u)
    # print('kap1_eps * kap2_eps may be different from Keps')

    kap1_eps = kap1_eps.detach().cpu().numpy()
    kap2_eps = kap2_eps.detach().cpu().numpy()

    return kap1_eps, kap2_eps


def curvhist(
    vol,
    kap1_eps,
    kap2_eps,
    lev=0,
    delta_x=0.01,
    show_figs=True,
    bins=100,
    save=False,
    save_name="default.png",
):

    Z, X, Y = vol.shape
    print("Analysis for level set u = {}".format(lev))
    verts, faces, normals, values_mc = skimage.measure.marching_cubes(
        vol, level=lev
    )  # , allow_degenerate = False)

    # compute cells areas
    As = verts[faces[:, 0]]
    Bs = verts[faces[:, 1]]
    Cs = verts[faces[:, 2]]
    Gs = (As + Bs + Cs) / 3  # barycenters of the mesh cells
    alen = np.linalg.norm(Cs - Bs, axis=1)
    blen = np.linalg.norm(Cs - As, axis=1)
    clen = np.linalg.norm(Bs - As, axis=1)
    midlen = 0.5 * (alen + blen + clen)
    areas = np.sqrt(
        midlen * (midlen - alen) * (midlen - blen) * (midlen - clen)
    )

    total_area = areas.sum()
    print(
        "Total area of the mesh is ", total_area * delta_x**2, " (real unit)"
    )

    total_volume = (vol > lev).sum()
    print(
        "Total volume of the solid: ",
        total_volume * delta_x**3,
        "approximately (real unit),",
    )
    print(
        "i.e. ", total_volume / (Z * X * Y), "% of the total domain volume. \n"
    )

    # interpolate the values of kap1_eps and kap2_eps on the barycenters of each mesh cell
    aux_grid = (np.arange(Z), np.arange(X), np.arange(Y))

    interpolator_kap1 = RGI(aux_grid, kap1_eps, method="linear")
    kap1_vals = interpolator_kap1(Gs)

    interpolator_kap2 = RGI(aux_grid, kap2_eps, method="linear")
    kap2_vals = interpolator_kap2(Gs)

    if show_figs:
        "plot histogram bins --- probability densities"

        fig, ax = plt.subplots(2, 2, figsize=(12, 12))

        heights_kap1, bins_kap1, _ = ax[1, 0].hist(
            kap1_vals, bins=bins, density=True, weights=areas
        )
        ax[1, 0].set_title("kap1", y=-0.15, fontsize=20)
        heights_kap2, bins_kap2, _ = ax[1, 1].hist(
            kap2_vals, bins=bins, density=True, weights=areas
        )
        ax[1, 1].set_title("kap2", y=-0.15, fontsize=20)

        if save:
            plt.savefig(save_name)
        plt.show()

    return kap1_vals, kap2_vals, areas


def density_scatter(
    x,
    y,
    areas,
    xlabel="",
    ylabel="",
    showid=False,
    showparab=False,
    equalaxis=False,
    size=0.1,
    bins=20,
    sort=True,
    showfig=True,
    save=False,
    save_name="histogram.png",
    **kwargs,
):
    "Plots a beautiful colored scatter plot of a 2D histogram (such as curvature diagrams)"
    # this code comes from StackOverflow

    fig, ax = plt.subplots(figsize=(8, 6))
    data, x_e, y_e = np.histogram2d(
        x, y, weights=areas, bins=bins, density=True
    )
    z = interpn(
        (0.5 * (x_e[1:] + x_e[:-1]), 0.5 * (y_e[1:] + y_e[:-1])),
        data,
        np.vstack([x, y]).T,
        method="splinef2d",
        bounds_error=False,
    )
    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    if showid:
        mini = min(x.min(), y.min())
        maxi = max(x.max(), y.max())
        ax.plot([mini, maxi], [mini, maxi], color="red", linewidth=2)
    if showparab:
        T = np.linspace(x.min(), x.max(), 100)
        ax.plot(T, (T**2) / 4, color="red", linewidth=2)

    ax.scatter(x, y, c=z, s=size, **kwargs)
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax)
    cbar.ax.set_ylabel("Density")
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    if equalaxis:
        ax.axis("equal")
    ax.grid(True, which="both")
    ax.axhline(y=0, color="k")
    ax.axvline(x=0, color="k")
    if save:
        plt.savefig(save_name)
    if showfig:
        plt.show()
    else:
        plt.close()
