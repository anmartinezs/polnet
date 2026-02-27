#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 21:49:22 2021

@author: Anna SONG

This script builds the phase-field energy Feps(u) used for the generation of 3D
tubular shapes inside the module polykap_deg2()

It also defines the Cahn-Hilliard energy, as well as the discrepancy and the
discrepancy normalized by Cahn-Hilliard.

Note: convention that u has shape (Z,X,Y)

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

import torch

dtype = (
    torch.cuda.FloatTensor
)  # if torch.cuda.is_available() else torch.FloatTensor

from .utils import (
    double_W,
    double_W_prime,
    manual_softplus,
)  # polnet: relative import (was ``from cvtub.utils ...``)

""" Energy (curvature functional) """


def polykap_deg2(u, params, delta_x, xi, GradHessConv_ZXY):
    """
    Phase-field Feps(u)

        Feps(u) = integral of (1/eps) (eps^2 |grad u|^2)  {
                    a[2,0] kap_{1,eps}^2 + a[1,1] kap_{1,eps} kap_2 + a[0,2] kap_{2,eps}^2
                    + b[1,0] kap_{1,eps} + b[0,1] kap_{2,eps} + c
                    } dx^3

    that is a diffuse approximation of the sharp-interface energy (see article)

        F(S)    = integral on S of {
                    a[2,0] kap_1^2 + a[1,1] kap_1 kap_2 + a[0,2] kap_2^2
                    + b[1,0] kap_1 + b[0,1] kap_2 + c
                    } dH^2

    with params = eps, a20, a11, a02, b10, b01, c, [mu, theta]    (mu & theta optional)
    each of these parameters can be spatialized, if you give torch.tensors.

    The diffuse principal curvatures kap_{1,eps} and kap_{2,eps} are such that

        2 eps |grad u| kap_{1,eps} = Tr Meps + sqrt{ (|Meps|^2 - (Tr Meps)^2 )^+ }
        2 eps |grad u| kap_{2,eps} = Tr Meps - sqrt{ (|Meps|^2 - (Tr Meps)^2 )^+ }

    where we use the matrix field

        Meps = - eps Hess u + W'(u) / eps
        # sign convention: trace is positive at {u = 0}, if u > 0 inside balls and < 0 outside

    """

    if len(params) == 9:
        eps, a20, a11, a02, b10, b01, c, mu, theta = params
        orientation = True
    elif len(params) == 7:
        eps, a20, a11, a02, b10, b01, c = params
        orientation = False

    grad_, Hess_diag_, Hess_off_ = GradHessConv_ZXY(u)

    grad = (
        grad_ / delta_x
    )  # we cannot perform in-place divisions, such as 'grad_ /= delta_x'
    Hess_diag = Hess_diag_ / delta_x**2
    Hess_off = Hess_off_ / delta_x**2

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
    Wprim = double_W_prime(u)

    # nHn is n^T H n where H is the hessian of u and n approximates the direction
    # of the gradient (not well defined when grad u ~ 0) up to a small offset xi
    nHn = (
        Hzz_u * dz_u**2
        + Hxx_u * dx_u**2
        + Hyy_u * dy_u**2
        + 2 * (Hzx_u * dz_u * dx_u + Hzy_u * dz_u * dy_u + Hxy_u * dx_u * dy_u)
    )
    nHn /= norm_grad_sq + xi**2

    Tra = -eps * (Hzz_u + Hxx_u + Hyy_u) + Wprim / eps  # Tr Meps

    Nor_pow2 = (
        (eps**2) * ((Hess_diag**2).sum(-1) + 2 * (Hess_off**2).sum(-1))
        + Wprim**2 / (eps**2)
        - 2 * Wprim * nHn
    )  # |Meps|^2 = the SQUARED norm of Meps

    positive_part = manual_softplus(
        2 * Nor_pow2 - Tra**2
    )  # with a small offset
    sqrt_part = torch.sqrt(positive_part)
    ngd_u = (norm_grad_sq + xi**2).sqrt()

    integ = 0

    integ += ((a20 + a02 - a11) / (2 * eps) * Nor_pow2).sum()
    integ += ((a20 - a02) / (2 * eps) * Tra * sqrt_part).sum()
    integ += (a11 / (2 * eps) * Tra**2).sum()

    integ += ((b10 + b01) / 2 * ngd_u * Tra).sum()
    integ += ((b10 - b01) / 2 * ngd_u * sqrt_part).sum()

    integ += (c * eps * norm_grad_sq).sum()

    if orientation and mu != 0 and theta is not None:  # add orientation energy
        sca_prod = dz_u * theta[0] + dx_u * theta[1] + dy_u * theta[2]
        integ += (mu * eps * sca_prod**2).sum()

    # integ += ... you can add anything you want in the energy

    integ *= delta_x**3

    return integ


""" Cahn-Hilliard and discrepancy """
from .filters import (
    my_custom_GradHess,
)  # polnet: relative import (was ``from cvtub.filters ...``)


def auxiliary_function(u, eps, delta_x, mode):
    if type(u) != torch.Tensor:
        u = torch.tensor(u).type(dtype)

    Z, X, Y = u.shape
    GradHessian = my_custom_GradHess(mode)
    if torch.cuda.is_available():
        GradHessConv_ZXY = GradHessian(Z, X, Y).cuda()
    else:
        GradHessConv_ZXY = GradHessian(Z, X, Y)

    grad_, Hess_diag_, Hess_off_ = GradHessConv_ZXY(u)
    del GradHessConv_ZXY

    grad = (
        grad_ / delta_x
    )  # we cannot perform in-place divisions, such as 'grad_ /= delta_x'
    Hess_diag = Hess_diag_ / delta_x**2
    Hess_off = Hess_off_ / delta_x**2

    return grad, Hess_diag, Hess_off


def CahnHilliard(u, eps=0.02, delta_x=0.01, mode="periodic"):
    """Cahn-Hilliard energy of u:
    integral of { eps / 2 |grad u|^2 + W(u) / eps } dx^3"""
    if type(u) != torch.Tensor:
        u = torch.tensor(u).type(dtype)

    grad, Hess_diag, Hess_off = auxiliary_function(u, eps, delta_x, mode)

    norm_grad_sq = (grad**2).sum(-1)

    CH = (eps / 2) * norm_grad_sq + double_W(u) / eps

    return torch.abs(CH).sum() * delta_x**3


def discrepancy(u, eps=0.02, delta_x=0.01, mode="periodic"):
    """Discrepancy of u:
    integral of { |eps / 2 |grad u|^2 - W(u) / eps| } dx^3"""
    if type(u) != torch.Tensor:
        u = torch.tensor(u).type(dtype)

    grad, Hess_diag, Hess_off = auxiliary_function(u, eps, delta_x, mode)

    norm_grad_sq = (grad**2).sum(-1)

    discr = (eps / 2) * norm_grad_sq - double_W(u) / eps

    return torch.abs(discr).sum() * delta_x**3


def ratio_discr(u, eps=0.02, delta_x=0.01, mode="periodic"):
    """Normalized discrepancy of u:
    integral of { |eps / 2 |grad u|^2 - W(u) / eps| } dx^3
                            /
    integral of { eps / 2 |grad u|^2 + W(u) / eps } dx^3"""

    if type(u) != torch.Tensor:
        u = torch.tensor(u).type(dtype)

    grad, Hess_diag, Hess_off = auxiliary_function(u, eps, delta_x, mode)

    norm_grad_sq = (grad**2).sum(-1)

    discr = (eps / 2) * norm_grad_sq - double_W(u) / eps
    CH = (eps / 2) * norm_grad_sq + double_W(u) / eps

    ratio = torch.abs(discr).sum() / CH.sum()

    return ratio
