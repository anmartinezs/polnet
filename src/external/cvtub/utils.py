#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:03:21 2020

@author: Anna SONG

Auxiliary functions:
    - for file conversion: from nii.gz to np.array (and conversely)
    - for display: show slices of u
    - for initialization: u = u0 (random noise, balls, cubes)
        (not so useful if using the change of variable u = div A + m0 )
    - functions on u: double-well W(u), impose the average
    - obtain the coeffs (a20,a11,a02,b10,b01,c) from the non-reduced polynomial
        h2 (H-H0)^2 + k1 K + s (kap1 - kap1_0)^2 + t (kap2 - kap2_0)^2

Source: https://github.com/annasongmaths/curvatubes
License: MIT — Copyright (c) 2021 Anna SONG
See LICENSE in this directory.

NOTE (polnet integration): Modifications from the upstream source:

1. **Code formatting** — reformatted with ``black`` (line-length 79) on
   2026-02-26 for consistency with the rest of the polnet codebase.
   No functional changes.
"""

import numpy as np
import matplotlib.pyplot as plt
import re

import torch

dtype = (
    torch.cuda.FloatTensor
)  # if torch.cuda.is_available() else torch.FloatTensor

import nibabel as nib
from skimage.morphology import ball
from scipy.ndimage import median_filter

"Convert nii.gz object <-> 3D np.array"


def save_nii(vol, filename, compressed=True):
    "Saves a 3D np.array into a .nii.gz or .nii file"
    img = nib.Nifti1Image(vol, np.eye(4))
    if compressed:
        img.to_filename(filename + ".nii.gz")
    else:
        img.to_filename(filename + ".nii")


def load_nii(nii_file):
    "Loads a .nii.gz or .nii file into a np.array"
    img = nib.load(nii_file)
    vol = np.asarray(img.get_data())
    return vol


"""Display slices of a volume"""


def slices(
    u, figsize=(12, 4), rescale=True, cmap="gray", save=False, title=""
):
    """Visualize three 2D slices of a 3D volume at z = Z//3, Z//2, or 2*Z//3
    rescale = True: grays between u.min and u.max
    rescale = False: grays between 0 and 1 (black/white beyond 0/1)"""
    if type(u) == torch.Tensor:
        u = u.detach().cpu().numpy()
    vmin = None
    vmax = None
    Z = u.shape[0]
    if not rescale:
        vmin = 0.0
        vmax = 1.0
    fig, ax = plt.subplots(1, 3, figsize=figsize, sharex=True, sharey=True)
    ax[0].imshow(
        np.asarray(u[Z // 3], dtype=np.float64),
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )
    ax[1].imshow(
        np.asarray(u[Z // 2], dtype=np.float64),
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )
    ax[2].imshow(
        np.asarray(u[2 * Z // 3], dtype=np.float64),
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )
    fig.tight_layout()
    if title != "":
        fig.suptitle(title)
        fig.savefig(title + ".png")
    plt.show()


def single(img, figsize=(10, 10), rescale=True):
    "Visualize a single 2D image"
    if type(img) == torch.Tensor:
        img = img.detach().cpu().numpy()
    vmin = None
    vmax = None
    if not rescale:
        vmin = 0.0
        vmax = 1.0
    plt.figure(figsize=figsize)
    plt.imshow(img, vmin=vmin, vmax=vmax, cmap="gray")
    plt.show()


def rescale(vol):
    "rescales between 0 and 1"
    vol = vol - vol.min()
    vol = vol / (vol.max())
    return vol


"""Functions for Initialization of the phase-field u"""


def random_init(size, prop=0.2):  # prop between 0 and 1
    """Initializes a 3D torch tensor with ratios prop of +1 pixels
    and (1-prop) of -1 pixels."""
    if isinstance(size, int):
        N = size
        v0 = torch.rand((N, N, N)).type(dtype)
    else:
        v0 = torch.rand(size).type(dtype)
    v0 = 1.0 * (v0 <= prop)
    return 2 * v0 - 1


def smallest_rem(K, n):
    "Used in place_balls() for periodicity. Returns K modulo n but close to zero."
    a = np.mod(K, n)
    return a * (a <= n // 2) + (a - n) * (a > n // 2)


def give_ZXY(size):
    "Retrieves the 3D size from an integer or a triplet"
    if isinstance(size, int):
        Z, X, Y = size, size, size
    else:
        Z, X, Y = size
    return Z, X, Y


def init_holed_cubes(
    size,
    N_cubes=1,
    option="",
    mode="replicate",
    R=20,
    r=10,
    med=False,
    med_rad=0,
):
    "Initializes u = holed cubes"
    Z, X, Y = give_ZXY(size)
    cubes = np.zeros((Z, X, Y))

    if N_cubes == 1:
        cubes[
            Z // 2 - R : Z // 2 + R,
            X // 2 - R : X // 2 + R,
            Y // 2 - R : Y // 2 + R,
        ] = 1
        cubes[:, X // 2 - r : X // 2 + r, Y // 2 - r : Y // 2 + r] = 0

    if N_cubes == 2:

        if option in ["", "very_close", "close", "apart"]:
            if option in ["very_close"]:
                a = 2 * Y // 5
                b = 3 * Y // 5
            if option in ["", "close"]:
                a = Y // 3
                b = 2 * Y // 3
            elif option == "apart":
                a = Y // 4
                b = 3 * Y // 4
            cubes[
                Z // 2 - R : Z // 2 + R, X // 2 - R : X // 2 + R, a - R : b + R
            ] = 1
            cubes[:, X // 2 - r : X // 2 + r, a - r : a + r] = 0
            cubes[:, X // 2 - r : X // 2 + r, b - r : b + r] = 0

        if option in [
            "Aa_very_close",
            "Aa_close",
            "Aa_apart",
        ]:  # not the same direction for the holes
            if option == "Aa_very_close":
                a = 2 * Y // 5
                b = 3 * Y // 5
            if option == "Aa_close":
                a = Y // 3
                b = 2 * Y // 3
            if option == "Aa_apart":
                a = Y // 4
                b = 3 * Y // 4
            cubes[
                Z // 2 - R : Z // 2 + R, X // 2 - R : X // 2 + R, a - R : b + R
            ] = 1
            cubes[:, X // 2 - r : X // 2 + r, a - r : a + r] = 0
            cubes[Z // 2 - r : Z // 2 + r, :, b - r : b + r] = 0

    if med == True:
        cubes = median_filter(cubes, footprint=ball(med_rad))

    cubes = 2.0 * cubes - 1
    v0 = torch.Tensor(cubes).type(dtype)
    return v0


def place_balls(
    Z, X, Y, ZZ, XX, YY, coords, r_pix, mode, dist, return_tensor=True
):
    """Used in init_balls(). Places balls with centers coords and radius r_pix"""

    def distance(z, x, y, dist):
        if dist == "L2":
            return np.sqrt(x**2 + y**2 + z**2)
        if dist == "L1":
            return np.abs(x) + np.abs(y) + np.abs(z)

    balls = np.zeros((Z, X, Y))
    if mode == "periodic":
        for j in range(coords.shape[1]):

            z = smallest_rem(ZZ - coords[0, j], Z)
            x = smallest_rem(XX - coords[1, j], X)
            y = smallest_rem(YY - coords[2, j], Y)
            balls += distance(z, x, y, dist) <= r_pix

    elif mode == "replicate":
        for j in range(coords.shape[1]):
            z = ZZ - coords[0, j]
            x = XX - coords[1, j]
            y = YY - coords[2, j]
            balls += distance(z, x, y, dist) <= r_pix

    balls = 2.0 * (balls != 0) - 1
    if return_tensor:
        balls = torch.Tensor(balls).type(dtype)
    return balls


def init_balls(
    size, N_balls=10, r_pix=10, dist_bord=0, mode="replicate", dist="L2"
):
    """Initializes random balls of radius r_pix (pixels) inside a cubic domain
    with mode = periodic or mode = replicate.
    dist_bord specifies a min distance of the centers to the borders"""

    Z, X, Y = give_ZXY(size)
    XX, ZZ, YY = np.meshgrid(range(X), range(Z), range(Y))

    lowZ, highZ = dist_bord, Z - dist_bord
    lowX, highX = dist_bord, X - dist_bord
    lowY, highY = dist_bord, Y - dist_bord
    coords = np.concatenate(
        (
            np.random.randint(lowZ, highZ, size=N_balls)[None],
            np.random.randint(lowX, highX, size=N_balls)[None],
            np.random.randint(lowY, highY, size=N_balls)[None],
        ),
        axis=0,
    )

    return place_balls(Z, X, Y, ZZ, XX, YY, coords, r_pix, mode, dist)


def coords_n_balls(size, n, option=""):
    """Used in init_custom_balls().
    Gives the centers of n balls in special configurations."""
    Z, X, Y = give_ZXY(size)
    if n == 1:
        return np.array([[Z // 2], [X // 2], [Y // 2]])
    if n == 2:
        if option == "apart":
            return np.array(
                [[Z // 2, Z // 2], [X // 2, X // 2], [Y // 4, 3 * Y // 4]]
            )
        if option == "close":
            return np.array(
                [[Z // 2, Z // 2], [X // 2, X // 2], [Y // 3, 2 * Y // 3]]
            )
    if n == 3:
        if option == "right":
            return np.array(
                [
                    [Z // 2, Z // 2, Z // 2],
                    [X // 4, X // 4, 3 * X // 4],
                    [Y // 4, 3 * Y // 4, Y // 4],
                ]
            )
    if n == 4:
        if option == "apart":
            return np.array(
                [
                    [Z // 2, Z // 2, Z // 2, Z // 2],
                    [X // 4, X // 4, 3 * X // 4, 3 * X // 4],
                    [Y // 4, 3 * Y // 4, Y // 4, 3 * Y // 4],
                ]
            )
        if option == "apart_up":
            return np.array(
                [
                    [Z // 3, Z // 3, Z // 3, Z // 3],
                    [X // 4, X // 4, 3 * X // 4, 3 * X // 4],
                    [Y // 4, 3 * Y // 4, Y // 4, 3 * Y // 4],
                ]
            )
        if option == "close":
            return np.array(
                [
                    [Z // 2, Z // 2, Z // 2, Z // 2],
                    [X // 3, X // 3, 2 * X // 3, 2 * X // 3],
                    [Y // 3, 2 * Y // 3, Y // 3, 2 * Y // 3],
                ]
            )


def init_custom_balls(
    size, N_balls=1, option="", mode="replicate", r_pix=20, dist="L2"
):
    "Initialize balls in a few special configurations (see coords_n_balls() )"

    Z, X, Y = give_ZXY(size)
    XX, ZZ, YY = np.meshgrid(range(X), range(Z), range(Y))

    coords = coords_n_balls(size, N_balls, option=option)

    return place_balls(Z, X, Y, ZZ, XX, YY, coords, r_pix, mode, dist)


def init_tubes_voronoi(size, N_pts, N_kmeans, rad_pix, rad_med=6):
    "Initialize random tubes along Voronoi segments."
    from scipy.spatial import Voronoi
    import skfmm

    Z, X, Y = give_ZXY(size)
    pts = np.concatenate(
        (
            np.random.randint(Z, size=N_pts)[None],
            np.random.randint(X, size=N_pts)[None],
            np.random.randint(Y, size=N_pts)[None],
        ),
        axis=0,
    )
    pts = pts.T  # shape (N,3)

    from sklearn.cluster import KMeans

    kmeans = KMeans(n_clusters=N_kmeans)
    kmeans.fit(pts)
    centers = kmeans.cluster_centers_

    def check_pt_in_domain(pt):
        return 0 <= pt[0] < Z and 0 <= pt[1] < X and 0 <= pt[2] < Y

    vor = Voronoi(centers)

    binvol = np.zeros((Z, X, Y))
    for ind_faces in vor.ridge_vertices:  # in 3D the "edges" are faces
        if min(ind_faces) >= 0:  # index -1 means a point at infinity
            for i in range(len(ind_faces) - 1):
                x = vor.vertices[ind_faces[i]]
                y = vor.vertices[ind_faces[i + 1]]

                if check_pt_in_domain(x) and check_pt_in_domain(y):
                    # adding a segment
                    N = np.abs(x - y).max()  # greatest coordinate difference
                    t = np.linspace(0, 1, num=N)
                    coords = x[:, None] + t * (y[:, None] - x[:, None])
                    coords = np.round(coords)  # shape (3,...)
                    binvol[tuple(np.round(coords).astype(int))] = 1

    binvol = skfmm.distance(1 - binvol, dx=1) <= rad_pix
    if rad_med > 0:
        binvol = median_filter(binvol, footprint=ball(rad_med))
    vol = 2.0 * (binvol > 0) - 1

    return vol


"""Functions for phase-field u: impose the mean, double-well, softplus (...)^+ """


def project_average(u, m=0):
    "Translates u in order to have [mean of u on the domain] = m"
    return u - u.mean() + m


def double_W(s):
    "Double-well centered at +1 and -1"
    return 0.25 * (1 - s**2) ** 2


def double_W_prime(s):
    "Derivative of the double-well"
    return s**3 - s  # = s * (s**2 - 1)


def manual_softplus(u, xi_bis=1e-6):
    """Softplus approximation of the positive part u^+ with
    xi * log(1 + exp(u / xi)).
    Requires splitting wrt the sign of u for stable computations."""
    return (
        xi_bis
        + u * (u > 0)
        + xi_bis * torch.log(1 + torch.exp(-torch.abs(u) / xi_bis))
    )


""" Convert polynomial coefficients: non-reduced --> reduced """


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):  # human order : exp 100 before exp 1000
    return [atoi(c) for c in re.split(r"(\d+)", text)]


def extract_exp_nb(text):  # text = 'exp 100 ......'
    return [atoi(c) for c in re.split(r"(\d+)", text)][1]


def convert_params(h2, H0, k1, s, kap1_0, t, kap2_0):
    """Retrieves the polynomial coefficients (a20,a11,a02,b10,b01,c) from
    the non-reduced polynomial
    h2 (H-H0)^2 + k1 K + s (kap1 - kap1_0)^2 + t (kap2 - kap2_0)^2"""
    a20 = h2 + s
    a11 = 2 * h2 + k1
    a02 = h2 + t
    b10 = -2 * h2 * H0 - 2 * s * kap1_0
    b01 = -2 * h2 * H0 - 2 * t * kap2_0
    c = h2 * H0**2 + s * kap1_0**2 + t * kap2_0**2
    coeffs = a20, a11, a02, b10, b01, c
    return np.array(coeffs)
