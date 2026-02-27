#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:06:53 2020

@author: Anna SONG

Convolutional filters PyTorch + GPU with buffered memory :
    - gaussian blur
    - differential operators: grad, Hess, Div
        (still without division by delta_x or delta_x^2)


Note: convention that u has shape (Z,X,Y)

Source: https://github.com/annasongmaths/curvatubes
License: MIT — Copyright (c) 2021 Anna SONG
See LICENSE in this directory.

NOTE (polnet integration): Modifications from the upstream source:

1. **Code formatting** — reformatted with ``black`` (line-length 79) on
   2026-02-26 for consistency with the rest of the polnet codebase.
   No functional changes.
"""

import torch

dtype = (
    torch.cuda.FloatTensor if torch.cuda.is_available() else torch.FloatTensor
)

# new syntax
import torch.fft

""" 3D Gaussian blur """


def separable_fourier_filter(N, sigma):
    "Returns the Fourier transform of a gaussian filter k of size sigma (vector of length 2N)"
    k = (torch.arange(-N, N).type(dtype) ** 2) / (2 * sigma**2)
    k = torch.cat((k[N:], k[:N]))  # Put 0 in the first cell ("fftshift")
    k = (-k).exp().view(1, 1, -1)  # (1,1,N)
    k = k / k.sum()

    fk_c = torch.fft.fft(k, dim=-1)  # (1,1,2*N) complex
    return fk_c


def gaussian_blur_3D_periodic(u, sigma):
    "Applies a Gaussian blur of size sigma on u, using Fourier transforms"
    if sigma == None:
        return u

    # input u is real
    N = u.shape[-1] // 2
    fk_c = separable_fourier_filter(N, sigma)  # (1,1,2*N) complex

    fu_c = torch.fft.fft(u, dim=-1)  # 1D Fourier Transform, (N,N,2*N) complex
    fu_c = fu_c * fk_c  # (N,N,2*N) complex
    u_c = torch.fft.ifft(fu_c, dim=-1).permute(
        0, 2, 1
    )  # (N,N,2*N) complex (with zero imag part)

    fu_c = torch.fft.fft(
        u_c, dim=-1
    )  # 1D Fourier Transform, (N,N,2*N) complex
    fu_c = fu_c * fk_c  # (N,N,2*N) complex
    u_c = torch.fft.ifft(fu_c, dim=-1).permute(
        2, 1, 0
    )  # (N,N,2*N) complex (with zero imag part)

    fu_c = torch.fft.fft(
        u_c, dim=-1
    )  # 1D Fourier Transform, (N,N,2*N) complex
    fu_c = fu_c * fk_c  # (N,N,2*N) complex
    u_c = torch.fft.ifft(fu_c, dim=-1).permute(
        2, 0, 1
    )  # (N,N,2*N) complex (with zero imag part)

    return u_c.real  # output is real


def general_gaussian_blur_3D_periodic(u, fk_c, fi_c, fj_c):
    "Applies a separable convolution whose Fourier filters are (fk,fi,fj) on u, typically anisotropic Gaussian blur"

    # input u is real
    fu_c = torch.fft.fft(u, dim=-1)  # 1D Fourier Transform, (Z,X,Y) complex
    fu_c = fu_c * fj_c  # (Z,X,Y) complex
    u_c = torch.fft.ifft(fu_c, dim=-1).permute(
        0, 2, 1
    )  # (Z,X,Y) complex (with zero imag part)

    fu_c = torch.fft.fft(u_c, dim=-1)  # 1D Fourier Transform, (Z,X,Y) complex
    fu_c = fu_c * fi_c  # (Z,X,Y) complex
    u_c = torch.fft.ifft(fu_c, dim=-1).permute(
        2, 1, 0
    )  # (Z,X,Y) complex (with zero imag part)

    fu_c = torch.fft.fft(u_c, dim=-1)  # 1D Fourier Transform, (Z,X,Y) complex
    fu_c = fu_c * fk_c  # (Z,X,Y) complex
    u_c = torch.fft.ifft(fu_c, dim=-1).permute(
        2, 0, 1
    )  # (Z,X,Y) complex (with zero imag part)

    return u_c.real  # output is real


class GeneralGaussianBlur3D_periodic(torch.nn.Module):
    "Applies a Gaussian blur of sizes (sig_z,sig_x,sig_y) in PERIODIC mode"

    def __init__(self, Z, X, Y, sig_z, sig_x, sig_y):
        super(GeneralGaussianBlur3D_periodic, self).__init__()

        if None in [sig_z, sig_x, sig_y]:
            self.trivial = True
        else:
            self.trivial = False
            self.fk_c = separable_fourier_filter(
                Z // 2, sig_z
            )  # (1,1,Z) complex
            self.fi_c = separable_fourier_filter(
                X // 2, sig_x
            )  # (1,1,X) complex
            self.fj_c = separable_fourier_filter(
                Y // 2, sig_y
            )  # (1,1,Y) complex

    def forward(self, u):

        if self.trivial:
            return u

        return general_gaussian_blur_3D_periodic(
            u, self.fk_c, self.fi_c, self.fj_c
        )  # (Z,X,Y) real


class GeneralGaussianBlur3D_notperiodic(torch.nn.Module):
    """Mimics a Gaussian blur of sizes (sig_z,sig_x,sig_y) in REPLICATE mode:
    - replicate padding at the borders with size ~ 3*(sig_z,sig_x,sig_y) pixels
    - periodic blurring on the extended array
    - crop down to initial size
    """

    def __init__(self, Z, X, Y, sig_z, sig_x, sig_y, padding_mode="replicate"):
        # Note: padding_mode can also be 'constant' (zero padding) or 'periodic'
        # (but then would be very close to the usual GeneralGaussianBlur3D_periodic).
        super(GeneralGaussianBlur3D_notperiodic, self).__init__()

        if None not in [sig_z, sig_x, sig_y]:

            self.trivial = False
            self.pad_z = int(3 * sig_z)
            self.pad_x = int(3 * sig_x)
            self.pad_y = int(3 * sig_y)
            self.padding_mode = padding_mode

            self.fk_pad_c = separable_fourier_filter(
                (Z + 2 * self.pad_z) // 2, sig_z
            )  # (1,1,Z + 2 pad_z) ? complex
            self.fi_pad_c = separable_fourier_filter(
                (X + 2 * self.pad_x) // 2, sig_x
            )  # (1,1,X + 2 pad_x) ? complex
            self.fj_pad_c = separable_fourier_filter(
                (Y + 2 * self.pad_y) // 2, sig_y
            )  # (1,1,Y + 2 pad_y) ? complex

        else:
            self.trivial = True

    def forward(self, u):

        if self.trivial:
            return u

        u_pad = torch.nn.functional.pad(
            u[None, None],
            (
                self.pad_y,
                self.pad_y,
                self.pad_x,
                self.pad_x,
                self.pad_z,
                self.pad_z,
            ),
            mode=self.padding_mode,
        )[0, 0]

        u_pad = general_gaussian_blur_3D_periodic(
            u_pad, self.fk_pad_c, self.fi_pad_c, self.fj_pad_c
        )
        u = u_pad[
            self.pad_z : -self.pad_z,
            self.pad_x : -self.pad_x,
            self.pad_y : -self.pad_y,
        ]

        return u


""" Differential operators: grad, Hess, Div (with buffered memory) """

"First, the padding modes"


def periodic_padding(u_pad, u):
    "Periodic padding of u (pad = 1 pix), with in-place copy of u into the buffer u_pad."

    c, C = 1, -1

    # Central chunk:
    u_pad[c:C, c:C, c:C] = u

    # The 6 2D panels:
    u_pad[0, c:C, c:C] = u[-1, :, :]
    u_pad[-1, c:C, c:C] = u[0, :, :]

    u_pad[c:C, 0, c:C] = u[:, -1, :]
    u_pad[c:C, -1, c:C] = u[:, 0, :]

    u_pad[c:C, c:C, 0] = u[:, :, -1]
    u_pad[c:C, c:C, -1] = u[:, :, 0]

    # The 12 1D lines:
    u_pad[0, 0, c:C] = u[-1, -1, :]
    u_pad[0, -1, c:C] = u[-1, 0, :]

    u_pad[-1, 0, c:C] = u[0, -1, :]
    u_pad[-1, -1, c:C] = u[0, 0, :]

    u_pad[0, c:C, 0] = u[-1, :, -1]
    u_pad[0, c:C, -1] = u[-1, :, 0]

    u_pad[-1, c:C, 0] = u[0, :, -1]
    u_pad[-1, c:C, -1] = u[0, :, 0]

    u_pad[c:C, 0, 0] = u[:, -1, -1]
    u_pad[c:C, 0, -1] = u[:, -1, 0]

    u_pad[c:C, -1, 0] = u[:, 0, -1]
    u_pad[c:C, -1, -1] = u[:, 0, 0]

    # The 8 corners: not useful
    # u_pad[0, 0, 0] = u[-1, -1, -1]

    # u_pad[-1, 0, 0] = u[0, -1, -1]
    # u_pad[0, -1, 0] = u[-1, 0, -1]
    # u_pad[0, 0, -1] = u[-1, -1, 0]

    # u_pad[-1, -1, 0] = u[0, 0, -1]
    # u_pad[0, -1, -1] = u[-1, 0, 0]
    # u_pad[-1, 0, -1] = u[0, -1, 0]

    # u_pad[-1, -1, -1] = u[0, 0, 0]


def replicate_padding(u_pad, u):
    "Replicate padding of u (pad = 1 pix), with in-place copy of u into the buffer u_pad."

    c, C = 1, -1

    # Central chunk:
    u_pad[c:C, c:C, c:C] = u

    # The 6 2D panels:
    u_pad[0, c:C, c:C] = u[0, :, :]
    u_pad[-1, c:C, c:C] = u[-1, :, :]

    u_pad[c:C, 0, c:C] = u[:, 0, :]
    u_pad[c:C, -1, c:C] = u[:, -1, :]

    u_pad[c:C, c:C, 0] = u[:, :, 0]
    u_pad[c:C, c:C, -1] = u[:, :, -1]

    # The 12 1D lines:
    u_pad[0, 0, c:C] = u[0, 0, :]
    u_pad[0, -1, c:C] = u[0, -1, :]

    u_pad[-1, 0, c:C] = u[-1, 0, :]
    u_pad[-1, -1, c:C] = u[-1, -1, :]

    u_pad[0, c:C, 0] = u[0, :, 0]
    u_pad[0, c:C, -1] = u[0, :, -1]

    u_pad[-1, c:C, 0] = u[-1, :, 0]
    u_pad[-1, c:C, -1] = u[-1, :, -1]

    u_pad[c:C, 0, 0] = u[:, 0, 0]
    u_pad[c:C, 0, -1] = u[:, 0, -1]

    u_pad[c:C, -1, 0] = u[:, -1, 0]
    u_pad[c:C, -1, -1] = u[:, -1, -1]

    # The 8 corners: not useful
    # u_pad[0, 0, 0] = u[0, 0, 0]

    # u_pad[-1, 0, 0] = u[-1, 0, 0]
    # u_pad[0, -1, 0] = u[0, -1, 0]
    # u_pad[0, 0, -1] = u[0, 0, -1]

    # u_pad[-1, -1, 0] = u[-1, -1, 0]
    # u_pad[0, -1, -1] = u[0, -1, -1]
    # u_pad[-1, 0, -1] = u[-1, 0, -1]


"Then, grad and Hess"


def grad_hessian(u_pad, grad, H_diag, H_off):
    "In-place filling of grad, H_diag (diag of Hess), H_off (off-diag of Hess) with u_pad"

    # Offset slices:
    m, M = 0, -2  # 0:-2
    c, C = 1, -1  # 1:-1
    p, P = 2, None  # 2:0

    # Center pixel:
    U = u_pad[c:C, c:C, c:C]

    # 6 neighbors:
    U_z, U_Z = u_pad[m:M, c:C, c:C], u_pad[p:, c:C, c:C]
    U_x, U_X = u_pad[c:C, m:M, c:C], u_pad[c:C, p:, c:C]
    U_y, U_Y = u_pad[c:C, c:C, m:M], u_pad[c:C, c:C, p:]

    # 12 2-neighbors:
    U_zx, U_zX = u_pad[m:M, m:M, c:C], u_pad[m:M, p:, c:C]
    U_Zx, U_ZX = u_pad[p:, m:M, c:C], u_pad[p:, p:, c:C]

    U_zy, U_zY = u_pad[m:M, c:C, m:M], u_pad[m:M, c:C, p:]
    U_Zy, U_ZY = u_pad[p:, c:C, m:M], u_pad[p:, c:C, p:]

    U_xy, U_xY = u_pad[c:C, m:M, m:M], u_pad[c:C, m:M, p:]
    U_Xy, U_XY = u_pad[c:C, p:, m:M], u_pad[c:C, p:, p:]

    # Gradient: [-1, 0, 1]/2
    grad[:, :, :, 0] = (U_Z - U_z) / 2  # d/dz
    grad[:, :, :, 1] = (U_X - U_x) / 2  # d/dx
    grad[:, :, :, 2] = (U_Y - U_y) / 2  # d/dy

    # Diagonal of the Hessian: [1, -2, 1]
    H_diag[:, :, :, 0] = U_Z + U_z - 2 * U  # d^2/dz^2
    H_diag[:, :, :, 1] = U_X + U_x - 2 * U  # d^2/dx^2
    H_diag[:, :, :, 2] = U_Y + U_y - 2 * U  # d^2/dy^2

    # Off-diagonal coefficients of the Hessian: [[-1, 0, 1], [0, 0, 0], [1, 0, -1]] / 4
    H_off[:, :, :, 0] = (U_ZX + U_zx - U_zX - U_Zx) / 4  # d^2/(dz*dx)
    H_off[:, :, :, 1] = (U_ZY + U_zy - U_zY - U_Zy) / 4  # d^2/(dz*dy)
    H_off[:, :, :, 2] = (U_XY + U_xy - U_xY - U_Xy) / 4  # d^2/(dx*dy)


def backward_grad(d_pad, d_u):
    "Manually encodes backward of the gradient, in-place filling of d_u with d_pad"

    # Offset slices:
    m, M = 0, -2  # 0:-2
    c, C = 1, -1  # 1:-1
    p, P = 2, None  # 2:0

    # Center pixel:
    # dG    = d_pad[c:C,c:C,c:C]

    # 6 neighbors:
    dG_z, dG_Z = d_pad[m:M, c:C, c:C, 0], d_pad[p:, c:C, c:C, 0]
    dG_x, dG_X = d_pad[c:C, m:M, c:C, 1], d_pad[c:C, p:, c:C, 1]
    dG_y, dG_Y = d_pad[c:C, c:C, m:M, 2], d_pad[c:C, c:C, p:, 2]

    # Convolution with [1, 0, -1]
    d_u[:, :, :] += (dG_z - dG_Z) / 2
    d_u[:, :, :] += (dG_x - dG_X) / 2
    d_u[:, :, :] += (dG_y - dG_Y) / 2


def backward_H_diag(d_pad, d_u):
    "Manually encodes backward of Hess diag, in-place filling of d_u with d_pad"

    # Offset slices
    m, M = 0, -2  # 0:-2
    c, C = 1, -1  # 1:-1
    p, P = 2, None  # 2:0

    # Center pixel:
    dHd = d_pad[c:C, c:C, c:C]

    # 6 neighbors:
    dHd_z, dHd_Z = d_pad[m:M, c:C, c:C, 0], d_pad[p:, c:C, c:C, 0]
    dHd_x, dHd_X = d_pad[c:C, m:M, c:C, 1], d_pad[c:C, p:, c:C, 1]
    dHd_y, dHd_Y = d_pad[c:C, c:C, m:M, 2], d_pad[c:C, c:C, p:, 2]

    # Convolution with [1, -2, 1]:
    d_u[:, :, :] += dHd_z + dHd_Z
    d_u[:, :, :] += dHd_x + dHd_X
    d_u[:, :, :] += dHd_y + dHd_Y
    d_u[:, :, :] -= 2 * dHd.sum(-1)


def backward_H_off(d_pad, d_u):
    "Manually encodes backward of Hess off-diag, in-place filling of d_u with d_pad"

    # Offset slices:
    m, M = 0, -2  # 0:-2
    c, C = 1, -1  # 1:-1
    p, P = 2, None  # 2:0

    # 12 2-neighbors:
    dHo_zx, dHo_zX = d_pad[m:M, m:M, c:C, 0], d_pad[m:M, p:, c:C, 0]
    dHo_Zx, dHo_ZX = d_pad[p:, m:M, c:C, 0], d_pad[p:, p:, c:C, 0]

    dHo_zy, dHo_zY = d_pad[m:M, c:C, m:M, 1], d_pad[m:M, c:C, p:, 1]
    dHo_Zy, dHo_ZY = d_pad[p:, c:C, m:M, 1], d_pad[p:, c:C, p:, 1]

    dHo_xy, dHo_xY = d_pad[c:C, m:M, m:M, 2], d_pad[c:C, m:M, p:, 2]
    dHo_Xy, dHo_XY = d_pad[c:C, p:, m:M, 2], d_pad[c:C, p:, p:, 2]

    # Convolution with [[-1, 0, 1], [0, 0, 0], [1, 0, -1]] / 4
    d_u[:, :, :] += (dHo_ZX + dHo_zx - dHo_zX - dHo_Zx) / 4  # d^2/(dz*dx)
    d_u[:, :, :] += (dHo_ZY + dHo_zy - dHo_zY - dHo_Zy) / 4  # d^2/(dz*dy)
    d_u[:, :, :] += (dHo_XY + dHo_xy - dHo_xY - dHo_Xy) / 4  # d^2/(dx*dy)


def my_custom_GradHess(mode):

    if mode == "periodic":
        padding = periodic_padding
    elif mode == "replicate":
        padding = replicate_padding

    class GradHessianFunc(torch.autograd.Function):

        @staticmethod
        def forward(ctx, u, u_pad, grad, H_diag, H_off, d_pad, d_u):

            ctx.save_for_backward(d_pad, d_u)

            padding(u_pad, u)  # Pad u
            grad_hessian(
                u_pad, grad, H_diag, H_off
            )  # Compute the gradient and Hessian in place

            return grad, H_diag, H_off

        @staticmethod
        def backward(ctx, d_grad, d_H_diag, d_H_off):

            d_pad, d_u = ctx.saved_tensors
            d_u.zero_()  # The backward_grad/H_diag/H_off act by increments

            # Backward from d_grad:
            padding(d_pad, d_grad)  # d_pad <- d_grad
            backward_grad(d_pad, d_u)  # d_pad -> d_u

            # Backward from d_H_diag:
            padding(d_pad, d_H_diag)  # d_pad <- d_H_diag
            backward_H_diag(d_pad, d_u)  # d_pad -> d_u

            # Backward from d_H_off:
            padding(d_pad, d_H_off)  # d_pad <- d_H_off
            backward_H_off(d_pad, d_u)  # d_pad -> d_u

            return d_u, None, None, None, None, None, None

    class GradHessian(torch.nn.Module):
        def __init__(self, Z, X, Y):
            super(GradHessian, self).__init__()

            # Pre-allocate all the memory
            self.register_buffer("u_pad", torch.zeros(Z + 2, X + 2, Y + 2))
            self.register_buffer("grad", torch.zeros(Z, X, Y, 3))
            self.register_buffer("H_diag", torch.zeros(Z, X, Y, 3))
            self.register_buffer("H_off", torch.zeros(Z, X, Y, 3))

            self.register_buffer("d_pad", torch.zeros(Z + 2, X + 2, Y + 2, 3))
            self.register_buffer("d_u", torch.zeros(Z, X, Y))

        def forward(self, u):

            return GradHessianFunc.apply(
                u,
                self.u_pad,
                self.grad,
                self.H_diag,
                self.H_off,
                self.d_pad,
                self.d_u,
            )

    return GradHessian


"Finally, divergence"


def fill_divergence(A_pad, div):

    # Offset slices:
    m, M = 0, -2  # 0:-2
    c, C = 1, -1  # 1:-1
    p, P = 2, None  # 2:0

    Az_pad, Ax_pad, Ay_pad = A_pad

    # 6 neighbors:
    Az_z, Az_Z = Az_pad[m:M, c:C, c:C], Az_pad[p:, c:C, c:C]
    Ax_x, Ax_X = Ax_pad[c:C, m:M, c:C], Ax_pad[c:C, p:, c:C]
    Ay_y, Ay_Y = Ay_pad[c:C, c:C, m:M], Ay_pad[c:C, c:C, p:]

    # Divergence : [-1, 0, 1] / 2
    div[...] = (Az_Z - Az_z) / 2 + (Ax_X - Ax_x) / 2 + (Ay_Y - Ay_y) / 2


def backward_divergence(d_pad, d_A):

    # Offset slices:
    m, M = 0, -2  # 0:-2
    c, C = 1, -1  # 1:-1
    p, P = 2, None  # 2:0

    # 6 neighbors:
    dD_z, dD_Z = d_pad[m:M, c:C, c:C], d_pad[p:, c:C, c:C]
    dD_x, dD_X = d_pad[c:C, m:M, c:C], d_pad[c:C, p:, c:C]
    dD_y, dD_Y = d_pad[c:C, c:C, m:M], d_pad[c:C, c:C, p:]

    # Convolution with [1, 0, -1] / 2
    d_A[0, :, :, :] += (dD_z - dD_Z) / 2
    d_A[1, :, :, :] += (dD_x - dD_X) / 2
    d_A[2, :, :, :] += (dD_y - dD_Y) / 2


def my_custom_Div(mode):

    if mode == "periodic":
        padding = periodic_padding
    elif mode == "replicate":
        padding = replicate_padding

    def A_padding(A_pad, A):
        padding(A_pad[0], A[0])
        padding(A_pad[1], A[1])
        padding(A_pad[2], A[2])

    class DivergenceFunc(torch.autograd.Function):

        @staticmethod
        def forward(ctx, A, A_pad, div, d_pad, d_A):

            ctx.save_for_backward(d_pad, d_A)

            A_padding(A_pad, A)  # Pad A
            fill_divergence(A_pad, div)  # Compute the divergence in place

            return div

        @staticmethod
        def backward(ctx, d_div):

            d_pad, d_A = ctx.saved_tensors
            d_A.zero_()  # The backward_divergence() acts by increments

            # Backward from d_div:
            padding(d_pad, d_div)  # d_pad <- d_div
            backward_divergence(d_pad, d_A)  # d_pad -> d_A

            return d_A, None, None, None, None

    class Divergence_class(torch.nn.Module):
        def __init__(self, Z, X, Y):
            super(Divergence_class, self).__init__()

            # Pre-allocate all the memory
            self.register_buffer("A_pad", torch.zeros(3, Z + 2, X + 2, Y + 2))
            self.register_buffer("div", torch.zeros(Z, X, Y))

            self.register_buffer("d_pad", torch.zeros(Z + 2, X + 2, Y + 2))
            self.register_buffer("d_A", torch.zeros(3, Z, X, Y))

        def forward(self, u):

            return DivergenceFunc.apply(
                u, self.A_pad, self.div, self.d_pad, self.d_A
            )

    return Divergence_class
