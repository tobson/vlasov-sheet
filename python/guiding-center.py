# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from DSolveGC import DSolveGC
import helpers
pylab = helpers.setup ('ps')
from matplotlib import rcParams, pyplot as plt
import os, sys

def left (axis):

    psi = np.pi/4
    bz = np.cos (psi)

    nkx = 8*1024
    kmax = np.sqrt (2*qshear)*Omega/va
    kz = .5*np.sqrt ((4-qshear)*qshear)*Omega/(va*bz)
    kx = np.sqrt (kmax**2/bz**2 - kz**2)*(np.arange (nkx) + .5)/nkx

    axis.text (.08, .38, labels[-1], fontsize = 8)
    axis.text (.08, .63, labels[0], fontsize = 8)
    axis.text (.08, .86, labels[1], fontsize = 8)
    axis.text (.08, 1.07, labels[2], fontsize = 8)
    axis.text (.08, 1.22, labels[3], fontsize = 8)
    axis.text (.08, 1.42, labels[4], fontsize = 8)

    gamma = np.empty (nkx)
    for level in [-1, 0, 1, 2, 3, 4]:
        beta = bz*bz*float (10**level)
        dsolve = DSolveGC (beta, psi, va, Omega, qshear)
        gamma[-1] = dsolve (kx[-1], kz, .1)
        for i in range (nkx-1)[::-1]:
            gamma[i] = dsolve (kx[i], kz, gamma[i+1])
        axis.plot (kx/kz, gamma/Omega, 'k')

    axis.set_xlim (0., 1.5)
    axis.set_ylim (0., 1.5)
    axis.set_xlabel (r'$k_x/k_z$')
    axis.set_ylabel (r'$\gamma/\Omega$')
    axis.text (.75, .85, r'$k_\parallel v_a/\Omega = \sqrt{15/16}$',
            horizontalalignment = 'center', verticalalignment = 'center',
            transform = axis.transAxes)

def right (axis):

    psi = np.pi/4.
    bz = np.cos (psi)

    nkz = 8*1024
    kmax = np.sqrt (2.*qshear)*Omega/(va*bz)
    kx = 0.
    kz = kmax*(np.arange (nkz) + .5)/nkz

    axis.text (.78, .1, labels[-4], fontsize = 8)
    axis.text (1.2, .22, labels[-2], fontsize = 8)
    axis.text (.88, .78, labels[0], fontsize = 8)
    axis.text (.73, 1.03, labels[1], fontsize = 8)
    axis.text (.53, 1.31, labels[2], fontsize = 8)
    axis.text (.33, 1.54, labels[4], fontsize = 8)
    axis.text (.12, 1.77, labels[8], fontsize = 8)

    gamma = np.empty (nkz)
    for level in [-4, -2, 0, 1, 2, 3, 8]:
        beta = bz*bz*float (10**level)
        dsolve = DSolveGC (beta, psi, va, Omega, qshear)
        if beta > 1.:
            gamma[-1] = dsolve (kx, kz[-1], .1)
            for i in range (nkz-1)[::-1]:
                gamma[i] = dsolve (kx, kz[i], gamma[i+1])
        else:
            gamma[0] = dsolve (kx, kz[0], 1e-3)
            for i in range (1, nkz):
                gamma[i] = dsolve (kx, kz[i], gamma[i-1])
        axis.plot (kz*bz*va/Omega, gamma/Omega, 'k')

    axis.set_xlim (0., 2.)
    axis.set_ylim (0., 2.)
    axis.set_xlabel (r'$k_\parallel v_A/\Omega$')
    axis.text (.75, .85, r'$k_x=0$', horizontalalignment = 'center',
            verticalalignment = 'center', transform = axis.transAxes)

# Units
Omega = 1.
va = 1.

# Logarithmic rate of shear
qshear = 1.5

# Create figure
width, height = plt.rcParams['figure.figsize']
height = .4*width
fig = plt.figure (num = 1, figsize = (width, height))
fig.clf ()
fig, axes = plt.subplots (num = 1, ncols = 2)

# Line labels
labels = {}
for n in range (-10, 10):
    labels[n] = r'$10^{%g}$' % n
labels[0] = r'$1$'
labels[1] = r'$10$'

# Plots
left (axes[0])
right (axes[1])

# Save figure
if not pylab:
    tight = dict (pad = 0., w_pad = .5, h_pad = .5,
            rect = [6e-3,.005,1-1e-3,1-5e-3])
    fig.set_tight_layout (tight)
    helpers.savefig (fig, __file__)
