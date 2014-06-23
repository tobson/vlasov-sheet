from __future__ import division
import numpy as np
import os, sys, time
from DSolveParallel import DSolveParallel
from findmax import findmax
import helpers
pylab = helpers.setup ('ps')
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from progressbar import ProgressBar

def compute ():

    nkz = 1024
    noc = 1024
    nbeta = 1024

    oc_p = np.logspace (np.log10 (.2), 4, noc)
    oc_m = -np.logspace (np.log10 (2.), 4, noc)
    oc = np.vstack ([oc_p, oc_m])

    beta = np.logspace (-1, 5, nbeta)

    gamma_max = np.empty (oc.shape + beta.shape, dtype = np.double)
    kz_opt = np.empty (oc.shape + beta.shape, dtype = np.double)

    beta_start = 1e-2

    pbar = ProgressBar ().start ()

    for ioc, oc1 in np.ndenumerate (oc):

        if ioc == (1,0): continue

        dsolve = DSolveParallel (oc1)
        kz = dsolve.kmax*(np.arange (nkz) + .5)/nkz
        gamma = dsolve.cold (kz)
        logspace = np.logspace (np.log10 (beta_start), np.log10 (beta[0]), 128)

        for beta1 in logspace:
            try:
                gamma = dsolve.warm (beta1, kz, gamma)
            except RuntimeError as e:
                print e.message
                print "oc = %g, beta = %g" % (oc1, beta1)
                raise SystemExit

        for ibeta, beta1 in np.ndenumerate (beta):
            try:
                gamma = dsolve.warm (beta1, kz, gamma)
            except RuntimeError as e:
                print e.message
                print "oc = %g, beta = %g" % (oc1, beta1)
                raise SystemExit
            ind = ioc + ibeta
            kz_opt[ind], gamma_max[ind] = findmax (kz, gamma)

        percentage = pbar.maxval*ioc[1]/(noc - 1)
        pbar.update (percentage)

    return oc, beta, kz_opt, gamma_max

# Units
Omega = 1.
va = 1.

# Logarithmic rate of shear
qshear = 1.5

# Obtain data
data = 'vf-mri-2D'
try: 
    oc, beta, kz_opt, gamma_max = helpers.load (data)
except:
    print 'Recomputing...'
    oc, beta, kz_opt, gamma_max = compute ()
    helpers.save (data, oc, beta, kz_opt, gamma_max)

oc_p, oc_m = oc
gamma_p, gamma_m = gamma_max
kz_p, kz_m = kz_opt

gamma_m = np.ma.masked_where (gamma_m > 1.2, gamma_m)
kz_m = np.ma.masked_where (kz_m > 2., kz_m)

cmap = plt.cm.jet
cmap.set_over ('w', 1.)

fig = plt.figure (num = 1)
fig.clf ()

grid = ImageGrid (fig, [.065, .09, .89, .855],
        nrows_ncols = (2, 2), axes_pad = .35,
        cbar_location = 'right', cbar_mode = 'each', cbar_pad = .1)

extent_p = np.array ([0,4.5,0,2.49], dtype = np.double)
extent_m = np.array ([0,3.5,0,2.49], dtype = np.double)

imshow_kwds = {'cmap': cmap, 'origin': 'lower'}

data = [gamma_p.transpose (), gamma_m.transpose (),
        kz_p.transpose (), kz_m.transpose ()]
extent = [extent_p, extent_m, extent_p, extent_m]

from matplotlib.ticker import MultipleLocator
for i in range (4):
    im = grid[i].imshow (data[i], extent = extent[i], **imshow_kwds)
    colorbar = grid.cbar_axes[i].colorbar (im)
    colorbar.ax.yaxis.set_major_locator (MultipleLocator (.2))
colorbar.ax.yaxis.set_major_locator (MultipleLocator (.4))

grid[1].invert_xaxis ()

beta_ticks = np.array ([1e-1, 1e1, 1e3, 1e5])
yticklabels = [r'$0.1$', r'$10$', r'$10^3$', r'$10^5$']
logrange = np.log (beta_ticks/beta.min ())/np.log (beta.max ()/beta.min ())
yticks = extent_p[2] + (extent_p[3] - extent_p[2])*logrange
for ax in grid:
    ax.set_yticks (yticks)
    ax.set_yticklabels (yticklabels)

oc_ticks_p = np.array ([.2, 2., 20., 200., 2000.])
xticklabels = [r'$0.1$', r'$1$', r'$10$', r'$100$', r'$1000$']
logrange = np.log (oc_ticks_p/oc_p.min ())/np.log (oc_p.max ()/oc_p.min ())
xticks = extent_p[0] + (extent_p[1] - extent_p[0])*logrange
for ax in grid[0::2]:
    ax.set_xticks (xticks)
    ax.set_xticklabels (xticklabels)

oc_ticks_m = np.array ([-2000., -200., -20., -2.])
xticklabels = [r'$-1000$', r'$-100$', r'$-10$', r'$-1$']
logrange = np.log (oc_ticks_m/oc_m.max ())/np.log (oc_m.min ()/oc_m.max ())
xticks = extent_m[0] + (extent_m[1] - extent_m[0])*logrange
for ax in grid[1::2]:
    ax.set_xticks (xticks)
    ax.set_xticklabels (xticklabels)

ratio = np.log (oc_p.max ()/oc_p.min ())/np.log (beta.max ()/beta.min ())
ngrid = 128
log_oc = np.linspace (extent_p[0], extent_p[1], ngrid)
log_beta = ratio*np.linspace (extent_p[2], extent_p[3], ngrid)
grid[0].plot (log_oc[ngrid/2:], log_beta[ngrid/2:], 'k--')

ratio = np.log (oc_m.min ()/oc_m.max ())/np.log (beta.max ()/beta.min ())
ngrid = 128
log_oc = np.linspace (extent_m[0], extent_m[1], ngrid)
log_beta = ratio*np.linspace (extent_m[2], extent_m[3], ngrid)
grid[1].plot (log_oc[ngrid/2:], log_beta[ngrid/2:] + .22*extent_m[3], 'k--')

for ax in grid[:2]:
    ax.set_title (r'$\gamma/\Omega$')
for ax in grid[2:]:
    ax.set_title (r'$k_z v_a/\Omega$')
    ax.set_xlabel (r'$\omega_cb_z/2\Omega$')
for ax in grid[::2]:
    ax.set_ylabel (r'$\beta$', labelpad = -.2)

if not pylab:
    helpers.savefig (fig, __file__)
