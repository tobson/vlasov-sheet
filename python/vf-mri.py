from __future__ import division
import numpy as np
import os, sys, time
from DSolveParallel import DSolveParallel
from findmax import findmax
import helpers
pylab = helpers.setup ('ps')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from progressbar import ProgressBar

def compute ():

    # How many wave numbers?
    nkz = 1024

    # Range of cyclotron frequencies
    noc = 512
    oc_p = np.logspace (np.log10 (.2), np.log10 (2e4), noc)
    oc_m = -np.logspace (np.log10 (2.), np.log10 (2e4), noc+1)[1:]
    oc = np.vstack ([oc_p, oc_m])

    # Range of plasma beta's
    beta = np.array ([1e0, 1e1, 1e2, 1e3, 1e4])

    # Arrays to hold optimal wave numbers and maximum growth rates
    kz_opt = np.empty (oc.shape + beta.shape, dtype = np.double)
    gamma_max = np.empty (oc.shape + beta.shape, dtype = np.double)

    pbar = ProgressBar ().start ()

    # Outer loop over cyclotron frequencies
    for ioc, oc1 in np.ndenumerate (oc):
        # Create dispersion relation solver
        dsolve = DSolveParallel (oc1)
        # Get range of wave numbers
        kz = dsolve.kmax*(np.arange (nkz) + .5)/nkz
        # Cold plasma growth rate. This is used as initial guess
        gamma = dsolve.cold (kz)
        # Start sweep in beta at this value
        beta_min = 1e-2
        # Increment log (beta) by this much each iteration
        dlogbeta = .005
        # Coarse loop over beta
        for j, beta1 in np.ndenumerate (beta):
            ratio = beta1/beta_min
            nbeta = np.log10 (ratio)/dlogbeta
            # Fine grained loop over beta
            for beta2 in beta_min*pow (ratio, np.arange (nbeta + 1)/nbeta):
                try:
                    gamma = dsolve.warm (beta2, kz, gamma)
                except RuntimeError as e:
                    print e.message
                    print "oc = %g, beta = %g" % (oc1, beta2)
                    raise SystemExit
            # Save result
            kz_opt[ioc + j], gamma_max[ioc + j] = findmax (kz, gamma)
            # Start next fine grained loop at this value
            beta_min = beta1

        percentage = pbar.maxval*ioc[1]/(noc - 1)
        pbar.update (percentage)

    return oc, beta, kz_opt, gamma_max

# Units
Omega = 1.
va = 1.

# Logarithmic rate of shear
qshear = 1.5

# Obtain data
data = 'vf-mri'
try: 
    oc, beta, kz_opt, gamma_max = helpers.load (data)
except:
    print 'Recomputing...'
    oc, beta, kz_opt, gamma_max = compute ()
    helpers.save (data, oc, beta, kz_opt, gamma_max)

fig = plt.figure (num = 1)
fig.clf ()

width_ratios = [np.log10 (oc[0,:].max ()/oc[0,:].min ()),
                np.log10 (oc[1,:].min ()/oc[1,:].max ())]
gs = gridspec.GridSpec (2, 2, width_ratios = width_ratios)

ax = [fig.add_subplot (gs1) for gs1 in gs]

fmt = mpl.ticker.LogFormatterMathtext ()

colors = ['b','g','r','c','m']

for i, (beta1, color) in enumerate (zip (beta, colors)):

    label = r'$\beta=%s$' % fmt.format_data (beta1)

    # bz = +1
    ax[0].semilogx (+oc[0,:]/(2*Omega), gamma_max[0,:,i], color)
    oc0 = oc[0,:]
    kz0 = kz_opt[0,:,i]
    arg = np.argwhere (abs (np.diff (kz0)) > .01).flatten () + 1
    last = 0
    for pos in arg:
        ax[2].semilogx (oc0[last:pos]/(2*Omega), kz0[last:pos], color)
        last = pos
    plot = ax[2].semilogx (oc0[last:-1]/(2*Omega), kz0[last:-1], color,
            label = label)

    # bz = -1
    ax[1].semilogx (-oc[1,:]/(2*Omega), gamma_max[1,:,i])
    oc1 = oc[1,:]
    kz1 = kz_opt[1,:,i]
    cutoff = 15
    arg = np.argwhere (abs (np.diff (kz1[cutoff:])) > .1).flatten () + 1
    arg += cutoff
    last = 0
    for pos in arg:
        ax[3].semilogx (-oc1[last:pos]/(2*Omega), kz1[last:pos], color)
        last = pos
    ax[3].semilogx (-oc1[last:-1]/(2*Omega), kz1[last:-1], color)

ax[0].set_ylabel (r'$\gamma/\Omega$')
ax[2].set_ylabel (r'$k_z v_a/\Omega$')

for ax1 in ax[:2]:
    ax1.set_xticklabels ('')
    ax1.set_ylim (0, 2)
for ax1 in ax[2:]:
    ax1.set_xlabel (r'$\omega_cb_z/2\Omega$')
    ax1.set_ylim (0, 3)
for ax1 in ax[1::2]:
    ax1.set_yticklabels ('')
    ax1.invert_xaxis ()

legend = ax[2].legend (loc = 'upper left', ncol = 2)

scientific = lambda x, pos: '$-10^{%.2g}$' % np.log10 (x)
formatter = mpl.ticker.FuncFormatter (scientific)
ax[3].xaxis.set_major_formatter (formatter)
labels = ax[3].get_xticklabels ()
for label in labels: label.set_position ((0,1e-3))

fig.set_tight_layout ({'pad': 0., 'rect': (6e-3,0,1-1e-3,1-5e-3)})
helpers.savefig (fig, __file__)
