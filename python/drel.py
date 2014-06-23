from __future__ import division
import numpy as np
from DSolveParallel import DSolveParallel
from findmax import findmax
import helpers
pylab = helpers.setup ('ps')
import matplotlib.pyplot as plt

def drel_oc (axis):

    # How many wave numbers?
    nkz = 256

    # Plasma beta
    beta = 1000.0

    # Cyclotron frequencies
    oc = 2.*(10.**np.array ([-1, 0, 1, 2, 4]))

    # Arrays to hold wave numbers and growth rates
    kz = np.empty (oc.shape + (nkz,), dtype = np.double)
    gamma = np.empty (oc.shape + (nkz,), dtype = np.double)

    for i, oc1 in np.ndenumerate (oc):

        # Dispersion relation solver
        dsolve = DSolveParallel (oc1)

        # Range of unstable wave numbers
        kz1 = np.linspace (3e-3*dsolve.kmax, (1-1e-5)*dsolve.kmax, nkz)

        # Cold plasma growth rate. Use this as initial guess
        cold = dsolve.cold (kz1)
        gamma1 = cold.copy ()

        # Start of beta sweep
        beta_start = 1e-2
        # Increment
        beta_incr = 1e-2
        ratio = beta/beta_start
        nbeta = np.log10 (ratio)/beta_incr
        # Fine grained loop over beta
        for beta1 in beta_start*pow (ratio, np.arange (nbeta + 1)/nbeta):
            try:
                gamma1 = dsolve.warm (beta1, kz1, gamma1)
            except RuntimeError as e:
                print e.message
                print "oc = %g, beta = %g" % (oc1, beta1)
                raise SystemExit
        kz[i] = kz1
        gamma[i] = gamma1

    for kz1, gamma1, oc1 in zip (kz, gamma, oc):
        label = r'$10^{%.2g}$' % np.log10 (oc1/2.)
        axis.plot (kz1, gamma1, label = label)

    axis.legend (loc = 'upper right', title = r'$\omega_cb_z/2\Omega$')
    axis.text (.775, 1.6, r'$\beta=%g$' % beta)
    axis.set_ylim (0, 1.8)
    axis.set_ylabel (r'$\gamma/\Omega$')

def drel_beta (axis):

    # How many wave numbers?
    nkz = 256

    # Plasma beta
    oc = 100.0

    # Dispersion relation solver
    dsolve = DSolveParallel (oc)

    # Range of unstable wave numbers
    kz = np.linspace (3e-3*dsolve.kmax, (1-1e-5)*dsolve.kmax, nkz)

    # Cold plasma growth rate. Use this as initial guess
    cold = dsolve.cold (kz)
    gamma1 = cold.copy ()

    # Cyclotron frequencies
    beta = 10.**np.arange (1, 5)

    # Array to hold growth rates
    gamma = np.empty (beta.shape + (nkz,), dtype = np.double)

    # Start of beta sweep
    beta_start = 1e-2
    # Increment
    beta_incr = 1e-2

    for j, beta1 in np.ndenumerate (beta):
        ratio = beta1/beta_start
        nbeta = np.log10 (ratio)/beta_incr
        # Fine grained loop over beta
        for beta2 in beta_start*pow (ratio, np.arange (nbeta + 1)/nbeta):
            try:
                gamma1 = dsolve.warm (beta2, kz, gamma1)
            except RuntimeError as e:
                print e.message
                print "oc = %g, beta = %g" % (oc, beta2)
                raise SystemExit
        # Save result
        gamma[j] = gamma1
        # Start next fine grained loop at this value
        beta_start = beta1

    # Right panel
    for gamma1, beta1 in zip (gamma, beta):
        label = r'$10^{%.2g}$' % np.log10 (beta1)
        axis.plot (kz, gamma1, label = label)
    axis.legend (loc = 'upper right', title = r'$\beta$')
    axis.text (.775, 1.6, r'$\omega_cb_z=%g\,\Omega$' % oc)

def drel_negative (axes):

    def compute ():

        # How many wave numbers?
        nkz = 2048

        beta = np.array ([1e-4, 1e-2, 5e-2, 1e-1, 1e0, 1e2, 1e3, 1e4, 1e6])
        dsolve = DSolveParallel (oc)

        # Range of unstable wave numbers
        kz_min = 1e-3
        kz_max = (1 - 1e-3)*dsolve.kmax
        kz = np.logspace (np.log10 (kz_min), np.log10 (kz_max), nkz)

        # Start of beta sweep
        beta_start = 1e-2
        # Increment
        beta_incr = 1e-3

        # Array to hold growth rate
        gamma = np.empty (beta.shape + kz.shape, dtype = np.double)
        # Dispersion relation solver
        dsolve = DSolveParallel (oc)
        # Cold plasma growth rate. Use this as initial guess
        cold = dsolve.cold (kz)
        guess = cold.copy ()
        # Coarse loop over beta
        for j, beta1 in np.ndenumerate (beta):
            ratio = beta1/beta_start
            nbeta = np.log10 (ratio)/beta_incr
            # Fine grained loop over beta
            for beta2 in beta_start*pow (ratio, np.arange (nbeta + 1)/nbeta):
                try:
                    guess = dsolve.warm (beta2, kz, guess)
                except RuntimeError as e:
                    print e.message
                    print "oc = %g, beta = %g" % (oc, beta2)
                    raise SystemExit
            # Save result
            gamma[j] = guess
            # Start next fine grained loop at this value
            beta_start = beta1

        return beta, kz, gamma, cold

    # Cyclotron frequency
    oc = -2.1

    data = 'drel-negative'
    try: 
        beta, kz, gamma, cold = helpers.load (data)
    except:
        print 'Recomputing...'
        beta, kz, gamma, cold = compute ()
        helpers.save (data, beta, kz, gamma, cold)

    ax = axes[0]
    labels = [r'$10^{%.2g}$' % np.log10 (beta1) for beta1 in beta[3:]]
    for beta1, gamma1 in zip (beta[3:], gamma[3:]):
        label = r'$10^{%.2g}$' % np.log10 (beta1)
        ax.semilogx (kz, gamma1, label = label)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend (handles[::-1], labels[::-1], loc = 'upper left',
            title = r'$\beta$')
    ax.text (.06, 13, r'$\omega_cb_z=%g\Omega$' % oc)

    ax = axes[1]
    ax.plot (kz, cold, 'k', label = r'$\mathrm{Cold}$')
    ax.plot (kz, gamma[0], 'r--', label = r'$10^{-4}$')
    for beta1, gamma1, col in zip (beta[1:4], gamma[1:4], ['b','g','m']):
        ax.plot (kz, gamma1, col, label = beta1)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend (handles[::-1], labels[::-1],
            bbox_to_anchor = (.075,-.01,1.,1.), loc = 'lower left',
            title = r'$\beta$')
    ax.text (5.1, .975, r'$\omega_cb_z=%g\Omega$' % oc)

    for ax in axes:
        ax.set_xlabel (r'$k_z v_a/\Omega$')
    axes[0].set_ylabel (r'$\gamma/\Omega$')

# Units
Omega = 1.
va = 1.

# Logarithmic rate of shear
qshear = 1.5

# Create figure
width, height = plt.rcParams['figure.figsize']
height = .8*width
fig = plt.figure (num = 1, figsize = (width, height))
fig.clf ()
fig, axes = plt.subplots (num = 1, nrows = 2, ncols = 2)

# Plots
drel_oc (axes[0,0])
drel_beta (axes[0,1])
drel_negative (axes[1,:])

# Save figure
if not pylab:
    tight = dict (pad = 0, h_pad = .5, w_pad = 1.,
            rect = [.0075,0,.9975,.995])
    fig.set_tight_layout (tight)
    helpers.savefig (fig, __file__)
