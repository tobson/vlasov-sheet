#!/usr/bin/env python

def setup (backend = "pdf"):

    from sys import modules

    settings = {}

    if 'matplotlib.backends' in modules:
        from matplotlib import get_backend
        if get_backend () != backend: settings['interactive'] = True
    else:
        from matplotlib import use
        from math import sqrt
        use (backend)

        # Font size
        settings['font.size'] = 10
        settings['legend.fontsize'] = 10
        settings['axes.titlesize'] = 10

        # Latex stuff
        settings['text.usetex'] = True
        settings['font.family'] = 'serif'
        settings['text.latex.preamble'] = [r'\usepackage[T1]{fontenc}',
            r'\usepackage{lmodern}', r'\usepackage{microtype}',
            r'\usepackage{amsmath}', r'\usepackage{bm}']

        # Figure size
        dpi = 72.27
        textwidth = 510.0
        width = textwidth/dpi
        height = width*(sqrt(5) - 1)/2

        settings['figure.dpi'] = dpi
        settings['figure.figsize'] = [width, height]

    from matplotlib import rcParams, is_interactive
    rcParams.update (settings)

    return is_interactive ()

def savefig (figure, scriptname):

    from matplotlib import get_backend
    from os import makedirs
    from os.path import dirname, realpath, exists, splitext, basename

    name, py = splitext (basename (scriptname))

    backend = get_backend ()

    ext = {'pdf': 'pdf', 'ps': 'eps'}[backend]

    figdir = dirname (realpath (__file__)) + '/../figures'
    if not exists (figdir):
        makedirs (figdir)

    figure.savefig (figdir + '/' + name + '.' + ext)

    if backend == 'ps':
        from os import system
        system ('epstopdf --hires ../figures/' + name + '.eps')

def load (name):

    from sys import argv
    from os.path import dirname, realpath
    from numpy import any, array, load

    if any (array (argv) == 'recompute'):
        raise

    datadir = dirname (realpath (__file__)) + '/../data'
    datafile = name + '.npz'

    try:
        data = []
        with load (datadir + '/' + datafile) as npz:
            for key, value in sorted (npz.iteritems ()):
                data.append (value)
        return data
    except:
        print ('Data file "' + datafile + '" not found.')
        raise

def save (name, *args):

    from os import makedirs
    from os.path import dirname, exists, realpath
    from numpy import savez

    datadir = dirname (realpath (__file__)) + '/../data'
    if not exists (datadir):
        makedirs (datadir)
    datafile = name + '.npz'

    savez (datadir + '/' + datafile, *args)
