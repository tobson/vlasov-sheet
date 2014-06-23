The paper can be compiled simply by running `make` in the top directory. This
uses [revtex] [1] and [latexmk] [2]. All other LaTeX packages needed such as
[cleveref] [3] and [tikz] [4] are listed in `preamble.tex`.

Every figure in the paper except for `stix.pdf` is produced by a corresponding
Python script. The required Python libraries such as [numpy] [5], [cython] [6],
and [matplotlib] [7] are listed in `python/requirements.txt`.

[1]: https://journals.aps.org/revtex  "REVTeX"
[2]: http://www.ctan.org/pkg/latexmk/ "latexmk"
[3]: http://www.ctan.org/pkg/cleveref "cleveref"
[4]: http://www.texample.net/tikz/    "TikZ"
[5]: http://www.numpy.org/            "NumPy"
[6]: http://cython.org/               "Cython"
[7]: http://matplotlib.org/           "matplotlib"
