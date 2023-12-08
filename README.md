samplelim : MCMC sampling algorithms for linear inverse models in R
================

![](https://img.shields.io/badge/lifecycle-maturing-blue.svg) [![Project
Status: Active – The project has reached a stable, usable state and is
being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: LGPL
v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/regexplain)](https://cran.r-project.org/package=samplelim)

The package `{samplelim}` for the R statistical software provides
performant implementations (C++ encoded) of Monte Carlo Markov Chains
(MCMC) algorithms for uniformly sampling hig-dimensional polytopes. It
is particularly convenient for solving linear inverse models (LIM) in
metabolic (trophic, biochemical or urban) networks. Particularly, some
support functions inspired by the package limsolve are designed to ease
its use by ecological practitioners.

## Objective

`{samplelim}` aims at providing performant implementations of two MCMC
algorithms for sampling high-dimensional polytopes; namely, the Mirror
Walk (MiW) introduced by Van Oevelen *et al.* (2010) and the Billard
Walk (BiW) introduced by Chalkis *et al.* (2014). It also provides
support functions to ease its use by ecologists practitioners interested
in solving LIM for trophic networks. Hence, `{samplelim}` can be viewed
as an updated, extended and low-level-encoded version of the R package
[`{limsolve}`](https://cran.r-project.org/web/packages/limSolve/index.html).

`{samplelim}` is built upon the C++ library
[`{volesti}`](https://github.com/GeomScale/volesti). Precisely, the
source code of the R package `{volesti}` 1.1.2-6 has been forked from
its [GitHub
repository](https://github.com/GeomScale/volesti/releases/tag/v1.1.2-6)
and has served as a basis for developping `{samplelim}`.

The C++ library `{volesti}` provides efficient implementations of
several MCMC algorithms for sampling high-dimensional polytopes and
estimating their volumes. Its R package counterpart includes a subset of
these algorithms among which, the BiW.

The modifications and additions made by `{samplelim}` (compared with
`{volesti}` 1.1.2-6) include:

- the removal of the bound on the number of reflections for the BiW ;
- the introduction of the MiW. More precisely, its implementation in
  `{samplelim}` is an optimized, C++ encoded version of the algorithm
  implemented in pure R programming in `{limsolve}` ;
- support functions inspired by `{limsolve}` making their users to be
  able to master `{samplelim}` very easily.

`{samplelim}` does not aim at replacing nor concurrencing `{volesti}`,
which offers high-quality contents that are not covered by
`{samplelim}`. We expect `{samplelim}` to ease the junction between
ecologists that are familiar with `{limsolve}` and other communities
interested in sampling high-dimensional polytopes.

## Installation

The package `{samplelim}` is not yet available on the CRAN mirrors. Its
source code is available on GitHub. It can be installed from this repo
by executing the following chunk in an R console.

``` r
if (! "remotes" %in% installed.packages()[,"Package"]) {install.packages("remotes")}
remotes::install_github("https://github.com/GrenteTheo/samplelim")
```

## Users’ guide

ADD 2 or 3 EXAMPLES.

A complete user guide has been included in the package, in the form of a
vignette.

A comparison study, including computation time and sampling quality,
between the implementations of the MiW of `{samplelim}`, the BiW of
`{volesti}` (1.1.2-6) and the Coordinate Hit-and-Run with Rounding
(CHRR) of the MatLab library `{COBRA}`, has been performed; see Girardin
*et al.* (2023).

## Credits

The modifications of the core C++ implementation of the BiW, so as the
C++ implementation of the MiW, have been performed by Matthieu DIEN and
Théo GRENTE.

The R packaging has been performed by Théo GRENTE and Philippe REGNAULT.

The Declaration Files in `inst/extdata` have been produced by Quentin
NOGUÈS; see Noguès *at al.* (2020) for details on the ecological network
they rely on.

We refer to the [`credits.md`
file](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/credits.md)
of the R package `{volesti}` for additional information on the
development of the core C++ implementation of the BiW (and other
algorithms implemented in `{volesti}`).

## Licensing

You may redistribute or modify the software under the [GNU Lesser
General Public License](LICENSE.md) as published by Free Software
Foundation, either version 3 of the License, or (at your option) any
later version. It is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY.

## References

V. Girardin, T. Grente, N. Niquil and P. Regnault, *Comparing and
updating R packages of MCMC Algorithms for Linear Inverse Modeling of
Metabolic Networks*, hal: (2023)

Q. Noguès, A. Raoux, E. Araignous, T. Hattab, B. Leroy, F. Ben Rais
Lasram, F. Le Loc’h, J. Dauvin and N. Niquil, *Cumulative effects of
marine renewable energy and climate change on ecosystem properties:
Sensitivity of ecological network analysis*, Ecological Indicators
**121**, 107128 (2020)

L. Cales, A. Chalkis, I.Z. Emiris and V. Fisikopoulos, *Practical volume
computation of structured convex bodies, and an application to modeling
portfolio dependencies and financial crises*, Proc. of Symposium on
Computational Geometry, Budapest, Hungary (2018).

B.T. Polyak and E.N. Gryazina, *Billiard walk - a new sampling algorithm
for control and optimization*, IFAC Proceedings Volumes, **47(3)**,
6123-6128 (2014)

D. Van Oevelen, K. Van den Meersche, F. J. R. Meysman, K. Soetaert, J.
J. Middelburg and A. F. Vézina, *Quantifying Food Web Flows Using Linear
Inverse Models*, Ecosystems **13**, 32-45 (2010)
