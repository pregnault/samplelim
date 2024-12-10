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
efficient implementations (C++ encoded) of Monte Carlo Markov Chains
(MCMC) algorithms for uniformly sampling high-dimensional polytopes. It
is particularly aimed at handling linear inverse models (LIM) in
metabolic (trophic, biochemical or urban) networks. Some support
functions are also included to facilitate its use by practitioners.

## Objective

`{samplelim}` provides efficient implementations of two MCMC algorithms
for sampling high-dimensional polytopes, the Mirror Walk (MiW)
introduced by Van Oevelen *et al.* (2010) and the Billard Walk (BiW)
introduced by Polyak and Gryazina (2014). Thanks to the inclusion of
easy to use support functions for linear inverse modeling of metabolic
networks, `{samplelim}` can be viewed as an updated, extended and
low-level-encoded version of the R package
[`{limsolve}`](https://cran.r-project.org/web/packages/limSolve/index.html).

`{samplelim}` is built upon the C++ library
[`{volesti}`](https://github.com/GeomScale/volesti): the source code of
the R package `{volesti}` 1.1.2-6 has been forked from its [GitHub
repository](https://github.com/GeomScale/volesti/releases/tag/v1.1.2-6)
as a basis for developing `{samplelim}`.

The C++ library `{volesti}` provides efficient implementations of
different MCMC algorithms for sampling high-dimensional polytopes and
estimating their volumes. Its R package counterpart proposes a subset of
these algorithms among which the BiW.

`{samplelim}` aims at combining the performance of `{volesti}` and the
convenient features of `{limsolve}`. Precisely :

- the MiW is implemented. It is an optimized, C++ encoded version of the
  algorithm coded in pure R programming language in `{limsolve}`;
- `{samplelim}` includes and slightly modifies the BiW from `{volesti}`
  1.1.2-6 – the uniform distribution of the path length in `{volesti}`
  is replaced by the exponential distribution, as originally suggested
  by Polyak and Gryazina (2014);
- support functions allowing an easy handling by users are included.
  These functions are updated/modified versions of the ones present in
  `{limsolve}`.

## Installation

The package `{samplelim}` is not yet available on CRAN mirrors. Its
source code is available on GitHub. It can be installed from its
repository by executing the following chunk in any R console.

``` r
# Install package remotes if not already installed
if (! "remotes" %in% installed.packages()[,"Package"]) {install.packages("remotes")}
# Install package samplelim from its GitHub repo
remotes::install_github("https://github.com/pregnault/samplelim")
```

## Typical workflow

The workflow of `{samplelim}` is greatly inspired by the one of
`{limsolve}`. The main front-end function of `{samplelim}` is `rlim()`,
performing uniform sampling in the polytope associated to LIM, by means
of MCMC algorithm. Its main, mandatory, argument is `lim`, a list or an
object of class `lim` (introduced in `{limsolve}`) encompassing the
description of the polytope to be sampled. This list or lim object can
be defined manually or, more suitably, from a description file, as
illustrated in the following chunk.

``` r
# Load package
library("samplelim")
# Find path to the example of declaration file (DF) included in samplelim
DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
# Read DF and create a lim object
BOWF <- df2lim(DF)
```

Then, sampling is performed by a simple call to the function `rlim()`,
as follows.

``` r
sample <- rlim(lim = BOWF, 
               seed = 123, # Set the seed of PRNG
               nsamp = 5000) # Number of points in the returned sample
# The points are presented in an N*n matrix, where
# N is the number of sampled points (here, 5000, default = 3000)
# n is the number of flows (ambiant space of the polytope)
dim(sample)
```

    ## [1] 5000   28

Diagnostics on sampling performances can then be performed using package
`{coda}`, e.g., the Raftery and Lewis diagnostics, as illustrated below.

``` r
coda::raftery.diag(data = sample)
```

    ## 
    ## Quantile (q) = 0.025
    ## Accuracy (r) = +/- 0.005
    ## Probability (s) = 0.95 
    ##                                                 
    ##           Burn-in  Total Lower bound  Dependence
    ##           (M)      (N)   (Nmin)       factor (I)
    ##  FIX->PHY 18       19902 3746         5.310     
    ##  IMP->FBF 2        3680  3746         0.982     
    ##  DET->BIV 2        3803  3746         1.020     
    ##  DET->ZOO 2        3803  3746         1.020     
    ##  DET->BAC 18       20799 3746         5.550     
    ##  PHY->RES 12       14451 3746         3.860     
    ##  PHY->DET 15       17676 3746         4.720     
    ##  PHY->BIV 2        3680  3746         0.982     
    ##  PHY->ZOO 2        3930  3746         1.050     
    ##  PHY->BAC 2        3741  3746         0.999     
    ##  PHY->LOS 2        3803  3746         1.020     
    ##  BAC->RES 8        10712 3746         2.860     
    ##  BAC->DET 4        5211  3746         1.390     
    ##  BAC->LOS 8        11704 3746         3.120     
    ##  ZOO->RES 3        4062  3746         1.080     
    ##  ZOO->DET 2        3803  3746         1.020     
    ##  ZOO->FBF 2        3741  3746         0.999     
    ##  ZOO->BIV 2        3680  3746         0.982     
    ##  ZOO->ZOO 2        3930  3746         1.050     
    ##  ZOO->LOS 2        3803  3746         1.020     
    ##  BIV->RES 2        3803  3746         1.020     
    ##  BIV->DET 2        3620  3746         0.966     
    ##  BIV->FBF 2        3741  3746         0.999     
    ##  BIV->LOS 2        3741  3746         0.999     
    ##  FBF->RES 2        3803  3746         1.020     
    ##  FBF->DET 2        3741  3746         0.999     
    ##  FBF->FBF 2        3680  3746         0.982     
    ##  FBF->LOS 2        3620  3746         0.966

A complete user guide will soon be included in the package, in the form
of a vignette.

A comparison study of computation time and sampling quality, between the
implementations of the MiW of `{samplelim}`, the BiW of `{volesti}`
(1.1.2-6) and the Coordinate Hit-and-Run with Rounding (CHRR) of the
MatLab library `{COBRA}` is available in Girardin *et al.* (2024).

## Credits

The modifications of the core C++ implementation of the BiW and the C++
implementation of the MiW, have been performed by Matthieu DIEN and Théo
GRENTE.

The R packaging has been performed by Jacques BRÉHÉLIN, Théo GRENTE and
Philippe REGNAULT.

The Declaration File in `inst/extdata` has been produced by Quentin
NOGUÈS; see Noguès *at al.* (2020) for details on the ecological network
it relies on.

We refer to the [`credits.md`
file](https://github.com/GeomScale/volesti/blob/v1.1.1/doc/credits.md)
of the R package `{volesti}` for additional information on the
development of the core C++ implementation of the BiW (and other
algorithms implemented in `{volesti}`).

The function `df2lim()` (reading and formatting a declaration file) is a
wrapped copy of functions `Read()` and `Setup()` from package `{LIM}`,
developped by Karline SOETAERT. It has been added to the present package
to limit its dependency tree.

## Licensing

You may redistribute or modify the software under the [GNU Lesser
General Public License](LICENSE.md) as published by Free Software
Foundation, either version 3 of the License, or any later version. It is
distributed in the hope that it may be helpful to the community, but
WITHOUT ANY WARRANTY.

## References

V. Girardin, T. Grente, N. Niquil and P. Regnault, *Comparing and
updating R packages using MCMC Algorithms for Linear Inverse Modeling of
Metabolic Networks*, [hal-04455831](https://hal.science/hal-04455831)
(2024).

Q. Noguès, A. Raoux, E. Araignous, T. Hattab, B. Leroy, F. Ben Rais
Lasram, F. Le Loc’h, J. Dauvin and N. Niquil, *Cumulative effects of
marine renewable energy and climate change on ecosystem properties:
Sensitivity of ecological network analysis*, Ecological Indicators
**121**, 107128 (2020).

L. Cales, A. Chalkis, I.Z. Emiris and V. Fisikopoulos, *Practical volume
computation of structured convex bodies, and an application to modeling
portfolio dependencies and financial crises*, Proc. of Symposium on
Computational Geometry, Budapest, Hungary (2018).

B.T. Polyak and E.N. Gryazina, *Billiard walk - a new sampling algorithm
for control and optimization*, IFAC Proceedings Volumes, **47(3)**,
6123-6128 (2014).

D. Van Oevelen, K. Van den Meersche, F. J. R. Meysman, K. Soetaert, J.
J. Middelburg and A. F. Vézina, *Quantifying Food Web Flows Using Linear
Inverse Models*, Ecosystems **13**, 32-45 (2010).
