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
is particularly convenient for solving linear inverse models (LIM) in
metabolic (trophic, biochemical or urban) networks. Particularly, some
support functions inspired by the package `{limsolve}` are designed to
facilitate its use by practitioners.

## Objective

`{samplelim}` aims at providing efficient implementations of two MCMC
algorithms for sampling high-dimensional polytopes, the Mirror Walk
(MiW) introduced by Van Oevelen *et al.* (2010) and the Billard Walk
(BiW) introduced by Polyak and Gryazina (2014). It also provides support
functions to facilitate its use by practitioners interested in solving
LIM for networks. In this regard, `{samplelim}` can be viewed as an
updated, extended and low-level-encoded version of the R package
[`{limsolve}`](https://cran.r-project.org/web/packages/limSolve/index.html).

`{samplelim}` is built upon the C++ library
[`{volesti}`](https://github.com/GeomScale/volesti). Precisely, the
source code of the R package `{volesti}` 1.1.2-6 has been forked from
its [GitHub
repository](https://github.com/GeomScale/volesti/releases/tag/v1.1.2-6)
as a basis for developing `{samplelim}`.

The C++ library `{volesti}` provides efficient implementations of
different MCMC algorithms for sampling high-dimensional polytopes and
estimating their volumes. Its R package counterpart includes a subset of
these algorithms including the BiW.

`{samplelim}` aims at combining the performance of `{volesti}` and the
convenient features of `{limsolve}`. More precisely :

- `{samplelim}` copies and slightly modifies the BiW from `{volesti}`
  1.1.2-6 – the uniform distribution of the path length in `{volesti}`
  is replaced by the exponential distribution, as originally suggested
  by Polyak and Gryazina (2014);
- the MiW is implemented. It is an optimized, C++ encoded version of the
  algorithm implemented in pure R programming language in `{limsolve}`;
- support functions inspired by `{limsolve}` add included, that allow
  users an easy handling.

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

The workflow of `{samplelim}` is greatly inspired by the package
`{limsolve}`. The main front-end function of `{samplelim}` is `rlim()`,
performing uniform sampling into the polytope associated to a LIM, by
means of MCMC algorithm. Its main, mandatory, argument is `lim`, a list
or an object of class `lim` (introduced in `{limsolve}`) encompassing
the description of the polytope to be sampled. This list or lim object
can be defined by hand or, more suitably, from a description file, as
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
    ##  Burn-in  Total Lower bound  Dependence
    ##  (M)      (N)   (Nmin)       factor (I)
    ##  15       19070 3746         5.090     
    ##  2        3803  3746         1.020     
    ##  2        3741  3746         0.999     
    ##  2        3803  3746         1.020     
    ##  16       18132 3746         4.840     
    ##  12       13308 3746         3.550     
    ##  16       17556 3746         4.690     
    ##  2        3741  3746         0.999     
    ##  2        3866  3746         1.030     
    ##  2        3680  3746         0.982     
    ##  2        3620  3746         0.966     
    ##  8        9730  3746         2.600     
    ##  8        10664 3746         2.850     
    ##  9        12927 3746         3.450     
    ##  2        3680  3746         0.982     
    ##  2        3803  3746         1.020     
    ##  2        3803  3746         1.020     
    ##  2        3680  3746         0.982     
    ##  2        3680  3746         0.982     
    ##  2        3930  3746         1.050     
    ##  2        3741  3746         0.999     
    ##  2        3680  3746         0.982     
    ##  2        3866  3746         1.030     
    ##  2        3803  3746         1.020     
    ##  2        3803  3746         1.020     
    ##  2        3620  3746         0.966     
    ##  2        3680  3746         0.982     
    ##  2        3866  3746         1.030

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

The R packaging has been performed by Théo GRENTE and Philippe REGNAULT.

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
updating R packages of MCMC Algorithms for Linear Inverse Modeling of
Metabolic Networks*, hal: (2024).

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
