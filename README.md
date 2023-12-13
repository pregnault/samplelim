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
ease its use by ecological practitioners.

## Objective

`{samplelim}` aims at providing efficient implementations of two MCMC
algorithms for sampling high-dimensional polytopes; namely, the Mirror
Walk (MiW) introduced by Van Oevelen *et al.* (2010) and the Billard
Walk (BiW) introduced by Polyak and Gryazina (2014). It also provides
support functions to ease its use by ecologists practitioners interested
in solving LIM for trophic networks. Hence, `{samplelim}` can be viewed
as an updated, extended and low-level-encoded version of the R package
[`{limsolve}`](https://cran.r-project.org/web/packages/limSolve/index.html).

`{samplelim}` is built upon the C++ library
[`{volesti}`](https://github.com/GeomScale/volesti). Precisely, the
source code of the R package `{volesti}` 1.1.2-6 has been forked from
its [GitHub
repository](https://github.com/GeomScale/volesti/releases/tag/v1.1.2-6)
and has served as a basis for developing `{samplelim}`.

The C++ library `{volesti}` provides efficient implementations of
several MCMC algorithms for sampling high-dimensional polytopes and
estimating their volumes. Its R package counterpart includes a subset of
these algorithms among which, the BiW.

The modifications and additions made by `{samplelim}` (compared to
`{volesti}` 1.1.2-6) include:

- the removal of the bound on the number of reflections for the BiW ;
- the introduction of the MiW. More precisely, its implementation in
  `{samplelim}` is an optimized, C++ encoded version of the algorithm
  implemented in pure R programming in `{limsolve}` ;
- support functions inspired by `{limsolve}` making their users to be
  able to master `{samplelim}` very easily.

`{samplelim}` does not aim at replacing nor competing `{volesti}`, which
offers high-quality contents that are not covered by `{samplelim}`. We
expect `{samplelim}` to ease the junction between ecologists that are
familiar with `{limsolve}` and other communities interested in sampling
high-dimensional polytopes.

## Installation

The package `{samplelim}` is not yet available on CRAN mirrors. Its
source code is available on GitHub. It can be installed from its
repository by executing the following chunk in an R console.

``` r
if (! "remotes" %in% installed.packages()[,"Package"]) {install.packages("remotes")}
remotes::install_github("https://github.com/GrenteTheo/samplelim")
```

## Users’ guide

The main front-end function of `{samplelim}` is `rlim()`, performing
uniform sampling into a polytope by means of MCMC algorithm. A polytope
$\mathcal{P}$ is an $n$-dimensional convex set defined as the
intersection of hyper-planes and half-spaces.
<!-- It is mathematically characterized by two couples of a matrix and a vector, $(A, B)$ and $(G,H)$, with A of dimension $d\times ..$, $B$ of dimension $..$, $G$ of dimension $d \times ..$ and $H$ of dimension $..$, describing respectively the hyper-planes and half-spaces.  -->
More precisely,
$\mathcal{P} = \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}$, where $A$
is an $m\times n$ matrix, with $m \leq n$, $B \in \mathbb{R}^m$, $G$ is
a $k \times n$ matrix, with $k \geq ..$ and $H \in \mathbb{R}^k$.
Inequality constraints, represented by the matrix $G$, involved in
linear inverse models for metabolic networks make the polytope bounded.

The main argument `lim` of `rlim()` is a list with (at least) four
components named `A`, `B`, `G` and `H`, the matrices and vectors
defining the polytope to be sampled. Right below is a minimal example.

``` r
library("samplelim")
# Define equality and inequality constraints through matrices A, B, G, H
A <- matrix(rep(1, 3L), nrow = 1, ncol = 3)
B <- 1
G <- -matrix(c(0,0,1,
              0, 0, -1,
              0, 1, 0, 
              0, -1, 0,
              1, 0, 0, 
              -1, 0, 0),
            byrow = TRUE,
            nrow = 6, ncol = 3)
H <- -matrix(c(0.7, 0, 0.8, 0, 0.8, 0), nrow = 6)
# Store into a list
lim_exm <- list(A = A, B = B,G = G,H = H)
# Sampling into the polytope defined by these constraints
sample <- rlim(lim_exm, seed = 123)
# Show first points of the sample
head(sample)
```

    ##           [,1]       [,2]      [,3]
    ## [1,] 0.2002727 0.13542660 0.6643007
    ## [2,] 0.1601968 0.20447857 0.6353246
    ## [3,] 0.2434451 0.26803863 0.4885162
    ## [4,] 0.2464662 0.19465646 0.5588773
    ## [5,] 0.3173816 0.12140515 0.5612132
    ## [6,] 0.3147804 0.02503919 0.6601804

In the previous – pedagogical, example, we can visualize the sample
generated by `rlim()` thanks to a 3D scatterplot.

``` r
# Preparing sample for 3D scatterplot
colnames(sample) <- c("X", "Y","Z")
sample <- as.data.frame(sample)
library("plotly")
plot_sample <- plot_ly(sample, x = ~X, y = ~Y, z = ~Z, size = 0.2)
plot_sample <- plot_sample %>% add_markers()
plot_sample
```

A complete user guide will be included soon in the package, in the form
of a vignette.

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
