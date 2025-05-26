#' MCMC sampling algorithms for linear inverse models
#' 
#' @description
#' The package provides efficient implementations (C++ encoded) of 
#' Monte Carlo Markov Chains (MCMC) algorithms for uniformly sampling high-dimensional polytopes. 
#' It is particularly aimed at handling linear inverse models (LIM) in metabolic 
#' (trophic, biochemical or urban) networks. 
#' Some support functions are also included to facilitate its use by practitioners.
#' @keywords 
#' polytope
#' sampling
#' mcmc
#' lim
#' mirror-walk
#' billard-walk
#' trophic-network
#' urban-network
#' biochemical-network
#' metabolic-network
#' @references {
#' D. Van Oevelen, K. Van den Meersche, F. J. R. Meysman, K. Soetaert, 
#' J. J. Middelburg and A. F. Vézina, 
#' \emph{Quantifying Food Web Flows Using Linear Inverse Models},
#'  Ecosystems \strong{13}, 32-45 (2010).
#' 
#' B.T. Polyak and E.N. Gryazina, 
#' \emph{Billiard walk - a new sampling algorithm for control and optimization}, 
#' IFAC Proceedings Volumes, \strong{47(3)}, 6123-6128 (2014).
#' 
#' V. Girardin, T. Grente, N. Niquil and P. Regnault,
#' \emph{Comparing and updating R packages of MCMC Algorithms for 
#' Linear Inverse Modeling of Metabolic Networks}, to appear in
#' Communications in Statistics - Simulation and Computation (2025)
#' }
"_PACKAGE"