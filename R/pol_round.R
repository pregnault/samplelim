#' Round a polytope
#'
#' The functions \code{pol.round()} and \code{lim.round()} round a given polytope
#' \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Gx \geq H \}}: they apply to it the affine change
#' of coordinates \eqn{x = x_0 + Z y} that maps the Dikin ellipsoid at its analytic center
#' \eqn{x_0} onto the unit ball, so that the rounded polytope is close to isotropic.
#'
#' The Dikin ellipsoid at \eqn{x_0} is \eqn{\{x : (x - x_0)^\top M (x - x_0) \leq 1\}},
#' with \eqn{M = \sum_i g_i g_i^\top / (\langle g_i, x_0 \rangle - H_i)^2}, \eqn{g_i} being
#' the \eqn{i}-th row of \code{G}, and \eqn{Z = M^{-1/2}}. The change of coordinates being
#' affine, \code{\link{red2full}()} maps a uniform sample of the rounded polytope onto a
#' uniform sample of \eqn{\mathcal{P}}. Redundant constraints distort the ellipsoid, hence
#' \code{\link{pol.exfoliate}()} should be applied first.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope \eqn{\mathcal{P}}.
#'
#' @return For \code{pol.round()}, a list with four components; namely:
#' \itemize{
#' \item \code{G} and \code{H}, describing the rounded polytope
#' \eqn{\{ y : GZy \geq H - Gx_0 \}};
#' \item \code{x0} and \code{Z}, the change of coordinates \eqn{x = x_0 + Z y}.
#' }
#' For \code{lim.round()}, \code{lim} with these four components replaced. When \code{lim}
#' is returned by \code{\link{lim.redpol}()}, \code{x0} and \code{Z} are composed with the
#' reduction, so that \code{\link{red2full}()} maps the rounded polytope onto the original
#' one, with its equality constraints.
#' @export
#'
#' @rdname pol.round
#' @references I. I. Dikin,
#' \emph{Iterative solution of problems of linear and quadratic programming},
#' Soviet Mathematics Doklady \strong{8}, 674-675 (1967).
pol.round <- function(G, H) {
  P <- .pol_check(G, H); G <- P$G; H <- P$H
  x0 <- pol.center(G = G, H = H)
  d <- as.numeric(G %*% x0) - H
  nz <- rowSums(G != 0) > 0                                # zero rows play no part
  e <- eigen(crossprod(G[nz, , drop = FALSE], G[nz, , drop = FALSE] / d[nz]^2),
             symmetric = TRUE)                             # M, positive definite
  Z <- e$vectors %*% (t(e$vectors) / sqrt(e$values))     # M^(-1/2)
  list(G = G %*% Z, H = -d, x0 = x0, Z = Z)
}


#' @param lim A list with components \code{G} and \code{H} describing a reduced polytope,
#'   as returned by \code{\link{lim.redpol}()}. A list with equality constraints
#'   (component \code{A}) is refused.
#' @export
#' @rdname pol.round
lim.round <- function(lim) {
  .lim_check_reduced(lim)
  rnd <- pol.round(G = lim$G, H = lim$H)
  if (!is.null(lim$x0) && !is.null(lim$Z)) {   # x = x0 + Z (x0' + Z' y)
    rnd$x0 <- as.numeric(lim$x0 + lim$Z %*% rnd$x0)
    rnd$Z <- lim$Z %*% rnd$Z
  }
  lim[c("G", "H", "x0", "Z")] <- rnd[c("G", "H", "x0", "Z")]
  lim
}
