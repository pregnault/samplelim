

#' Round a polytope using the Dikin ellipsoid
#'
#' A very elongated polytope is disastrous for a random walk: a billiard trajectory
#' bounces without advancing, and an isotropic Gaussian jump is accepted along the short
#' directions and rejected along the long ones. \emph{Rounding} changes coordinates so
#' that the body becomes as isotropic as possible.
#'
#' The Dikin ellipsoid at an interior point \eqn{c} is
#' \eqn{E(c) = \{x : (x-c)^\top \mathcal{H}(c) (x-c) \leq 1\}} with
#' \eqn{\mathcal{H}(c) = \sum_i g_i g_i^\top / d_i(c)^2} and
#' \eqn{d_i(c) = \langle g_i, c \rangle - H_i}. It is inscribed in the polytope for
#' \emph{any} interior \eqn{c}. Applying \eqn{x = c + T x'} with
#' \eqn{T = \mathcal{H}(c)^{-1/2}} maps it onto the unit ball, which rounds the body.
#'
#' \strong{Why this is legitimate.} The map is affine, so its Jacobian \eqn{|\det T|} is
#' constant. The image of the uniform distribution on the rounded polytope is therefore
#' exactly the uniform distribution on the original one: sample in the rounded space, then
#' apply \code{back}, and no correction or reweighting is needed. This would be false for
#' a non-uniform target, where the density itself would have to be transported.
#'
#' \strong{Where to take the centre matters far more than one would think.} Since
#' \eqn{E(c)} is inscribed for any interior \eqn{c}, one could round from the Chebyshev
#' centre, which is readily available. Measured on the exfoliated BOWF-short model, in
#' minimum effective sample size per second:
#'
#' \tabular{lrrr}{
#'   \tab BiW \tab MiW \tab BPS \cr
#'   no rounding \tab 51 \tab 77 \tab 44 \cr
#'   rounded at the Chebyshev centre \tab 74 \tab 125 \tab 107 \cr
#'   rounded at the analytic centre \tab 14353 \tab 12262 \tab 19272
#' }
#'
#' Rounding from the Chebyshev centre buys a factor of 1.5 to 2.4; from the analytic
#' centre, a factor of 160 to 440. The reason is geometric: the Chebyshev centre only sees
#' the nearest faces, and the ellipsoid taken there is small and nearly round, so it hugs
#' nothing. The ellipsoid at the analytic centre is elongated \emph{like the body}, which
#' is precisely what lets it straighten it out.
#'
#' Apply \code{\link{pol.exfoliate}} first: \eqn{\mathcal{H}} is a sum over the
#' constraints, so a redundant one adds a term that pulls the ellipsoid for no geometric
#' reason. On BOWF-short, exfoliating first gains a factor of 55 in ellipsoid volume.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope
#'   \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Gx \geq H \}}. Exfoliate first.
#' @param H A numeric vector corresponding to \code{H}.
#' @param center Either \code{"analytic"} (the default, strongly recommended) or
#'   \code{"chebyshev"}. The latter exists only so that the comparison above can be
#'   reproduced; it should not be used in production.
#' @param x0 A strictly interior starting point, passed to \code{\link{pol.center}}.
#'
#' @return A list with components:
#'   \item{G, H}{the rounded polytope \eqn{\{x' : G' x' \geq H'\}}, with
#'     \eqn{G' = G T} and \eqn{H' = H - G c}. Rows are \strong{not} normalised;}
#'   \item{center}{the centre used, in the original coordinates;}
#'   \item{T, Tinv}{the transformation and its inverse, \eqn{x = center + T x'};}
#'   \item{forth}{a function mapping original coordinates to rounded ones;}
#'   \item{back}{a function mapping rounded coordinates back to the original space.
#'     \strong{This is the one to apply to a sample.} It accepts a vector or a matrix with
#'     one point per row, that is, the output of \code{rlim()} as is;}
#'   \item{axis.ratio}{the ratio of the semi-axes of the Dikin ellipsoid. Descriptive
#'     only: it measures the \emph{shape of the ellipsoid}, not the quality of the
#'     rounding. Judge the latter with \code{pol.ranges()} on the rounded polytope;}
#'   \item{log.volume}{the log-volume of the ellipsoid, up to an additive constant. This
#'     one does measure quality: the larger, the better the ellipsoid hugs the body.}
#'
#' @export
#'
#' @rdname pol.round
#' @examples
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' red <- lim.redpol(BOWF)
#' exf <- pol.exfoliate(G = red$G, H = red$H)
#' rnd <- pol.round(G = exf$G, H = exf$H)
#'
#' # anisotropy, before and after
#' rg0 <- pol.ranges(G = exf$G, H = exf$H)[, 3]
#' rg1 <- pol.ranges(G = rnd$G, H = rnd$H)[, 3]
#' c(before = max(rg0) / min(rg0), after = max(rg1) / min(rg1))
#'
#' # Sample in the rounded space, then come back. Note that Hpolytope uses the
#' # convention Ax <= b, hence the change of sign, and that rlim() requires `lim`
#' # to be passed explicitly even when a Hpolytope is supplied.
#' smp <- rlim(lim = NULL, Hpol = Hpolytope(A = -rnd$G, b = -rnd$H), nsamp = 100)
#' pts <- rnd$back(smp)
pol.round <- function(G, H, center = c("analytic", "chebyshev"), x0 = NULL) {
  center <- match.arg(center)
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(G)) G <- t(G)
  if (is.null(G)) stop("G is NULL, the polytope has 0 dimensions.")
  H <- as.numeric(H)
  if (nrow(G) != length(H)) stop("G and H have incompatible dimensions.")
  if (any(!is.finite(G)) || any(!is.finite(H))) stop("G and H must be finite.")

  # Internally: Gx >= H becomes (-G) x <= (-H). Rows are normalised so that no
  # constraint of large norm dominates numerically; the Dikin ellipsoid itself is
  # invariant under row rescaling, since g_i -> s g_i implies d_i -> s d_i.
  A <- -G; b <- -H
  s <- sqrt(rowSums(A^2))
  if (any(s <= 0)) stop("Degenerate constraint: some row of G is zero.")
  An <- A / s; bn <- b / s
  d <- ncol(An)

  ctr <- if (identical(center, "analytic")) {
    pol.center(G = -An, H = -bn, type = "analytic", x0 = x0)
  } else {
    .pol_chebyshev(An, bn)$center
  }
  dc <- bn - as.numeric(An %*% ctr)
  if (any(dc <= 0)) stop("The chosen centre is not interior.")

  Hess <- crossprod(An, An * (1 / dc^2))
  eig <- eigen(Hess, symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  Tm   <- eig$vectors %*% (t(eig$vectors) / sqrt(vals))   # H^{-1/2}: rounds
  Tinv <- eig$vectors %*% (t(eig$vectors) * sqrt(vals))   # H^{+1/2}: comes back
  semi <- 1 / sqrt(vals)

  # x = ctr + T x'  turns  Gx >= H  into  (G T) x' >= H - G ctr.
  Gr <- G %*% Tm
  Hr <- as.numeric(H - as.numeric(G %*% ctr))

  as_mat <- function(z) if (is.matrix(z)) z else matrix(as.numeric(z), nrow = 1L)
  list(
    G = Gr, H = Hr, center = ctr, T = Tm, Tinv = Tinv,
    forth = function(x)  { M <- as_mat(x); R <- sweep(M, 2L, ctr, "-") %*% t(Tinv)
                           if (is.matrix(x)) R else as.numeric(R) },
    back  = function(xp) { M <- as_mat(xp); R <- sweep(M %*% t(Tm), 2L, ctr, "+")
                           if (is.matrix(xp)) R else as.numeric(R) },
    axis.ratio = max(semi) / min(semi),
    log.volume = sum(log(semi)),
    center.type = center
  )
}


#' @param lim A list with four components \code{A}, \code{B}, \code{G} and \code{H}
#'   representing the polytope, already reduced by \code{lim.redpol()}.
#' @export
#' @rdname pol.round
lim.round <- function(lim, center = c("analytic", "chebyshev"), x0 = NULL) {
  pol.round(G = lim$G, H = lim$H, center = center, x0 = x0)
}
