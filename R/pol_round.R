#' Rounding a polytope by means of the Dikin ellipsoid
#'
#' The functions \code{pol.round()} and \code{lim.round()} round a given polytope
#' \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Gx \geq H \}}, that is, apply to it an affine
#' change of coordinates making it as isotropic as possible.
#'
#' A strongly elongated polytope is difficult to explore for a random walk: a billiard
#' trajectory reflects without progressing, and an isotropic Gaussian jump is accepted
#' along the short directions and rejected along the long ones. Rounding removes this
#' anisotropy.
#'
#' The Dikin ellipsoid at an interior point \eqn{c} is
#' \eqn{E(c) = \{x : (x-c)^\top \nabla^2 F(c) (x-c) \leq 1\}}, with
#' \eqn{\nabla^2 F(c) = \sum_i g_i g_i^\top / d_i(c)^2} the Hessian of the logarithmic
#' barrier, \eqn{g_i} the \eqn{i}-th row of \code{G} and
#' \eqn{d_i(c) = \langle g_i, c \rangle - H_i}. It is inscribed in \eqn{\mathcal{P}} for
#' \emph{any} interior point \eqn{c}. Applying \eqn{x = c + T x'} with
#' \eqn{T = \nabla^2 F(c)^{-1/2}} maps it onto the unit ball, hence rounds the polytope.
#'
#' This transformation is affine, so that its Jacobian \eqn{|\det T|} is constant. The
#' image of the uniform distribution on the rounded polytope is therefore exactly the
#' uniform distribution on the original one: it suffices to sample in the rounded space
#' and to apply \code{back}, with neither correction nor reweighting. This would not hold for
#' a non-uniform target distribution, whose density would have to be transported as well.
#'
#' The choice of the center is of practical importance. Since \eqn{E(c)} is inscribed for
#' any interior \eqn{c}, one may round from the Chebyshev center, which is readily
#' available; see \code{\link{pol.center}()}. Measured on the exfoliated BOWF-short
#' model, in minimum effective sample size per second:
#'
#' \tabular{lrr}{
#'   \tab BiW \tab MiW \cr
#'   no rounding \tab 51 \tab 77 \cr
#'   rounded at the Chebyshev center \tab 74 \tab 125 \cr
#'   rounded at the analytic center \tab 14353 \tab 12262
#' }
#'
#' Rounding from the Chebyshev center yields a factor 1.5 to 1.6, and from the analytic
#' center a factor 160 to 280. The reason is geometric: the Chebyshev center involves the
#' nearest faces only, so that the ellipsoid taken there is small and nearly spherical,
#' and fits the polytope loosely. The ellipsoid at the analytic center is elongated in
#' the same way as the polytope, which is what allows it to straighten it.
#'
#' The function \code{\link{pol.exfoliate}()} should be applied first, since
#' \eqn{\nabla^2 F} is a sum over the inequality constraints: a redundant constraint
#' adds a term that distorts the ellipsoid for no geometric reason. On the BOWF-short
#' model, exfoliating first yields a factor 55 in ellipsoid volume.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope
#'   \eqn{\mathcal{P}}. Redundant constraints should be removed beforehand.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param center A character string specifying the center at which the Dikin ellipsoid is
#'   taken, whether \code{"analytic"} (the default) or \code{"chebyshev"}. The latter is
#'   provided so that the comparison reported in the section \emph{Details} below can be
#'   reproduced, and is not recommended for actual sampling.
#' @param x0 A numeric vector giving the coordinates of a strictly interior point, passed
#'   to \code{\link{pol.center}()} as starting point of the Newton iterations. It is
#'   unrelated to the component \code{x0} returned by \code{\link{lim.redpol}()}.
#'
#' @return A list with ten components; namely:
#' \itemize{
#'   \item \code{G} and \code{H}, describing the rounded polytope
#'   \eqn{\{x' : G' x' \geq H'\}}, with \eqn{G' = G T} and \eqn{H' = H - G c}. Their
#'   rows are \emph{not} normalised;
#'   \item \code{center}, the center used, in the original coordinates;
#'   \item \code{T} and \code{Tinv}, the transformation and its inverse, so that
#'   \eqn{x = c + T x'}, \eqn{c} being \code{center};
#'   \item \code{forth}, a function mapping original coordinates to rounded ones;
#'   \item \code{back}, a function mapping rounded coordinates back to the original
#'   space. This is the one to be applied to a sample. It accepts a numeric vector or a
#'   matrix with one point per row, that is, the output of \code{\link{rlim}()} as is;
#'   \item \code{axis.ratio}, the ratio of the semi-axes of the Dikin ellipsoid. It is
#'   descriptive only, as it measures the shape of the ellipsoid and not the quality of
#'   the rounding; the latter is best assessed with \code{\link{pol.ranges}()} on the
#'   rounded polytope;
#'   \item \code{log.volume}, the log-volume of the ellipsoid, up to an additive
#'   constant. This one does measure quality: the larger, the more closely the ellipsoid
#'   fits the polytope;
#'   \item \code{center.type}, the value of the argument \env{center}.
#' }
#'
#' @export
#'
#' @rdname pol.round
#' @examples
#' # Create a lim object from a Description file
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' # These functions operate on the reduced polytope, exfoliated beforehand
#' red <- lim.redpol(BOWF)
#' exf <- pol.exfoliate(G = red$G, H = red$H)
#' rnd <- pol.round(G = exf$G, H = exf$H)
#'
#' # Anisotropy, before and after rounding
#' rg0 <- pol.ranges(G = exf$G, H = exf$H)[, 3]
#' rg1 <- pol.ranges(G = rnd$G, H = rnd$H)[, 3]
#' c(before = max(rg0) / min(rg0), after = max(rg1) / min(rg1))
#'
#' # Sample in the rounded space, then come back. Hpolytope uses the convention
#' # Ax <= b, hence the change of sign.
#' smp <- rlim(lim = NULL, Hpol = Hpolytope(A = -rnd$G, b = -rnd$H), nsamp = 100,
#'             seed = 123)
#' pts <- rnd$back(smp)
#' @seealso \code{\link{lim.redpol}()} for reducing a polytope,
#' \code{\link{pol.exfoliate}()} for removing its redundant constraints,
#' \code{\link{pol.center}()} for its centers, \code{\link{rlim}()} for sampling it.
#' @references {
#' I. I. Dikin,
#' \emph{Iterative solution of problems of linear and quadratic programming},
#' Soviet Mathematics Doklady \strong{8}, 674-675 (1967).
#'
#' S. Boyd and L. Vandenberghe,
#' \emph{Convex Optimization},
#' Cambridge University Press (2004).
#' }
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
  if (any(dc <= 0)) stop("The chosen center is not interior.")

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


#' @param lim A list with at least two components \code{G} and \code{H} representing
#'   the \emph{reduced} polytope, as returned by \code{\link{lim.redpol}()}. A list
#'   still carrying equality constraints in its component \code{A} is rejected.
#' @export
#' @rdname pol.round
lim.round <- function(lim, center = c("analytic", "chebyshev"), x0 = NULL) {
  .lim_check_reduced(lim)
  pol.round(G = lim$G, H = lim$H, center = center, x0 = x0)
}
