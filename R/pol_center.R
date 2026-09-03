#' Centers of a polytope
#'
#' The functions \code{pol.center()} and \code{lim.center()} compute a center of a given
#' polytope \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Gx \geq H \}}. Two notions of center
#' are available; they do not coincide and do not share the same properties.
#'
#' Writing \eqn{d_i(x) = \langle g_i, x \rangle - H_i \geq 0} for the margin of
#' inequality constraint \eqn{i}, the two centers are the following.
#' \describe{
#'   \item{\code{"chebyshev"}}{the center of the largest ball inscribed in
#'     \eqn{\mathcal{P}}, obtained by a single linear program. It involves only the
#'     \emph{nearest} faces and ignores all the others; on an elongated polytope, it may
#'     lie anywhere along the long axis, and it need not be unique, since a rectangle
#'     admits a whole segment of Chebyshev centers.}
#'   \item{\code{"analytic"}}{the minimiser of the logarithmic barrier
#'     \eqn{F(x) = - \sum_i \log d_i(x)}, equivalently the point maximising the product
#'     of the margins. It involves \emph{all} the faces, and it is unique.}
#' }
#'
#' The two centers may lie far apart. On the BOWF-short model, the distance between them
#' is 1024, for a Chebyshev radius of 0.868. Rounding at the Chebyshev center instead of
#' the analytic one costs a factor 100 to 200 in effective sample size per second;
#' see \code{\link{pol.round}()}.
#'
#' The analytic center is defined by a sum over the inequality constraints. It therefore
#' depends on the \emph{description} of the polytope, and not on its geometry only:
#' adding a redundant constraint moves it. The function \code{\link{pol.exfoliate}()}
#' should be applied first.
#'
#' The analytic center is computed by damped Newton iterations. The step is bounded to
#' 95\% of the distance to the boundary, so that the iterate remains strictly interior,
#' and an Armijo line search ensures an actual decrease of \eqn{F}. The stopping rule
#' relies on the Newton decrement
#' \eqn{\lambda^2 = \nabla F^\top (\nabla^2 F)^{-1} \nabla F}, half of which bounds
#' \eqn{F(x) - F(c)} and which is invariant under affine changes of coordinates, unlike
#' the norm of the gradient.
#'
#' This criterion bounds the gap on the \emph{objective}, not on the \emph{position}. On
#' an ill-conditioned polytope, the barrier is very flat along the long directions, so
#' that two correct implementations may return noticeably different points while both
#' having a small gradient. Convergence should not be assessed on the position.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param type A character string specifying the notion of center to be computed, whether
#'   \code{"analytic"} (the default) or \code{"chebyshev"}; see the section
#'   \emph{Details} above.
#' @param x0 A numeric vector giving the coordinates of a strictly interior point, used
#'   as starting point for the Newton iterations. If \code{NULL} (the default), the
#'   Chebyshev center is used.
#' @param max_iter An integer giving the maximum number of Newton iterations. It is a
#'   safeguard, never reached in practice.
#' @param tol A numeric value specifying the threshold on half the Newton decrement.
#'
#' @return For \code{type = "analytic"}, a numeric vector of length \eqn{n}, the
#' dimension of the polytope. For \code{type = "chebyshev"}, a list with two components;
#' namely:
#' \itemize{
#'   \item \code{center}
#'   \item \code{radius}
#' }
#'
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
#'
#' @rdname pol.center
#' @examples
#' # Create a lim object from a Description file
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' # These functions operate on the reduced polytope, exfoliated beforehand
#' red <- lim.redpol(BOWF)
#' exf <- pol.exfoliate(G = red$G, H = red$H)
#' ctr <- pol.center(G = exf$G, H = exf$H, type = "analytic")
#' # The gradient of the logarithmic barrier vanishes at the analytic center
#' max(abs(colSums(exf$G / as.numeric(exf$G %*% ctr - exf$H))))
#' @seealso \code{\link{lim.redpol}()} for reducing a polytope,
#' \code{\link{pol.exfoliate}()} for removing its redundant constraints,
#' \code{\link{pol.round}()} for rounding it.
#' @references {
#' S. Boyd and L. Vandenberghe,
#' \emph{Convex Optimization},
#' Cambridge University Press (2004).
#'
#' Y. Nesterov and A. Nemirovskii,
#' \emph{Interior-Point Polynomial Algorithms in Convex Programming},
#' SIAM (1994).
#' }
pol.center <- function(G, H, type = c("analytic", "chebyshev"),
                       x0 = NULL, max_iter = 200L, tol = 1e-10) {
  type <- match.arg(type)
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(G)) G <- t(G)
  if (is.null(G)) stop("G is NULL, the polytope has 0 dimensions.")
  H <- as.numeric(H)
  if (nrow(G) != length(H)) stop("G and H have incompatible dimensions.")
  # Same guards as pol.exfoliate() and pol.round(): without them a NA would surface
  # as an opaque "no interior point found" from the linear program.
  if (any(!is.finite(G)) || any(!is.finite(H))) stop("G and H must be finite.")
  if (any(sqrt(rowSums(G^2)) <= 0)) stop("Degenerate constraint: some row of G is zero.")

  # Gx >= H is rewritten as (-G) x <= (-H) for the internal computation.
  A <- -G; b <- -H

  if (identical(type, "chebyshev")) return(.pol_chebyshev(A, b))
  if (is.null(x0)) x0 <- .pol_chebyshev(A, b)$center
  x <- as.numeric(x0)
  if (any(as.numeric(A %*% x) >= b)) stop("x0 must be strictly interior (Gx > H).")

  converged <- FALSE
  for (it in seq_len(max_iter)) {
    d <- b - as.numeric(A %*% x)                 # margins
    if (any(d <= 0)) stop("Analytic center: the current point is not interior.")
    g <- as.numeric(crossprod(A, 1 / d))         # gradient  sum a_i / d_i
    Hess <- crossprod(A, A * (1 / d^2))          # Hessian   sum a_i a_i^T / d_i^2
    step <- tryCatch(solve(Hess, g),
                     error = function(e) solve(Hess + diag(1e-12, ncol(Hess)), g))
    decrement2 <- sum(g * step)                  # squared Newton decrement
    if (is.finite(decrement2) && decrement2 / 2 <= tol) { converged <- TRUE; break }

    # Longest step that stays interior: margin i becomes d_i + t * (A step)_i.
    Astep <- as.numeric(A %*% step)
    tmax <- 1.0
    neg <- Astep < 0
    if (any(neg)) tmax <- min(1.0, 0.95 * min(d[neg] / (-Astep[neg])))

    f0 <- -sum(log(d)); t <- tmax                # Armijo line search
    repeat {
      cand <- x - t * step
      dc <- b - as.numeric(A %*% cand)
      if (all(dc > 0) && -sum(log(dc)) <= f0 - 0.01 * t * decrement2) break
      t <- t / 2
      if (t < 1e-14) stop("Analytic center: line search stalled.")
    }
    x <- cand
    if (t * max(abs(step)) < tol) { converged <- TRUE; break }
  }
  if (!converged) warning("Analytic center: maximum number of iterations reached.")
  x
}


#' @param lim A list with at least two components \code{G} and \code{H} representing
#'   the \strong{reduced} polytope, as returned by \code{\link{lim.redpol}()}. A list
#'   still carrying equality constraints in its component \code{A} is rejected.
#' @export
#' @rdname pol.center
lim.center <- function(lim, type = c("analytic", "chebyshev"), x0 = NULL,
                       max_iter = 200L, tol = 1e-10) {
  .lim_check_reduced(lim)
  pol.center(G = lim$G, H = lim$H, type = type, x0 = x0,
             max_iter = max_iter, tol = tol)
}


# Chebyshev center of {x : A x <= b}, by one linear program:
#   maximise r subject to  a_i . x + ||a_i|| r <= b_i,
# which states that the ball of center x and radius r fits in half-space i.
.pol_chebyshev <- function(A, b) {
  s <- sqrt(rowSums(A^2)); d <- ncol(A)
  sol <- Rglpk_solve_LP(
    obj = c(numeric(d), 1), mat = cbind(A, s),
    dir = rep("<=", nrow(A)), rhs = b,
    bounds = list(lower = list(ind = seq_len(d + 1L), val = c(rep(-Inf, d), 0))),
    max = TRUE)
  if (sol$status != 0L || sol$solution[d + 1L] <= 0)
    stop("No interior point found: is the polytope empty or unbounded?")
  list(center = sol$solution[seq_len(d)], radius = sol$solution[d + 1L])
}
