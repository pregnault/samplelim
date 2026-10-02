#' Centers of a polytope
#'
#' The functions \code{pol.center()} and \code{lim.center()} compute a center of a given
#' polytope \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}}, the equality
#' constraints being optional. Two notions of center are available; they do not coincide
#' and do not share the same properties.
#'
#' Writing \eqn{d_i(x) = \langle g_i, x \rangle - H_i \geq 0} for the margin of
#' inequality constraint \eqn{i}, the two centers are the following.
#' \describe{
#'   \item{\code{"chebyshev"}}{the center of a largest ball inscribed in
#'     \eqn{\mathcal{P}}, the ball being drawn in the affine space \eqn{\{Ax = B\}};
#'     obtained by a single linear program. It involves only the \emph{nearest} faces.
#'     The radius is unique, the center is not in general: the set of centers may be
#'     large, and the point returned depends on the solver and on the formulation. This
#'     set depends on the polytope only, not on its description.}
#'   \item{\code{"analytic"}}{the minimiser of the logarithmic barrier
#'     \eqn{F(x) = - \sum_i \log d_i(x)} on the affine space, equivalently the point
#'     maximising the product of the margins. It involves \emph{all} the faces. It is
#'     unique for a given description of \eqn{\mathcal{P}}, but depends on this
#'     description: a redundant constraint moves it, and a duplicated one counts twice.
#'     The function \code{\link{pol.exfoliate}()} should be applied first.}
#' }
#'
#' With equality constraints, the analytic center is computed on the reduced polytope
#' returned by \code{\link{lim.redpol}()}, then brought back: the barrier being invariant
#' under affine changes of coordinates, it is the same point. The inequalities that are
#' constant on the affine hull of \eqn{\mathcal{P}} are left out of the barrier.
#'
#' The two centers may lie far apart. On the BOWF-short model, reduced and exfoliated,
#' the distance between them is of the order of a thousand (874 for the Chebyshev
#' center returned by the linear program), for a Chebyshev radius of 0.868. Rounding at
#' the Chebyshev center instead of the analytic one costs a factor 100 to 200 in
#' effective sample size per second; see \code{\link{pol.round}()}.
#'
#' The analytic center is computed by damped Newton iterations. The step is bounded to
#' 95\% of the distance to the boundary, so that the iterate remains strictly interior,
#' and an Armijo line search ensures an actual decrease of \eqn{F}. The stopping rule
#' relies on the Newton decrement
#' \eqn{\lambda^2 = \nabla F^\top (\nabla^2 F)^{-1} \nabla F}, half of which bounds
#' \eqn{F(x) - F(c)} and which is invariant under affine changes of coordinates, unlike
#' the norm of the gradient. It bounds the gap on the \emph{objective}, not on the
#' \emph{position}: on an ill-conditioned polytope, two correct implementations may
#' return noticeably different points.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param type A character string specifying the notion of center to be computed, whether
#'   \code{"analytic"} (the default) or \code{"chebyshev"}; see the section
#'   \emph{Details} above.
#' @param x0 A numeric vector giving the coordinates of a strictly interior point, used
#'   as starting point for the Newton iterations; with equality constraints, it must
#'   satisfy them. If \code{NULL} (the default), a Chebyshev center is used.
#' @param max_iter An integer giving the maximum number of Newton iterations. It is a
#'   safeguard, never reached in practice.
#' @param tol A numeric value specifying the threshold on half the Newton decrement.
#' @param A,B Optional matrix and numeric vector corresponding to \code{A} and \code{B}
#'   in the description of the polytope \eqn{\mathcal{P}}. If \code{NULL} (the default),
#'   the polytope is \eqn{\{ x : Gx \geq H \}}.
#'
#' @return For \code{type = "analytic"}, a numeric vector of length \eqn{n}, the number
#' of columns of \code{G}. For \code{type = "chebyshev"}, a list with two components;
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
#' # A Chebyshev center of the full polytope, in the space of the flows
#' ctr <- lim.center(BOWF, type = "chebyshev")
#' ctr$radius
#' # The analytic center is better computed once the redundant constraints are removed
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
                       x0 = NULL, max_iter = 200L, tol = 1e-10, A = NULL, B = NULL) {
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
  if (!is.null(A)) return(.pol_center_full(A, B, G, H, type, x0, max_iter, tol))

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


#' @param lim A list with at least two components \code{G} and \code{H}, and possibly
#'   \code{A} and \code{B}: either a full \code{lim} object, as returned by
#'   \code{\link{df2lim}()}, or a reduced polytope, as returned by
#'   \code{\link{lim.redpol}()}. The center is given in the coordinates of \code{lim}.
#' @export
#' @rdname pol.center
lim.center <- function(lim, type = c("analytic", "chebyshev"), x0 = NULL,
                       max_iter = 200L, tol = 1e-10) {
  if (!is.list(lim) || is.null(lim$G) || is.null(lim$H))
    stop("`lim` must be a list with components G and H.", call. = FALSE)
  full <- !is.null(lim$A) && length(lim$A) > 0L
  pol.center(G = lim$G, H = lim$H, type = type, x0 = x0, max_iter = max_iter,
             tol = tol, A = if (full) lim$A, B = if (full) lim$B)
}


# Center of {x : Ax = B, Gx >= H}. Chebyshev: one linear program in the full space, the
# ball being drawn in the affine hull of the polytope (see .chebyshev_hull() in
# redpol.R). Analytic: on the affine space, the barrier is that of the reduced
# polytope, so the center is computed there, from the origin (a Chebyshev center of the
# reduced polytope), and brought back.
.pol_center_full <- function(A, B, G, H, type, x0, max_iter, tol) {
  if (is.data.frame(A)) A <- as.matrix(A)
  if (is.vector(A)) A <- t(A)
  B <- as.numeric(B)
  if (ncol(A) != ncol(G) || nrow(A) != length(B))
    stop("A, B and G have incompatible dimensions.")
  if (any(!is.finite(A)) || any(!is.finite(B))) stop("A and B must be finite.")
  eps <- sqrt(.Machine$double.eps)

  if (identical(type, "chebyshev")) {
    ctr <- .chebyshev_hull(A, B, G, H, eps)
    return(list(center = ctr$center, radius = ctr$radius))
  }
  red <- lim.redpol(list(A = A, B = B, G = G, H = H))
  z0 <- numeric(ncol(red$Z))
  if (!is.null(x0)) {
    x0 <- as.numeric(x0)
    if (max(abs(A %*% x0 - B)) > eps * max(1, abs(B)))
      stop("x0 must satisfy the equality constraints (Ax = B).")
    z0 <- as.numeric(crossprod(red$Z, x0 - red$x0))
  }
  z <- pol.center(red$G, red$H, type = "analytic", x0 = z0, max_iter = max_iter,
                  tol = tol)
  as.numeric(red$x0 + red$Z %*% z)
}


# Chebyshev center of {x : A x <= b}, by one linear program:
#   maximise r subject to  a_i . x + ||a_i|| r <= b_i,
# which states that the ball of center x and radius r fits in half-space i.
.pol_chebyshev <- function(A, b) {
  s <- sqrt(rowSums(A^2)); d <- ncol(A)
  # A polytope smaller than 1 is scaled up to unit size; a larger one is left as it is,
  # see .chebyshev_hull() in redpol.R.
  sc <- max(abs(b)); if (!is.finite(sc) || sc == 0) sc <- 1
  sc <- min(1, sc)
  sol <- Rglpk_solve_LP(
    obj = c(numeric(d), 1), mat = cbind(A, s),
    dir = rep("<=", nrow(A)), rhs = b / sc,
    bounds = list(lower = list(ind = seq_len(d + 1L), val = c(rep(-Inf, d), 0))),
    max = TRUE)
  if (sol$status != 0L || sol$solution[d + 1L] <= 0)
    stop("No interior point found: is the polytope empty or unbounded?")
  list(center = sol$solution[seq_len(d)] * sc, radius = sol$solution[d + 1L] * sc)
}
