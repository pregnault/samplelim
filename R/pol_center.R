

#' Centres of a polytope
#'
#' Several points claim the title of "centre" of a polytope
#' \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}}. They do not coincide,
#' and they do not have the same properties. \code{pol.center()} computes either of two.
#'
#' Writing \eqn{d_i(x) = \langle g_i, x \rangle - H_i \geq 0} for the margin of
#' constraint \eqn{i}:
#'
#' \describe{
#'   \item{\code{"chebyshev"}}{the centre of the largest inscribed ball, obtained by one
#'     linear program. It only looks at the \emph{nearest} faces and ignores all the
#'     others; on an elongated body it may sit anywhere along the long axis, and it need
#'     not even be unique (a rectangle has a whole segment of them).}
#'   \item{\code{"analytic"}}{the minimiser of the logarithmic barrier
#'     \eqn{F(x) = - \sum_i \log d_i(x)}, that is the point maximising the product of the
#'     margins. It takes \emph{all} the faces into account, and it is unique.}
#' }
#'
#' The difference is not cosmetic. On the BOWF-short model the two centres lie 1024 apart,
#' for a Chebyshev radius of 0.868 — and rounding at the Chebyshev centre rather than the
#' analytic one costs a factor of 100 to 200 in effective sample size per second (see
#' \code{\link{pol.round}}).
#'
#' The analytic centre is defined by a \emph{sum over the constraints}: it therefore
#' depends on the \emph{list} describing the polytope, not only on its geometry. Adding a
#' redundant constraint moves it. Apply \code{\link{pol.exfoliate}} first.
#'
#' The analytic centre is computed by damped Newton iterations. The step is bounded to 95%
#' of the distance to the boundary, so that the iterate stays strictly interior, and an
#' Armijo line search guarantees an actual decrease. The stopping rule uses the Newton
#' decrement \eqn{\lambda^2 = \nabla F^\top (\nabla^2 F)^{-1} \nabla F}, whose half bounds
#' \eqn{F(x) - F(c)} and which is invariant under affine changes of coordinates, unlike
#' the norm of the gradient.
#'
#' Note that the criterion bounds the gap on the \emph{objective}, not on the
#' \emph{position}. On an ill-conditioned body the barrier is very flat along the long
#' directions, so two correct implementations may differ noticeably in position while both
#' having a tiny gradient. Do not test convergence on the position.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param type Either \code{"analytic"} (the default) or \code{"chebyshev"}.
#' @param x0 A strictly interior starting point for the Newton iterations. If \code{NULL},
#'   the Chebyshev centre is used.
#' @param max_iter Maximum number of Newton iterations. A safeguard; never reached in
#'   practice.
#' @param tol Threshold on half the Newton decrement.
#'
#' @return For \code{type = "analytic"}, a numeric vector of length \eqn{n}. For
#'   \code{type = "chebyshev"}, a list with components \code{center} and \code{radius}.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
#'
#' @rdname pol.center
#' @examples
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' red <- lim.redpol(BOWF)
#' exf <- pol.exfoliate(G = red$G, H = red$H)
#' ctr <- pol.center(G = exf$G, H = exf$H, type = "analytic")
#' # the gradient of the barrier vanishes there
#' max(abs(colSums(exf$G / as.numeric(exf$G %*% ctr - exf$H))))
pol.center <- function(G, H, type = c("analytic", "chebyshev"),
                       x0 = NULL, max_iter = 200L, tol = 1e-10) {
  type <- match.arg(type)
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(G)) G <- t(G)
  if (is.null(G)) stop("G is NULL, the polytope has 0 dimensions.")
  H <- as.numeric(H)
  if (nrow(G) != length(H)) stop("G and H have incompatible dimensions.")

  # Gx >= H is rewritten as (-G) x <= (-H) for the internal computation.
  A <- -G; b <- -H

  if (identical(type, "chebyshev")) return(.pol_chebyshev(A, b))
  if (is.null(x0)) x0 <- .pol_chebyshev(A, b)$center
  x <- as.numeric(x0)
  if (any(as.numeric(A %*% x) >= b)) stop("x0 must be strictly interior (Gx > H).")

  converged <- FALSE
  for (it in seq_len(max_iter)) {
    d <- b - as.numeric(A %*% x)                 # margins
    if (any(d <= 0)) stop("Analytic centre: the current point is not interior.")
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
      if (t < 1e-14) stop("Analytic centre: line search stalled.")
    }
    x <- cand
    if (t * max(abs(step)) < tol) { converged <- TRUE; break }
  }
  if (!converged) warning("Analytic centre: maximum number of iterations reached.")
  x
}


#' @param lim A list with four components \code{A}, \code{B}, \code{G} and \code{H}
#'   representing the polytope.
#' @export
#' @rdname pol.center
lim.center <- function(lim, type = c("analytic", "chebyshev"), x0 = NULL,
                       max_iter = 200L, tol = 1e-10) {
  pol.center(G = lim$G, H = lim$H, type = type, x0 = x0,
             max_iter = max_iter, tol = tol)
}


# Chebyshev centre of {x : A x <= b}, by one linear program:
#   maximise r subject to  a_i . x + ||a_i|| r <= b_i,
# which states that the ball of centre x and radius r fits in half-space i.
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
