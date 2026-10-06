#' Centers of a polytope
#'
#' The functions \code{pol.center()} and \code{lim.center()} compute the Chebyshev center
#' or the analytic center of a given polytope
#' \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}}, the equality constraints
#' being optional.
#'
#' The Chebyshev center is the center of the largest ball inscribed in \eqn{\mathcal{P}},
#' the ball being taken in the affine space \eqn{\{Ax = B\}}. It is computed by one linear
#' program, and need not be unique.
#'
#' The analytic center maximises the product of the margins
#' \eqn{\langle g_i, x \rangle - H_i}, \eqn{g_i} being the \eqn{i}-th row of \code{G}. It
#' is unique, but depends on the description of the polytope: a redundant constraint moves
#' it, hence \code{\link{pol.exfoliate}()} should be applied first. It is computed by damped
#' Newton iterations, in the coordinates of the affine space \eqn{\{Ax = B\}}.
#'
#' @param A A matrix corresponding to \code{A} in the description of the polytope
#'   \eqn{\mathcal{P}}, or \code{NULL} (the default) if there are no equality constraints.
#' @param B A numeric vector corresponding to \code{B} in the description of the polytope
#'   \eqn{\mathcal{P}}, or \code{NULL} (the default).
#' @param G A matrix corresponding to \code{G} in the description of the polytope \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope \eqn{\mathcal{P}}.
#' @param type A character string specifying the center, whether \code{"analytic"} (the
#'   default) or \code{"chebyshev"}.
#' @param x0 A numeric vector giving the coordinates of a strictly interior point, used as
#'   starting point of the Newton iterations. If \code{NULL} (the default), the Chebyshev
#'   center.
#' @param max_iter An integer giving the maximum number of Newton iterations.
#' @param tol A numeric value specifying the tolerance on half the squared Newton
#'   decrement, which stops the iterations.
#'
#' @return For \code{type = "analytic"}, a numeric vector of length \eqn{n}. For
#' \code{type = "chebyshev"}, a list with two components; namely:
#' \itemize{
#' \item \code{center}, a numeric vector of length \eqn{n};
#' \item \code{radius}, the radius of the inscribed ball.
#' }
#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom MASS Null
#' @export
#'
#' @rdname pol.center
#' @references S. Boyd and L. Vandenberghe,
#' \emph{Convex Optimization},
#' Cambridge University Press (2004), Section 8.5.
pol.center <- function(A = NULL, B = NULL, G, H, type = c("analytic", "chebyshev"),
                       x0 = NULL, max_iter = 200L, tol = 1e-10) {
  type <- match.arg(type)
  P <- .pol_check(G, H); G <- P$G; H <- P$H
  if (length(A) == 0L) A <- NULL
  if (!is.null(A)) {
    if (is.data.frame(A)) A <- as.matrix(A)
    if (is.vector(A)) A <- t(A)
    B <- as.numeric(B)
    if (ncol(A) != ncol(G) || nrow(A) != length(B))
      stop("A, B and G have incompatible dimensions.")
  }

  # x = c + Z z, the columns of Z being an orthonormal basis of the directions of
  # {Ax = B}; s_i is the norm of g_i within these directions, zero for an inequality
  # constant on {Ax = B}.
  Z <- if (is.null(A)) diag(ncol(G)) else Null(t(A))
  s <- sqrt(rowSums((G %*% Z)^2))
  s[s <= sqrt(.Machine$double.eps) * sqrt(rowSums(G^2))] <- 0
  if (identical(type, "chebyshev")) return(.pol_chebyshev(A, B, G, H, s))

  c0 <- if (is.null(x0)) .pol_chebyshev(A, B, G, H, s)$center else as.numeric(x0)
  if (!is.null(A) && max(abs(A %*% c0 - B)) > sqrt(.Machine$double.eps) * max(1, abs(B)))
    stop("x0 must satisfy the equality constraints (Ax = B).")
  Gz <- (G %*% Z)[s > 0, , drop = FALSE]
  Hz <- (H - as.numeric(G %*% c0))[s > 0]
  if (any(Hz >= 0)) stop("x0 must be strictly interior (Gx > H).")
  as.numeric(c0 + Z %*% .pol_analytic(Gz, Hz, max_iter, tol))
}


#' @param lim A list with components \code{G} and \code{H}, and possibly \code{A} and
#'   \code{B}: a \code{lim} object, as returned by \code{\link{df2lim}()}, or a reduced
#'   polytope, as returned by \code{\link{lim.redpol}()}.
#' @export
#' @rdname pol.center
lim.center <- function(lim, type = c("analytic", "chebyshev"), x0 = NULL,
                       max_iter = 200L, tol = 1e-10) {
  if (!is.list(lim) || is.null(lim$G) || is.null(lim$H))
    stop("`lim` must be a list with components G and H.", call. = FALSE)
  pol.center(A = lim$A, B = lim$B, G = lim$G, H = lim$H, type = type, x0 = x0,
             max_iter = max_iter, tol = tol)
}


# Chebyshev center, by one linear program: maximise r subject to Ax = B and
# g_i x - s_i r >= H_i, that is, the ball of center x and radius r within {Ax = B}
# lies in every half-space.
.pol_chebyshev <- function(A, B, G, H, s) {
  n <- ncol(G)
  mat <- cbind(G, -s); dir <- rep(">=", nrow(G)); rhs <- H
  if (!is.null(A)) {
    mat <- rbind(cbind(A, 0), mat); dir <- c(rep("==", nrow(A)), dir); rhs <- c(B, rhs)
  }
  sol <- Rglpk_solve_LP(obj = c(numeric(n), 1), mat = mat, dir = dir, rhs = rhs,
                        bounds = list(lower = list(ind = seq_len(n), val = rep(-Inf, n))),
                        max = TRUE)
  if (sol$status != 0L || sol$solution[n + 1L] <= sqrt(.Machine$double.eps))
    stop("No interior point found: the polytope is empty, unbounded, or flat ",
         "(some inequality holds as an equality on the whole polytope).")
  list(center = sol$solution[seq_len(n)], radius = sol$solution[n + 1L])
}


# Analytic center of {z : Gz >= H}, the origin being strictly interior: damped Newton
# iterations on F(z) = -sum(log(Gz - H)) from the origin, each step staying within the
# polytope and decreasing F (Armijo), until half the squared Newton decrement is
# below tol.
.pol_analytic <- function(G, H, max_iter, tol) {
  z <- numeric(ncol(G))
  for (it in seq_len(max_iter)) {
    d <- as.numeric(G %*% z) - H
    grad <- -as.numeric(crossprod(G, 1 / d))
    step <- -solve(crossprod(G, G / d^2), grad)
    dec2 <- -sum(grad * step)
    if (dec2 / 2 <= tol) return(z)
    Gs <- as.numeric(G %*% step)
    t <- if (any(Gs < 0)) min(1, 0.95 * min(-d[Gs < 0] / Gs[Gs < 0])) else 1
    repeat {
      dn <- as.numeric(G %*% (z + t * step)) - H
      if (all(dn > 0) && -sum(log(dn)) <= -sum(log(d)) - 0.01 * t * dec2) break
      t <- t / 2
      if (t < 1e-14) stop("Analytic center: the line search stalled.")
    }
    z <- z + t * step
  }
  warning("Analytic center: maximum number of iterations reached.")
  z
}
