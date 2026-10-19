#' Exfoliate a polytope: implicit equalities and redundant constraints
#'
#' The functions \code{pol.exfoliate()} and \code{lim.exfoliate()} clean the description of
#' a given polytope
#' \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}}, the equality constraints
#' being optional, without changing the polytope, up to the tolerances given below:
#' \enumerate{
#' \item the implicit equalities, that is, the inequality constraints holding as equalities
#' on the whole of \eqn{\mathcal{P}}, are moved from \code{G} to \code{A};
#' \item the redundant inequality constraints, that is, those implied by the others, are
#' removed.
#' }
#' This is the first step of the pre-processing, before the reduction by
#' \code{\link{lim.redpol}()} and the rounding by \code{\link{lim.round}()}.
#'
#' All the tests are made in orthonormal coordinates of the affine space
#' \eqn{\{Ax = B\}}, each row of \code{G} being normalised within this space, so that the
#' tolerances are distances. A constraint constant on \eqn{\{Ax = B\}} is compared with
#' \code{H} directly: an implicit equality if it holds as an equality, redundant if it
#' holds strictly, the polytope being empty otherwise.
#'
#' Constraint \eqn{i} is an implicit equality if its gap, the largest distance from a
#' point of \eqn{\mathcal{P}} to its hyperplane \eqn{\langle g_i, x \rangle = H_i}
#' within \eqn{\{Ax = B\}}, is zero, \eqn{g_i} being the \eqn{i}-th row of \code{G}.
#' There is none when \eqn{\mathcal{P}} contains a ball of \eqn{\{Ax = B\}} of radius
#' larger than \eqn{\varepsilon} = \code{2 * sqrt(.Machine$double.eps)}, twice the
#' threshold of \code{\link{pol.center}()}, which one linear program checks. Otherwise,
#' one linear program per constraint gives its gap, and the constraint of smallest gap,
#' if at most \eqn{2\varepsilon}, is moved to \code{A}, with the middle of the range of
#' \eqn{\langle g_i, x \rangle} over \eqn{\mathcal{P}} as value (\eqn{H_i} itself for a
#' gap at the level of rounding errors). The gaps of the other constraints change with
#' \eqn{\{Ax = B\}}: they are computed again, one implicit equality being moved at a
#' time, until \eqn{\mathcal{P}} contains such a ball. A polytope of Chebyshev radius at
#' most \eqn{\varepsilon} without a constraint of gap at most \eqn{2\varepsilon} is
#' refused.
#'
#' Constraint \eqn{i} is redundant if the minimum of \eqn{\langle g_i, x \rangle} under
#' \eqn{Ax = B} and the other inequality constraints is at least \eqn{H_i}, up to
#' \code{tol}: one linear program per constraint. The constraints are tested in turn
#' against those still kept, so that of two identical constraints only one is removed. A
#' linear program that fails keeps the constraint, with a warning.
#'
#' The tests being made within \eqn{\{Ax = B\}}, exfoliating \eqn{\mathcal{P}} and then
#' reducing it removes, up to the tolerances, the same inequality constraints as reducing
#' it and then exfoliating the reduced polytope. The tolerances are absolute, in the units
#' of the unknowns: a polytope much smaller than 1 in some direction should be rescaled.
#' Moving a constraint of nonzero gap to \code{A} keeps a slice of \eqn{\mathcal{P}}
#' only, and removing a constraint redundant up to \code{tol} may add to
#' \eqn{\mathcal{P}} a slice of thickness at most \code{tol}.
#'
#' @param A A matrix corresponding to \code{A} in the description of the polytope
#'   \eqn{\mathcal{P}}, or \code{NULL} (the default) if there are no equality constraints.
#' @param B A numeric vector corresponding to \code{B} in the description of the polytope
#'   \eqn{\mathcal{P}}, or \code{NULL} (the default).
#' @param G A matrix corresponding to \code{G} in the description of the polytope \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope \eqn{\mathcal{P}}.
#' @param tol A numeric value specifying the tolerance of the redundancy test, as a
#'   distance within \eqn{\{Ax = B\}}.
#'
#' @return For \code{pol.exfoliate()}, a list with six components; namely:
#' \itemize{
#' \item \code{A} and \code{B}, with the implicit equalities appended (\code{NULL} if there
#'   are no equality constraints at all);
#' \item \code{G} and \code{H}, without the implicit equalities and the redundant
#'   constraints;
#' \item \code{implicit}, a logical vector, \code{TRUE} for the implicit equalities, removed
#'   from \code{G} (and appended to \code{A} unless implied by it);
#' \item \code{redundant}, a logical vector, \code{TRUE} for the redundant constraints,
#'   removed from \code{G}.
#' }
#' For \code{lim.exfoliate()}, \code{lim} with its components \code{A}, \code{B}, \code{G}
#' and \code{H} replaced. Its other components, such as the names and numbers of the
#' declared constraints, are left unchanged.
#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom MASS Null ginv
#' @export
#'
#' @rdname pol.exfoliate
#' @references J. Telgen,
#' \emph{Identifying redundant constraints and implicit equalities in systems of linear
#' constraints},
#' Management Science \strong{29(10)}, 1209-1222 (1983).
pol.exfoliate <- function(A = NULL, B = NULL, G, H, tol = 1e-9) {
  P <- .pol_check(G, H); G <- P$G; H <- P$H
  E <- .pol_check_eq(A, B, G); A <- E$A; B <- E$B
  # Flatness threshold: twice that of pol.center(), so that an exfoliated polytope has,
  # in exact arithmetic, an interior point for pol.center().
  eps <- 2 * sqrt(.Machine$double.eps)
  implicit <- logical(nrow(G)); redundant <- logical(nrow(G))
  live <- rowSums(G^2) > 0                       # a zero row always holds
  redundant[!live] <- TRUE
  aff <- .pol_affine(A, B, ncol(G))              # {Ax = B} = {x0 + Z z}
  folds <- 0L                                    # implicit equalities appended so far

  # 1. Implicit equalities, one at a time, until the polytope contains a ball of {Ax = B}.
  repeat {
    rz <- .pol_rows(G, H, aff, live, slack = folds * eps)
    implicit[rz$tight] <- TRUE; redundant[rz$const & !rz$tight] <- TRUE
    live[rz$const] <- FALSE
    if (!any(live) || ncol(aff$Z) == 0L) break
    rad <- .pol_lp(c(numeric(ncol(aff$Z)), 1), cbind(rz$G[live, , drop = FALSE], -1),
                   rz$h[live], max = TRUE, free = ncol(aff$Z))
    if (rad$status == "infeasible")
      stop("The polytope is empty, or too small for the tolerances of the solver: rescale it.")
    if (rad$status == "unbounded") break
    if (rad$status != "optimal") stop("A linear program failed; the polytope may be too small: rescale it.")
    if (rad$optimum > eps) {
      # Origin moved to the Chebyshev center: the linear programs below start inside.
      aff$x0 <- aff$x0 + as.numeric(aff$Z %*% rad$solution[seq_len(ncol(aff$Z))])
      break
    }
    # Gap of each constraint, zero for an implicit equality; below the tolerance of the
    # solver, it may come out negative, hence the absolute value.
    top <- gap <- rep(Inf, nrow(G))
    for (i in which(live)) {
      hi <- .pol_lp(rz$G[i, ], rz$G[live, , drop = FALSE], rz$h[live], max = TRUE)
      if (hi$status == "optimal") { top[i] <- hi$optimum; gap[i] <- abs(hi$optimum - rz$h[i]) }
    }
    i <- which.min(gap)
    if (gap[i] > 2 * eps)
      stop("No interior point found: the Chebyshev radius of the polytope within {Ax = B} is ",
           signif(rad$optimum, 3), ", at most ", signif(eps, 2), ", and no inequality ",
           "constraint is within ", signif(2 * eps, 2), " of an implicit equality: rescale it.")
    # Constraint i appended to A, as g_i x = b_i, b_i the middle of the range of g_i x over
    # the polytope; exactly H_i for a gap at the level of rounding errors.
    if (gap[i] <= 1e3 * .Machine$double.eps * max(1, abs(rz$h[i]))) {
      b <- H[i]
    } else {
      lo <- .pol_lp(rz$G[i, ], rz$G[live, , drop = FALSE], rz$h[live])
      if (lo$status != "optimal") stop("A linear program failed; the polytope may be too small: rescale it.")
      b <- as.numeric(G[i, ] %*% aff$x0) + rz$sp[i] * (top[i] + lo$optimum) / 2
    }
    implicit[i] <- TRUE; live[i] <- FALSE
    u <- as.numeric(G[i, ] %*% aff$Z) / rz$sp[i]   # unit normal within {Ax = B}
    c0 <- (b - as.numeric(G[i, ] %*% aff$x0)) / rz$sp[i]
    A <- rbind(A, G[i, ]); B <- c(B, b); folds <- folds + 1L
    aff$x0 <- aff$x0 + as.numeric(aff$Z %*% (c0 * u))
    aff$Z <- aff$Z %*% Null(u)
  }

  # 2. Redundant constraints, tested in turn against those still kept, within {Ax = B}.
  rz <- .pol_rows(G, H, aff, live, slack = folds * eps)
  implicit[rz$tight] <- TRUE; redundant[rz$const & !rz$tight] <- TRUE
  live[rz$const] <- FALSE
  failed <- 0L
  for (i in which(live)) {
    idx <- setdiff(which(live), i)
    if (!length(idx)) next
    sol <- .pol_lp(rz$G[i, ], rz$G[idx, , drop = FALSE], rz$h[idx])
    if (sol$status == "optimal") {
      if (sol$optimum >= rz$h[i] - tol) { live[i] <- FALSE; redundant[i] <- TRUE }
    } else if (sol$status != "unbounded") failed <- failed + 1L
  }
  if (failed > 0L)
    warning(failed, " linear program(s) failed: the corresponding constraints are kept.")
  list(A = A, B = B, G = G[live, , drop = FALSE], H = H[live],
       implicit = implicit, redundant = redundant)
}


#' @param lim A list with components \code{G} and \code{H}, and possibly \code{A} and
#'   \code{B}: a \code{lim} object, as returned by \code{\link{df2lim}()}, or a reduced
#'   polytope, as returned by \code{\link{lim.redpol}()}.
#' @export
#' @rdname pol.exfoliate
lim.exfoliate <- function(lim, tol = 1e-9) {
  if (!is.list(lim) || is.null(lim$G) || is.null(lim$H))
    stop("`lim` must be a list with components G and H.", call. = FALSE)
  exf <- pol.exfoliate(A = lim$A, B = lim$B, G = lim$G, H = lim$H, tol = tol)
  if (any(exf$implicit) && length(lim$A) == 0L && !is.null(lim$Z))
    stop("The reduced polytope is flat: exfoliate the full polytope with lim.exfoliate(), ",
         "then reduce it.", call. = FALSE)
  lim$G <- exf$G
  lim$H <- exf$H
  if (NROW(exf$A) > NROW(lim$A)) {               # implicit equalities appended
    lim$A <- exf$A
    lim$B <- exf$B
  }
  lim
}


# Orthonormal coordinates of the affine space {Ax = B} of R^n: x = x0 + Z z, x0 the
# solution of least norm (A may have dependent rows), the columns of Z an orthonormal
# basis of the kernel of A; R^n itself if A is NULL.
.pol_affine <- function(A, B, n) {
  if (is.null(A)) return(list(x0 = numeric(n), Z = diag(n)))
  x0 <- as.numeric(ginv(A) %*% B)
  if (max(abs(A %*% x0 - B)) > sqrt(.Machine$double.eps) * max(1, abs(B)))
    stop("The equality constraints Ax = B have no solution.")
  list(x0 = x0, Z = Null(t(A)))
}

# Rows of G in the coordinates z, normalised by their norm within {Ax = B}: g z >= h. A
# row constant on {Ax = B} (const) holds as an equality (tight) or strictly, or the
# polytope is empty; `slack` is the distance by which the appended implicit equalities
# may have moved {Ax = B}.
.pol_rows <- function(G, H, aff, live, slack = 0) {
  Gz <- G %*% aff$Z
  hz <- H - as.numeric(G %*% aff$x0)
  sp <- sqrt(rowSums(Gz^2))
  s <- sqrt(rowSums(G^2))
  const <- live & sp <= 1e3 * .Machine$double.eps * s
  margin <- sqrt(.Machine$double.eps) * pmax(1, abs(H), s * sqrt(sum(aff$x0^2))) + slack * s
  if (any(hz[const] > margin[const]))
    stop("The polytope is empty, or too small for the tolerances: an inequality constant ",
         "on {Ax = B} does not hold (rescale it?).")
  tight <- const & hz >= -margin
  sp[const | sp == 0] <- 1
  list(G = Gz / sp, h = hz / sp, sp = sp, const = const, tight = tight)
}

# Linear program min (or max) obj'z subject to Mz >= rhs, the first `free` variables
# free, the others nonnegative; status "optimal", "unbounded", "infeasible" or "failed".
# The simplex may cycle on a degenerate polytope at the scale of the tolerances: it is
# stopped after 10 s, and the program counts as failed.
.pol_lp <- function(obj, M, rhs, max = FALSE, free = length(obj)) {
  sol <- Rglpk_solve_LP(obj = obj, mat = M, dir = rep(">=", nrow(M)), rhs = rhs,
                        bounds = list(lower = list(ind = seq_len(free), val = rep(-Inf, free))),
                        max = max, control = list(canonicalize_status = FALSE, tm_limit = 10000L))
  st <- switch(as.character(sol$status), "5" = "optimal", "6" = "unbounded",
               "4" = "infeasible", "failed")
  list(status = st, optimum = sol$optimum, solution = sol$solution)
}
