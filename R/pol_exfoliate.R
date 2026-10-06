#' Remove the redundant constraints of a polytope
#'
#' The functions \code{pol.exfoliate()} and \code{lim.exfoliate()} remove the redundant
#' inequality constraints of a given polytope
#' \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Gx \geq H \}}, that is, the constraints implied
#' by the others. The polytope itself is unchanged.
#'
#' Constraint \eqn{i} is redundant if the minimum of \eqn{\langle g_i, x \rangle} under the
#' other constraints is at least \eqn{H_i}, \eqn{g_i} being the \eqn{i}-th row of \code{G}:
#' one linear program per constraint. The constraints are tested in turn against those
#' still kept, so that of two identical constraints only one is removed. A zero row of
#' \code{G}, as left by \code{\link{lim.redpol}()} for an inequality constant on
#' \eqn{\{Ax = B\}}, is redundant.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope \eqn{\mathcal{P}}.
#' @param tol A numeric value specifying the tolerance of the test, the rows of \code{G}
#'   being normalised.
#'
#' @return For \code{pol.exfoliate()}, a list with three components; namely:
#' \itemize{
#' \item \code{G} and \code{H}, without the redundant constraints;
#' \item \code{redundant}, a logical vector, \code{TRUE} for the removed constraints.
#' }
#' For \code{lim.exfoliate()}, \code{lim} with its components \code{G} and \code{H} replaced.
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
#'
#' @rdname pol.exfoliate
#' @references J. Telgen,
#' \emph{Identifying redundant constraints and implicit equalities in systems of linear
#' constraints},
#' Management Science \strong{29(10)}, 1209-1222 (1983).
pol.exfoliate <- function(G, H, tol = 1e-9) {
  P <- .pol_check(G, H); G <- P$G; H <- P$H
  s <- sqrt(rowSums(G^2)); Gn <- G / s; Hn <- H / s
  d <- ncol(G)
  free <- list(lower = list(ind = seq_len(d), val = rep(-Inf, d)))
  keep <- s > 0                                  # a zero row always holds
  for (i in which(keep)) {
    idx <- setdiff(which(keep), i)
    if (!length(idx)) next
    sol <- Rglpk_solve_LP(obj = Gn[i, ], mat = Gn[idx, , drop = FALSE],
                          dir = rep(">=", length(idx)), rhs = Hn[idx], bounds = free)
    # A failed or unbounded program keeps the constraint: keeping one too many is
    # harmless, removing one too many would change the polytope.
    if (sol$status == 0L && sol$optimum >= Hn[i] - tol) keep[i] <- FALSE
  }
  list(G = G[keep, , drop = FALSE], H = H[keep], redundant = !keep)
}


#' @param lim A list with components \code{G} and \code{H} describing a reduced polytope,
#'   as returned by \code{\link{lim.redpol}()}. A list with equality constraints
#'   (component \code{A}) is refused.
#' @export
#' @rdname pol.exfoliate
lim.exfoliate <- function(lim, tol = 1e-9) {
  .lim_check_reduced(lim)
  exf <- pol.exfoliate(G = lim$G, H = lim$H, tol = tol)
  lim$G <- exf$G
  lim$H <- exf$H
  lim
}
