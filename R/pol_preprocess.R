

#' Remove redundant inequality constraints of a polytope
#'
#' The functions \code{pol.exfoliate()} and \code{lim.exfoliate()} remove the redundant
#' inequality constraints of a polytope
#' \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Gx \geq H \}}.
#'
#' An inequality constraint is \emph{redundant} when it is implied by the others:
#' removing it leaves the polytope \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}}
#' unchanged. The function \code{pol.exfoliate()} detects and removes all of them,
#' yielding a description from which no further constraint can be dropped.
#'
#' Constraint \eqn{i} is redundant if and only if
#' \eqn{\min \{ \langle g_i, x \rangle : \langle g_j, x \rangle \geq H_j, j \neq i \} \geq H_i},
#' which is one linear program per constraint.
#'
#' The sweep is \strong{sequential}: each constraint is tested against the constraints
#' still retained, not against the full set. This matters. Two identical constraints are
#' each redundant \emph{in the presence of the other}; testing both against the full set
#' would drop both, and change the polytope. The price is that the outcome depends on the
#' sweep order, so the minimal description is not unique in general. The polytope is.
#'
#' Removing redundant constraints leaves the polytope unchanged but modifies everything
#' that depends on its \emph{description} rather than on its geometry, in particular the
#' analytic centre and the Dikin ellipsoid, which are defined by sums over the constraints.
#' It should therefore be applied \strong{before} \code{pol.center()} and
#' \code{pol.round()}. On the BOWF-short model, 28 of the 72 inequality constraints of the
#' reduced polytope are redundant.
#'
#' @param G A matrix corresponding to \code{G} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param H A numeric vector corresponding to \code{H} in the description of the polytope
#'   \eqn{\mathcal{P}}.
#' @param tol A positive tolerance. Constraint \eqn{i} is declared redundant when the
#'   attained minimum is at least \eqn{H_i - tol}. Rows are normalised internally, so that
#'   \code{tol} is comparable across constraints.
#'
#' @return For \code{pol.exfoliate()}, a list with three components:
#'   \item{G}{the matrix \code{G} deprived of its redundant rows;}
#'   \item{H}{the corresponding right-hand side;}
#'   \item{redundant}{a logical vector of length \code{nrow(G)}, \code{TRUE} for the rows
#'     that were removed.}
#'   For \code{lim.exfoliate()}, the input list with \code{G} and \code{H} replaced by
#'   their exfoliated counterparts, every other component left untouched.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
#'
#' @rdname pol.exfoliate
#' @examples
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' red <- lim.redpol(BOWF)
#' exf <- pol.exfoliate(G = red$G, H = red$H)
#' # 28 of the 72 inequality constraints are redundant
#' sum(exf$redundant)
pol.exfoliate <- function(G, H, tol = 1e-9) {
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(G)) G <- t(G)
  if (is.null(G)) stop("G is NULL, the polytope has 0 dimensions.")
  H <- as.numeric(H)
  if (nrow(G) != length(H)) stop("G and H have incompatible dimensions.")
  if (any(!is.finite(G)) || any(!is.finite(H))) stop("G and H must be finite.")

  # Internally we work with Gx >= H rewritten as (-G) x <= (-H), rows normalised so
  # that `tol` means the same thing for every constraint.
  A <- -G; b <- -H
  s <- sqrt(rowSums(A^2))
  if (any(s <= 0)) stop("Degenerate constraint: some row of G is zero.")
  A <- A / s; b <- b / s
  m <- nrow(A); d <- ncol(A)

  keep <- rep(TRUE, m)
  # Variables are free: without this, Rglpk would silently impose x >= 0 and the
  # redundancy test would be run on the wrong domain.
  bnd <- list(lower = list(ind = seq_len(d), val = rep(-Inf, d)))

  for (i in seq_len(m)) {
    idx <- setdiff(which(keep), i)
    if (!length(idx)) next
    sol <- tryCatch(Rglpk_solve_LP(obj = A[i, ], mat = A[idx, , drop = FALSE],
                                   dir = rep("<=", length(idx)), rhs = b[idx],
                                   bounds = bnd, max = TRUE),
                    error = function(e) NULL)
    # A non-zero status means unbounded or failed. Unbounded means that, without
    # constraint i, direction g_i is no longer bounded: i is therefore supporting.
    # On failure we also keep it: keeping one constraint too many is harmless,
    # removing one too many changes the polytope.
    if (!is.null(sol) && identical(as.integer(sol$status), 0L) &&
        sol$optimum <= b[i] + tol) {
      keep[i] <- FALSE
    }
  }
  list(G = G[keep, , drop = FALSE], H = H[keep], redundant = !keep)
}


#' @param lim A list describing the \strong{reduced} polytope, with at least the
#'   components \code{G} and \code{H}, as returned by \code{\link{lim.redpol}}. An
#'   object still carrying equality constraints in \code{lim$A} is rejected: see
#'   \sQuote{Details}.
#' @export
#' @rdname pol.exfoliate
lim.exfoliate <- function(lim, tol = 1e-9) {
  .lim_check_reduced(lim)
  exf <- pol.exfoliate(G = lim$G, H = lim$H, tol = tol)
  lim$G <- exf$G
  lim$H <- exf$H
  lim
}


# The three pre-processing steps only make sense on the REDUCED polytope, the one
# lim.redpol() returns, whose description is purely {x : Gx >= H}. Fed a full lim
# object, they would silently ignore the equality constraints Ax = B and return a
# point outside the model -- observed at 1e8 away on BOWF-short. Hence this guard:
# better a clear refusal than a wrong answer.
.lim_check_reduced <- function(lim) {
  if (!is.list(lim) || is.null(lim$G) || is.null(lim$H))
    stop("`lim` must be a list with components G and H.", call. = FALSE)
  if (!is.null(lim$A) && length(lim$A) > 0L)
    stop("Equality constraints found in `lim$A`. These functions operate on the ",
         "REDUCED polytope {x : Gx >= H}: apply lim.redpol() first.", call. = FALSE)
  invisible(TRUE)
}
