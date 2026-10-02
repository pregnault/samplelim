#' Projection of full polytope into the reduced polytope
#'
#' The function \code{lim.redpol()} takes as input a polytope  and returns its projection into the non-empty reduced polytope.
#' Precisely, taking the polytope \eqn{\mathcal{P}= \{ x \in \mathbb{R}^n: Ax = B, Gx \geq H \}} as input, the function returns 
#' \itemize{
#' \item the matrix \code{Z}, orthonormal basis of the right null space of A, and \eqn{x_0} a point of \eqn{\mathcal{P}}, used for the reduction;
#' \item the matrix \eqn{G'=GZ} and the vector \eqn{H'=H-Gx_0} describing the reduced polytope \eqn{\mathcal{P'}= \{ x \in \mathbb{R}^{n-k}: G'x \geq H'\}} with \eqn{k=\mathtt{rank}(A)}, the rows flagged by \code{dropped} being removed.
#'  }
#'
#' The point \eqn{x_0} is a Chebyshev center of \eqn{\mathcal{P}}: the center of a largest
#' ball contained in \eqn{\mathcal{P}} within the affine space \eqn{\{Ax = B\}}, obtained
#' by a single linear program. It is strictly interior, so that the origin is a Chebyshev
#' center of \eqn{\mathcal{P'}}. It need not be unique.
#'
#' A ball of radius zero means that some inequalities hold as equalities on the whole
#' polytope. With \code{test = TRUE}, these implicit equalities are detected, by one linear
#' program per inequality, and added to the equalities, whose rank \eqn{k} increases
#' accordingly. The inequalities that are then
#' constant on the affine hull of \eqn{\mathcal{P}} give zero rows in \eqn{G'}; they are
#' removed, and flagged by the component \code{dropped}.
#'  
#'  The function \code{full2red()} (resp. \code{red2full()}) turns a sample of points inside the full polytope \eqn{\mathcal{P}} (resp. inside the reduced polytope \eqn{\mathcal{P'}})
#'  into the sample of corresponding points inside the reduced polytope \eqn{\mathcal{P'}} (resp. the full polytope \eqn{\mathcal{P}}).
#'
#'
#' @param lim A list with four components \code{A}, \code{B}, \code{G} and \code{H} representing
#' the polytope to be reduced.
#' @param test A boolean. If \code{TRUE} (the default), implicit equalities hidden in the
#' inequalities are detected and added to the equalities; if \code{FALSE}, they raise an error.
#'
#' @return A list with five components; namely:
#' \itemize{
#'   \item \code{G}
#'   \item \code{H}
#'   \item \code{x0}
#'   \item \code{Z}
#'   \item \code{dropped}, a logical vector flagging the inequalities of \code{lim}
#'     removed because constant on the affine hull of \eqn{\mathcal{P}}.
#' }
#' 
#' @rdname lim.redpol
#' @importFrom MASS Null
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
#'
#' @examples
#' DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
#' BOWF <- df2lim(DF)
#' BOWFred <- lim.redpol(BOWF)
#' str(BOWFred, max.length = 1)
lim.redpol <- function(lim, test = TRUE) {
  A <- lim$A
  B <- lim$B
  G <- lim$G
  H <- lim$H
  tol <- sqrt(.Machine$double.eps) #the smallest positive floating-point number x such that 1 + x != 1 on the current machine
  ## 0. Setup problem

  if (is.null(A) || length(A) == 0L) {
    stop("no equalities found")
  }
  if (is.null(G) || length(G) == 0L) {
    stop("no inequalities found")
  }
  if (is.vector(A)) A <- t(A)
  if (is.vector(G)) G <- t(G)
  A <- as.matrix(A); G <- as.matrix(G)
  B <- as.numeric(B); H <- as.numeric(H)

  ## 1. Reference point x0: a Chebyshev center; implicit equalities detected on the way
  ctr <- .chebyshev_hull(A, B, G, H, tol, test)
  x0 <- ctr$center
  Z <- ctr$Z

  ## 2. Projection of G and H onto reduced space
  g <- G %*% Z
  h <- H - G %*% x0
  g[abs(g) < tol * sqrt(rowSums(G^2))] <- 0          # relative to the norm of row i
  # h is not thresholded: zeroing a small margin would put the origin on the boundary.
  h <- as.numeric(h)

  ## 3. Drop the inequalities constant on the affine hull (zero rows of g)
  dropped <- rowSums(g != 0) == 0

  return(list("G" = g[!dropped, , drop = FALSE], "H" = h[!dropped], "x0" = x0, "Z" = Z,
              "dropped" = dropped))
}


# Chebyshev center of {x : Ax = B, Gx >= H} within its affine hull: a radius zero in
# {Ax = B} reveals implicit equalities, added to A before a second linear program.
.chebyshev_hull <- function(A, B, G, H, tol, test = TRUE) {
  # Only a polytope smaller than 1 is scaled up: scaling a larger one down would
  # inflate GLPK's feasibility tolerance (1e-7 relative to 1 + |bound|), and the
  # center of a thin polytope would no longer be strictly interior.
  sc <- max(abs(c(B, H)))
  if (!is.finite(sc) || sc == 0) sc <- 1
  sc <- min(1, sc)
  B <- B / sc
  H <- H / sc
  ctr <- .chebyshev_affine(A, B, G, H, tol)
  if (ctr$radius <= tol) {
    if (!test) {
      stop("Some inequalities hold as equalities on the whole polytope. ",
           "Use test = TRUE to detect them.")
    }
    implicit <- .implicit_equalities(A, B, G, H, ctr$s, tol)
    A <- rbind(A, G[implicit, , drop = FALSE])
    B <- c(B, H[implicit])
    ctr <- .chebyshev_affine(A, B, G, H, tol)
    if (ctr$radius <= tol) {
      stop("The polytope has an empty interior, even after adding the equalities ",
           "hidden in the inequalities.")
    }
  }
  ctr$center <- ctr$center * sc
  ctr$radius <- ctr$radius * sc
  ctr
}


# Chebyshev center of {x : Ax = B, Gx >= H} within {Ax = B}, by one linear program:
#   maximise r subject to  Ax = B  and  g_i . x - s_i r >= H_i,
# with s_i = ||Z^T g_i||, Z an orthonormal basis of ker(A): the norm of g_i within
# the affine space (s_i = 0 for an inequality constant on it).
.chebyshev_affine <- function(A, B, G, H, tol) {
  Z <- Null(t(A)); Z[abs(Z) < tol] <- 0 #x=x0+Zq ; AZ=0
  if (ncol(Z) == 0L) stop("The equalities fix all the unknowns: the polytope is a point.")
  s <- sqrt(rowSums((G %*% Z)^2))
  s[s <= tol * sqrt(rowSums(G^2))] <- 0
  n <- ncol(A)
  sol <- Rglpk_solve_LP(
    obj = c(numeric(n), 1),
    mat = rbind(cbind(A, 0), cbind(G, -s)),
    dir = c(rep("==", nrow(A)), rep(">=", nrow(G))),
    rhs = c(B, H),
    bounds = list(lower = list(ind = seq_len(n + 1L), val = c(rep(-Inf, n), 0))),
    max = TRUE)
  if (sol$status != 0L) stop("No point found: is the polytope empty or unbounded?")
  list(center = sol$solution[seq_len(n)], radius = sol$solution[n + 1L], Z = Z, s = s)
}


# Implicit equalities of {x : Ax = B, Gx >= H}: the inequalities whose largest margin
# over the polytope is zero, by one linear program each (inequalities constant on
# {Ax = B} excepted); the margin is divided by s_i, a distance in the affine space.
.implicit_equalities <- function(A, B, G, H, s, tol) {
  n <- ncol(A)
  mat <- rbind(A, G)
  dir <- c(rep("==", nrow(A)), rep(">=", nrow(G)))
  rhs <- c(B, H)
  bnd <- list(lower = list(ind = seq_len(n), val = rep(-Inf, n)))
  implicit <- logical(nrow(G))
  for (i in which(s > 0)) {
    sol <- Rglpk_solve_LP(obj = G[i, ], mat = mat, dir = dir, rhs = rhs,
                          bounds = bnd, max = TRUE)
    # A non-zero status means that g_i . x is unbounded above: not an equality.
    if (sol$status == 0L && (sol$optimum - H[i]) / s[i] <= tol) implicit[i] <- TRUE
  }
  implicit
}


#' @rdname lim.redpol
#' @param sample  A matrix where each row corresponds to a point inside either the full or the reduced polytope.
#' @param x0 A numeric vector of size \eqn{n}, the reference point used during the reduction (returned by \code{lim.redpol()}).
#' @param Z The matrix used during the reduction of the polytope (returned by \code{lim.redpol()}).
#' @export

red2full<- function(sample,x0,Z){
  res<-x0+Z%*%t(sample)
  x<-t(res)
  
  return(x)
}

#' @rdname lim.redpol
#' @export

full2red<- function(sample,x0,Z){
  res<-solve(t(Z)%*%Z)%*%t(Z)%*%(t(sample)-x0)
  return(t(res))

}
