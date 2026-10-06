# Checks shared by the pre-processing functions.

# Description {x : Gx >= H} of a polytope: G as a matrix, H as a vector, with matching
# dimensions and finite entries. A zero row of G, as left by lim.redpol() for an
# inequality constant on {Ax = B}, always holds when H <= 0; with H > 0, the polytope
# is empty.
.pol_check <- function(G, H) {
  if (is.null(G)) stop("G is NULL, the polytope has 0 dimensions.")
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.vector(G)) G <- t(G)
  H <- as.numeric(H)
  if (nrow(G) != length(H)) stop("G and H have incompatible dimensions.")
  if (any(!is.finite(G)) || any(!is.finite(H))) stop("G and H must be finite.")
  if (any(rowSums(G != 0) == 0 & H > 0)) stop("The polytope is empty: 0 >= H with H > 0.")
  list(G = G, H = H)
}

# Exfoliation and rounding apply to a reduced polytope {x : Gx >= H}: a list still
# carrying equality constraints is refused, rather than its equalities ignored.
.lim_check_reduced <- function(lim) {
  if (!is.list(lim) || is.null(lim$G) || is.null(lim$H))
    stop("`lim` must be a list with components G and H.", call. = FALSE)
  if (length(lim$A) > 0L)
    stop("Equality constraints found in `lim$A`: apply lim.redpol() first.", call. = FALSE)
  invisible(TRUE)
}
