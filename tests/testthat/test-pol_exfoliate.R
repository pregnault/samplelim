# Tests on the behaviour of the pol.exfoliate function ----

# The square [-1, 1]^2, as Gx >= H.
.square <- function() list(G = rbind(c(-1, 0), c(0, -1), c(1, 0), c(0, 1)),
                           H = c(-1, -1, -1, -1))

# Support function of {Gx >= H} in direction u: two descriptions define the same
# polytope if and only if they have the same support function in every direction.
.support <- function(G, H, u) {
  d <- ncol(G)
  sol <- Rglpk::Rglpk_solve_LP(u, G, rep(">=", nrow(G)), H, max = TRUE,
                               bounds = list(lower = list(ind = seq_len(d),
                                                          val = rep(-Inf, d))))
  if (sol$status != 0L) NA_real_ else sol$optimum
}

.same_polytope <- function(G1, H1, G2, H2, n = 50, seed = 1) {
  set.seed(seed)
  U <- matrix(rnorm(n * ncol(G1)), ncol = ncol(G1))
  max(vapply(seq_len(n), function(i)
    abs(.support(G1, H1, U[i, ]) - .support(G2, H2, U[i, ])), numeric(1)))
}

test_that("pol.exfoliate removes an implied and a duplicated constraint", {
  sq <- .square()
  # x >= -5 is implied by x >= -1 (row 3), which is also written twice
  G <- rbind(sq$G, c(1, 0), c(1, 0)); H <- c(sq$H, -5, -1)
  exf <- pol.exfoliate(G = G, H = H)
  expect_true(exf$redundant[5])
  expect_equal(sum(exf$redundant[c(3, 6)]), 1)   # one of the two copies only
  expect_equal(nrow(exf$G), 4)
  expect_lt(.same_polytope(G, H, exf$G, exf$H), 1e-9)
})

test_that("pol.exfoliate leaves a description without redundancy untouched", {
  sq <- .square()
  exf <- pol.exfoliate(G = sq$G, H = sq$H)
  expect_false(any(exf$redundant))
  expect_false(any(exf$implicit))
  expect_null(exf$A)
  expect_equal(exf$G, sq$G)
})

test_that("pol.exfoliate finds the 28 redundant constraints of BOWF-short", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  BOWF <- df2lim(DF)
  red <- lim.redpol(BOWF)
  exf <- pol.exfoliate(G = red$G, H = red$H)
  expect_equal(nrow(red$G), 72)
  expect_equal(sum(exf$redundant), 28)
  expect_lt(.same_polytope(red$G, red$H, exf$G, exf$H, n = 30), 1e-6)
  # on the full polytope, with Ax = B in every program: the same 28 constraints
  full <- pol.exfoliate(A = BOWF$A, B = BOWF$B, G = BOWF$G, H = BOWF$H)
  expect_identical(full$redundant, exf$redundant)
  expect_false(any(full$implicit))
  expect_identical(full$A, BOWF$A)
  # without the equalities, the programs miss some of them
  expect_lt(sum(pol.exfoliate(G = BOWF$G, H = BOWF$H)$redundant), 28)
})

test_that("pol.exfoliate moves the implicit equalities to A", {
  # x1 + x2 + x3 + x4 = 2, x1 + x2 >= 1, x3 + x4 >= 1, 0 <= x <= 1 force x1 + x2 = 1
  # and x3 + x4 = 1, a fixed combination of unknowns, none of them fixed
  A <- matrix(1, 1, 4); B <- 2
  G <- rbind(c(1, 1, 0, 0), c(0, 0, 1, 1), diag(4), -diag(4))
  H <- c(1, 1, rep(0, 4), rep(-1, 4))
  exf <- pol.exfoliate(A = A, B = B, G = G, H = H)
  expect_equal(exf$implicit, c(TRUE, TRUE, rep(FALSE, 8)))
  # x3 + x4 = 1 follows from the two others: only x1 + x2 = 1 is appended
  expect_equal(exf$A, rbind(A, G[1, ]))
  expect_equal(exf$B, c(B, 1))
  # then x1 <= 1 follows from x1 + x2 = 1 and x2 >= 0, and so on: 4 of the 8 bounds
  expect_equal(sum(exf$redundant), 4)
  expect_equal(nrow(exf$G), 4)
})

test_that("pol.exfoliate moves an equality written as two inequalities", {
  # the segment {x2 = 0, -1 <= x1 <= 1}, written x2 >= 0 and -x2 >= 0
  G <- rbind(c(-1, 0), c(1, 0), c(0, 1), c(0, -1)); H <- c(-1, -1, 0, 0)
  exf <- pol.exfoliate(G = G, H = H)
  expect_equal(exf$implicit, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(exf$A, matrix(c(0, 1), 1))     # one of the two, A keeping full rank
  expect_equal(exf$B, 0)
  expect_equal(nrow(exf$G), 2)
})

test_that("pol.exfoliate removes a zero row, which always holds", {
  sq <- .square()
  exf <- pol.exfoliate(G = rbind(sq$G, c(0, 0)), H = c(sq$H, -1))
  expect_equal(exf$redundant, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_error(pol.exfoliate(G = rbind(sq$G, c(0, 0)), H = c(sq$H, 1)), "empty")
})

test_that("pol.exfoliate rejects malformed input", {
  sq <- .square()
  expect_error(pol.exfoliate(G = sq$G, H = sq$H[-1]), "incompatible")
  expect_error(pol.exfoliate(G = NULL, H = sq$H), "0 dimensions")
  expect_error(pol.exfoliate(A = matrix(1, 1, 3), B = 1, G = sq$G, H = sq$H), "incompatible")
})

# Tests on the behaviour of the lim.exfoliate function ----

test_that("lim.exfoliate keeps the other components of a reduced polytope", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  red <- lim.redpol(df2lim(DF))
  out <- lim.exfoliate(red)
  expect_equal(nrow(out$G), 44)
  expect_named(out, names(red))
  expect_identical(out$x0, red$x0)
  expect_identical(out$Z, red$Z)
})

test_that("lim.exfoliate works on the full polytope, before the reduction", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  BOWF <- df2lim(DF)
  exf <- lim.exfoliate(BOWF)
  expect_equal(nrow(exf$G), 44)
  expect_identical(exf$A, BOWF$A)
  expect_named(exf, names(BOWF))
  # exfoliate, then reduce: the 44 constraints of reduce, then exfoliate, no zero row
  red <- lim.redpol(exf)
  expect_equal(nrow(red$G), 44)
  expect_true(all(rowSums(red$G != 0) > 0))
  expect_equal(red$G, lim.exfoliate(lim.redpol(BOWF))$G, tolerance = 1e-8)
  expect_error(lim.exfoliate(list(x0 = 1)), "components G and H")
})

test_that("lim.exfoliate refuses a reduced polytope that is flat", {
  flat <- list(G = rbind(c(-1, 0), c(1, 0), c(0, 1), c(0, -1)), H = c(-1, -1, 0, 0),
               x0 = c(0, 0), Z = diag(2))
  expect_error(lim.exfoliate(flat), "exfoliate the full polytope")
})

# Ill-conditioned descriptions ----

test_that("pol.exfoliate measures within {Ax = B}: a row almost along A is kept", {
  # simplex x1 + x2 + x3 = 1, x >= 0, and g x >= M - 0.9 with g = M (1, 1, 1) + (1, -1, 0),
  # that is x1 - x2 >= -0.9 on the simplex: not redundant, whatever M
  for (M in c(0, 1e3, 1e8)) {
    G <- rbind(diag(3), M * c(1, 1, 1) + c(1, -1, 0)); H <- c(0, 0, 0, M - 0.9)
    exf <- pol.exfoliate(A = matrix(1, 1, 3), B = 1, G = G, H = H)
    expect_false(exf$redundant[4])
  }
})

test_that("pol.exfoliate appends an implicit equality whatever the scale of the rows", {
  # 1e4 (x1 + x2) = 2e4, 1e-4 (x1 - x2) >= 0 and 1e-4 (x2 - x1) >= 0: x1 = x2 = 1
  G <- rbind(1e-4 * c(1, -1, 0), 1e-4 * c(-1, 1, 0), rbind(diag(3), -diag(3)))
  H <- c(0, 0, 0, 0, 0, -3, -3, -1)
  exf <- pol.exfoliate(A = 1e4 * matrix(c(1, 1, 0), 1), B = 2e4, G = G, H = H)
  expect_equal(exf$implicit[1:2], c(TRUE, TRUE))
  expect_equal(qr(exf$A / sqrt(rowSums(exf$A^2)))$rank, 2)
  expect_gt(max(abs(exf$A %*% c(2, 0, 0.5) - exf$B)), 1e-6)   # (2, 0, .) is not in P
})

test_that("pol.exfoliate keeps a row almost constant on {Ax = B} that cuts the polytope", {
  # x1 + x3 = 0 and (1 + 1e-8) x1 + x3 >= 0, that is 1e-8 x1 >= 0
  G <- rbind(c(1 + 1e-8, 0, 1), rbind(diag(3), -diag(3))); H <- c(0, rep(-1, 6))
  exf <- pol.exfoliate(A = matrix(c(1, 0, 1), 1), B = 0, G = G, H = H)
  expect_false(exf$redundant[1])
})

test_that("pol.exfoliate takes as value of an implicit equality the middle of its range", {
  # 1e-10 <= x1 <= 2e-10, 0 <= x2 <= 1, and the redundant x1 >= 0 first
  G <- rbind(c(1, 0), c(1, 0), c(-1, 0), c(0, 1), c(0, -1)); H <- c(0, 1e-10, -2e-10, 0, -1)
  exf <- pol.exfoliate(G = G, H = H)
  expect_equal(nrow(exf$A), 1)
  v <- exf$B / exf$A[1, 1]
  expect_gte(v, 1e-10); expect_lte(v, 2e-10)
})

test_that("pol.exfoliate moves one implicit equality at a time", {
  # x2 = 0, written x2 >= 0 and -x2 >= 0, 0 <= x1 <= 1 and x1 + M x2 >= 0.05: the segment
  # [0.05, 1] x {0}. Within R^2, the last row is almost along x2 = 0, and nearly tight on
  # the polytope; once x2 = 0 is in A, it is the bound x1 >= 0.05, far from tight.
  G <- rbind(c(0, 1), c(0, -1), c(1, 0), c(-1, 0), c(1, 1e8)); H <- c(0, 0, 0, -1, 0.05)
  exf <- pol.exfoliate(G = G, H = H)
  expect_equal(exf$implicit, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(exf$redundant, c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(exf$A, matrix(c(0, 1), 1))
  expect_equal(exf$B, 0)
  expect_equal(exf$G, G[4:5, ])
})

test_that("pol.exfoliate refuses a polytope too small for the tolerances", {
  # equilateral triangle of inradius 2.5e-8: no ball of radius 3e-8, and no constraint
  # within 6e-8 of an implicit equality
  G <- rbind(c(0, 1), c(-sqrt(3) / 2, -1 / 2), c(sqrt(3) / 2, -1 / 2)); H <- rep(-2.5e-8, 3)
  expect_error(pol.exfoliate(G = G, H = H), "Chebyshev radius")
  expect_false(any(pol.exfoliate(G = G, H = 1e6 * H)$implicit))     # rescaled
})

test_that("pol.exfoliate handles unbounded polyhedra and refuses empty ones", {
  exf <- pol.exfoliate(G = diag(2), H = c(0, 0))                    # the orthant
  expect_false(any(exf$redundant | exf$implicit))
  expect_error(pol.exfoliate(A = matrix(c(1, 1), 1), B = 1, G = rbind(c(1, 1)), H = 2),
               "empty")
})

