# Tests on the behaviour of the pol.exfoliate function ----

# A square {-1 <= x <= 1, -1 <= y <= 1}, written as Gx >= H.
.square <- function() list(G = rbind(c(-1, 0), c(0, -1), c(1, 0), c(0, 1)),
                           H = c(-1, -1, -1, -1))

# Support function of {Gx >= H} in direction u, by linear programming. Two
# descriptions define the same polytope if and only if they have the same support
# function in every direction. Counting rows proves nothing.
.support <- function(G, H, u) {
  d <- ncol(G)
  sol <- Rglpk::Rglpk_solve_LP(u, G, rep(">=", nrow(G)), H,
                               bounds = list(lower = list(ind = seq_len(d),
                                                          val = rep(-Inf, d))),
                               max = TRUE)
  if (sol$status != 0L) NA_real_ else sol$optimum
}

test_that("pol.exfoliate removes an implied constraint and keeps the polytope", {
  sq <- .square()
  # Row 3 of .square() is x >= -1. We append x >= -5, which it implies, and a second
  # copy of x >= -1, so that rows 3 and 6 are strictly identical.
  G <- rbind(sq$G, c(1, 0), c(1, 0))
  H <- c(sq$H, -5, -1)

  exf <- pol.exfoliate(G, H)

  # the useless constraint is removed
  expect_true(exf$redundant[5])
  # of the two copies of x >= -1 exactly one survives; the sequential sweep drops
  # the first one it meets, hence the test on the count and not on which one
  expect_equal(sum(exf$redundant[c(3, 6)]), 1)
  expect_equal(nrow(exf$G), 4)

  set.seed(1)
  U <- matrix(rnorm(50 * 2), ncol = 2); U <- U / sqrt(rowSums(U^2))
  ecart <- vapply(seq_len(nrow(U)), function(i)
    abs(.support(G, H, U[i, ]) - .support(exf$G, exf$H, U[i, ])), numeric(1))
  expect_lt(max(ecart), 1e-9)
})

test_that("pol.exfoliate leaves a minimal description untouched", {
  sq <- .square()
  exf <- pol.exfoliate(sq$G, sq$H)
  expect_equal(sum(exf$redundant), 0)
  expect_equal(dim(exf$G), dim(sq$G))
})

test_that("pol.exfoliate finds the 28 redundant constraints of BOWF-short", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  red <- lim.redpol(df2lim(DF))
  exf <- pol.exfoliate(G = red$G, H = red$H)

  expect_equal(nrow(red$G), 72)
  expect_equal(sum(exf$redundant), 28)
  expect_equal(nrow(exf$G), 44)

  # and the body is unchanged
  set.seed(2)
  U <- matrix(rnorm(30 * ncol(red$G)), ncol = ncol(red$G))
  U <- U / sqrt(rowSums(U^2))
  ecart <- vapply(seq_len(nrow(U)), function(i)
    abs(.support(red$G, red$H, U[i, ]) - .support(exf$G, exf$H, U[i, ])), numeric(1))
  expect_lt(max(ecart), 1e-6)
})

test_that("pol.exfoliate rejects malformed input", {
  sq <- .square()
  expect_error(pol.exfoliate(sq$G, sq$H[-1]), "incompatible")
  expect_error(pol.exfoliate(NULL, sq$H), "0 dimensions")
  expect_error(pol.exfoliate(rbind(sq$G, c(0, 0)), c(sq$H, 0)), "Degenerate")
})

test_that("lim.exfoliate keeps the lim structure", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  red <- lim.redpol(df2lim(DF))
  out <- lim.exfoliate(red)
  expect_equal(nrow(out$G), 44)
  expect_equal(length(out$H), 44)
  # only G and H change; every other component of the lim object is untouched
  expect_true(identical(out$A, red$A))
  expect_true(identical(out$B, red$B))
  expect_equal(out$NConstraints, 44)
})
