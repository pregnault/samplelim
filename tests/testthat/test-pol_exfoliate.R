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
  exf <- pol.exfoliate(G, H)
  expect_true(exf$redundant[5])
  expect_equal(sum(exf$redundant[c(3, 6)]), 1)   # one of the two copies only
  expect_equal(nrow(exf$G), 4)
  expect_lt(.same_polytope(G, H, exf$G, exf$H), 1e-9)
})

test_that("pol.exfoliate leaves a description without redundancy untouched", {
  sq <- .square()
  exf <- pol.exfoliate(sq$G, sq$H)
  expect_false(any(exf$redundant))
  expect_equal(exf$G, sq$G)
})

test_that("pol.exfoliate finds the 28 redundant constraints of BOWF-short", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  red <- lim.redpol(df2lim(DF))
  exf <- pol.exfoliate(G = red$G, H = red$H)
  expect_equal(nrow(red$G), 72)
  expect_equal(sum(exf$redundant), 28)
  expect_lt(.same_polytope(red$G, red$H, exf$G, exf$H, n = 30), 1e-6)
})

test_that("pol.exfoliate removes a zero row, which always holds", {
  sq <- .square()
  exf <- pol.exfoliate(rbind(sq$G, c(0, 0)), c(sq$H, -1))
  expect_equal(exf$redundant, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_error(pol.exfoliate(rbind(sq$G, c(0, 0)), c(sq$H, 1)), "empty")
})

test_that("pol.exfoliate rejects malformed input", {
  sq <- .square()
  expect_error(pol.exfoliate(sq$G, sq$H[-1]), "incompatible")
  expect_error(pol.exfoliate(NULL, sq$H), "0 dimensions")
})

# Tests on the behaviour of the lim.exfoliate function ----

test_that("lim.exfoliate replaces G and H only", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  red <- lim.redpol(df2lim(DF))
  out <- lim.exfoliate(red)
  expect_equal(nrow(out$G), 44)
  expect_named(out, names(red))
  expect_identical(out$x0, red$x0)
  expect_identical(out$Z, red$Z)
})

test_that("lim.exfoliate refuses a polytope that is not reduced", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  expect_error(lim.exfoliate(df2lim(DF)), "lim.redpol")
  expect_error(lim.exfoliate(list(x0 = 1)), "components G and H")
})
