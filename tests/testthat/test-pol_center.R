# Tests on the behaviour of the pol.center function ----

.square <- function() list(G = rbind(c(-1, 0), c(0, -1), c(1, 0), c(0, 1)),
                           H = c(-1, -1, -1, -1))

# Gradient of the log-barrier: it vanishes at the analytic center, by definition.
.grad <- function(G, H, x) {
  d <- as.numeric(G %*% x) - H
  sqrt(sum(as.numeric(crossprod(G, 1 / d))^2))
}

test_that("pol.center returns the origin as analytic center of the square", {
  sq <- .square()
  ctr <- pol.center(sq$G, sq$H, type = "analytic", x0 = c(0.4, -0.3))
  # Test the gradient rather than the position: the stopping rule bounds the gap on
  # the objective, not on the position (see the documentation).
  expect_lt(.grad(sq$G, sq$H, ctr), 1e-6)
  expect_lt(max(abs(ctr)), 1e-6)
})

test_that("pol.center returns an analytic center invariant under row scaling", {
  sq <- .square()
  c1 <- pol.center(sq$G, sq$H, x0 = c(0.4, -0.3))
  s <- c(1, 7, 0.2, 30)
  c2 <- pol.center(sq$G * s, sq$H * s, x0 = c(0.4, -0.3))
  expect_lt(max(abs(c1 - c2)), 1e-7)
})

test_that("pol.center returns an analytic center independent of the starting point", {
  set.seed(3); d <- 5
  G <- matrix(rnorm(40 * d), ncol = d); G <- G / sqrt(rowSums(G^2))
  H <- rep(-1, 40)                                # G x >= -1 contains the origin
  c1 <- pol.center(G, H, x0 = numeric(d))
  c2 <- pol.center(G, H, x0 = c1 * 0.1)
  expect_lt(.grad(G, H, c1), 1e-6)
  expect_lt(max(abs(c1 - c2)), 1e-6)
})

test_that("pol.center returns an analytic center that a redundant constraint moves", {
  sq <- .square()
  c1 <- pol.center(sq$G, sq$H, x0 = c(0, 0))
  # the same half-space written twice: the polytope is unchanged, the center is not
  G2 <- rbind(sq$G, sq$G[1, , drop = FALSE]); H2 <- c(sq$H, sq$H[1])
  c2 <- pol.center(G2, H2, x0 = c(0, 0))
  expect_gt(max(abs(c1 - c2)), 1e-3)
})

test_that("pol.center returns the origin and radius 1 as Chebyshev center of the square", {
  sq <- .square()
  ch <- pol.center(sq$G, sq$H, type = "chebyshev")
  expect_lt(max(abs(ch$center)), 1e-8)
  expect_equal(ch$radius, 1, tolerance = 1e-8)
})

test_that("pol.center rejects a non-interior starting point", {
  sq <- .square()
  expect_error(pol.center(sq$G, sq$H, x0 = c(1, 0)), "strictly interior")
})

test_that("pol.center returns two distant centers on BOWF-short", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  exf <- lim.exfoliate(lim.redpol(df2lim(DF)))
  ch <- pol.center(exf$G, exf$H, type = "chebyshev")
  an <- pol.center(exf$G, exf$H, type = "analytic")
  expect_lt(.grad(exf$G, exf$H, an), 1e-5)
  # the two centers are separated by far more than the inscribed radius
  expect_gt(sqrt(sum((ch$center - an)^2)), 100 * ch$radius)
})
