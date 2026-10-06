# Tests on the behaviour of the pol.center function ----

# The square [-1, 1]^2, as Gx >= H.
.square <- function() list(G = rbind(c(-1, 0), c(0, -1), c(1, 0), c(0, 1)),
                           H = c(-1, -1, -1, -1))

# The same square in the plane x3 = 1/2 of R^3, with 0 <= x3 <= 1: these two bounds are
# constant on {Ax = B}, and must be ignored.
.square3 <- function() list(A = matrix(c(0, 0, 1), 1), B = 0.5,
                            G = rbind(cbind(.square()$G, 0), c(0, 0, 1), c(0, 0, -1)),
                            H = c(.square()$H, 0, -1))

# Norm of the gradient of the log-barrier, which vanishes at the analytic center.
.grad <- function(G, H, x) sqrt(sum(crossprod(G, 1 / (as.numeric(G %*% x) - H))^2))

test_that("pol.center finds the centers of the square", {
  sq <- .square()
  ch <- pol.center(G = sq$G, H = sq$H, type = "chebyshev")
  expect_equal(ch$center, c(0, 0), tolerance = 1e-8)
  expect_equal(ch$radius, 1, tolerance = 1e-8)
  expect_equal(pol.center(G = sq$G, H = sq$H, x0 = c(0.4, -0.3)), c(0, 0), tolerance = 1e-6)
})

test_that("pol.center returns the same analytic center from any starting point", {
  set.seed(3)
  G <- matrix(rnorm(40 * 5), ncol = 5); G <- G / sqrt(rowSums(G^2)); H <- rep(-1, 40)
  c1 <- pol.center(G = G, H = H)
  c2 <- pol.center(G = G, H = H, x0 = c1 / 10)
  expect_lt(.grad(G, H, c1), 1e-6)
  expect_lt(max(abs(c1 - c2)), 1e-6)
  # and the same after scaling the rows
  expect_lt(max(abs(c1 - pol.center(G = G * (1:40), H = H * (1:40)))), 1e-6)
})

test_that("pol.center returns an analytic center that a redundant constraint moves", {
  sq <- .square()
  G2 <- rbind(sq$G, sq$G[1, ]); H2 <- c(sq$H, sq$H[1])   # same polytope
  expect_gt(max(abs(pol.center(G = G2, H = H2))), 1e-3)
})

test_that("pol.center ignores a zero row, which always holds", {
  sq <- .square()
  G <- rbind(sq$G, c(0, 0)); H <- c(sq$H, -1)
  expect_equal(pol.center(G = G, H = H, type = "chebyshev")$radius, 1, tolerance = 1e-8)
  expect_equal(pol.center(G = G, H = H), c(0, 0), tolerance = 1e-6)
})

test_that("pol.center finds the centers within the affine space {Ax = B}", {
  sq <- .square3()
  ch <- pol.center(A = sq$A, B = sq$B, G = sq$G, H = sq$H, type = "chebyshev")
  expect_equal(ch$center, c(0, 0, 0.5), tolerance = 1e-8)
  expect_equal(ch$radius, 1, tolerance = 1e-8)
  expect_equal(pol.center(A = sq$A, B = sq$B, G = sq$G, H = sq$H), c(0, 0, 0.5),
               tolerance = 1e-6)
})

test_that("pol.center rejects a bad starting point", {
  sq <- .square()
  expect_error(pol.center(G = sq$G, H = sq$H, x0 = c(1, 0)), "strictly interior")
  sq <- .square3()
  expect_error(pol.center(A = sq$A, B = sq$B, G = sq$G, H = sq$H, x0 = c(0, 0, 0.6)),
               "Ax = B")
})

test_that("pol.center refuses a polytope without interior point", {
  # x1 + x2 + x3 + x4 = 2, x1 + x2 >= 1, x3 + x4 >= 1, 0 <= x <= 1 force x1 + x2 = 1
  G <- rbind(c(1, 1, 0, 0), c(0, 0, 1, 1), diag(4), -diag(4))
  H <- c(1, 1, rep(0, 4), rep(-1, 4))
  expect_error(pol.center(A = matrix(1, 1, 4), B = 2, G = G, H = H, type = "chebyshev"),
               "No interior point")
})

# Tests on the behaviour of the lim.center function ----

test_that("lim.center gives on the full polytope the centers of the reduced one", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  BOWF <- df2lim(DF)
  red <- lim.redpol(BOWF)
  # Chebyshev: same radius, at a point of {Ax = B}
  ch <- lim.center(BOWF, type = "chebyshev")
  expect_equal(ch$radius, lim.center(red, type = "chebyshev")$radius, tolerance = 1e-8)
  expect_lt(max(abs(BOWF$A %*% ch$center - BOWF$B)), 1e-8)
  # Analytic: the same point, the analytic center being affinely invariant
  an <- lim.center(BOWF)
  expect_lt(max(abs(BOWF$A %*% an - BOWF$B)), 1e-8)
  expect_equal(an, as.numeric(red2full(t(lim.center(red)), red$x0, red$Z)),
               tolerance = 1e-6)
  expect_error(lim.center(list(x0 = 1)), "components G and H")
})
