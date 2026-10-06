# Tests on the behaviour of the pol.round function ----

# A 40:1 rectangle rotated by 0.6 rad, as Gx >= H.
.sliver <- function(theta = 0.6) {
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2)
  list(G = rbind(c(-1, 0), c(0, -1), c(1, 0), c(0, 1)) %*% t(R), H = c(-40, -1, -40, -1))
}

.anisotropy <- function(G, H) { r <- pol.ranges(G = G, H = H)[, 3]; max(r) / min(r) }

test_that("pol.round maps the Dikin ellipsoid onto the unit ball", {
  sl <- .sliver()
  # before: the Dikin ellipsoid at the analytic center has semi-axes in ratio 40
  d <- as.numeric(sl$G %*% pol.center(G = sl$G, H = sl$H)) - sl$H
  ev <- eigen(crossprod(sl$G, sl$G / d^2))$values
  expect_equal(sqrt(max(ev) / min(ev)), 40, tolerance = 1e-6)
  # after: the Hessian of the log-barrier at the origin is the identity
  rnd <- pol.round(sl$G, sl$H)
  expect_equal(crossprod(rnd$G, rnd$G / rnd$H^2), diag(2), tolerance = 1e-8)
})

test_that("pol.round maps the rounded polytope onto the original one", {
  sl <- .sliver()
  rnd <- pol.round(sl$G, sl$H)
  set.seed(1)
  Y <- matrix(runif(2000, -1, 1), ncol = 2)
  inside <- apply(rnd$G %*% t(Y) >= rnd$H, 2, all)
  X <- red2full(Y, rnd$x0, rnd$Z)
  expect_identical(apply(sl$G %*% t(X) >= sl$H - 1e-9, 2, all), inside)
  expect_lt(max(abs(full2red(X, rnd$x0, rnd$Z) - Y)), 1e-9)
})

test_that("pol.round collapses the anisotropy of BOWF-short", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  exf <- lim.exfoliate(lim.redpol(df2lim(DF)))
  rnd <- pol.round(exf$G, exf$H)
  expect_gt(.anisotropy(exf$G, exf$H), 100)
  expect_lt(.anisotropy(rnd$G, rnd$H), 6)
})

test_that("pol.round ignores a zero row, which always holds", {
  sl <- .sliver()
  rnd <- pol.round(rbind(sl$G, c(0, 0)), c(sl$H, -1))
  expect_equal(rnd$Z, pol.round(sl$G, sl$H)$Z, tolerance = 1e-8)
})

test_that("pol.round rejects malformed input", {
  sl <- .sliver()
  expect_error(pol.round(sl$G, sl$H[-1]), "incompatible")
  expect_error(pol.round(NULL, sl$H), "0 dimensions")
})

# Tests on the behaviour of the lim.round function ----

test_that("lim.round composes the rounding with the reduction", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  BOWF <- df2lim(DF)
  rnd <- lim.round(lim.exfoliate(lim.redpol(BOWF)))
  # the origin of the rounded polytope is the analytic center, mapped onto the unknowns
  x <- as.numeric(red2full(t(numeric(ncol(rnd$Z))), rnd$x0, rnd$Z))
  expect_lt(max(abs(BOWF$A %*% x - BOWF$B)), 1e-8)
  expect_gt(min(BOWF$G %*% x - BOWF$H), 0)
  expect_error(lim.round(BOWF), "lim.redpol")
})
