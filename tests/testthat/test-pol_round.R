# Tests on the behaviour of the pol.round function ----

# A 40:1 rectangle rotated by 0.6 rad, as Gx >= H: the case rounding must repair.
.sliver <- function(theta = 0.6) {
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2)
  G <- rbind(c(-1, 0), c(0, -1), c(1, 0), c(0, 1)) %*% t(R)
  list(G = G, H = c(-40, -1, -40, -1))
}

.inside <- function(G, H, X) all(G %*% t(as.matrix(X)) >= H - 1e-9)

test_that("pol.round brings the Dikin ellipsoid back to a ball", {
  sl <- .sliver()
  rnd <- pol.round(sl$G, sl$H)
  expect_gt(rnd$axis.ratio, 30)                        # the body IS anisotropic
  # rounding by a method makes its own ellipsoid spherical, by construction
  rnd2 <- pol.round(rnd$G, rnd$H)
  expect_equal(rnd2$axis.ratio, 1, tolerance = 1e-6)
})

test_that("pol.round returns forth and back inverse of one another", {
  sl <- .sliver()
  rnd <- pol.round(sl$G, sl$H)
  set.seed(1)
  X <- matrix(runif(200, -0.4, 0.4), ncol = 2)
  expect_lt(max(abs(rnd$forth(rnd$back(X)) - X)), 1e-9)
  # a lone vector is accepted like a one-row matrix
  expect_lt(max(abs(rnd$back(X[1, ]) - as.numeric(rnd$back(X[1, , drop = FALSE])))), 1e-12)
})

test_that("pol.round returns points that land inside the original polytope", {
  sl <- .sliver()
  rnd <- pol.round(sl$G, sl$H)
  set.seed(2)
  X <- matrix(runif(2000, -0.5, 0.5), ncol = 2)
  X <- X[apply(rnd$G %*% t(X) >= rnd$H, 2, all), , drop = FALSE]
  expect_gt(nrow(X), 50)
  expect_true(.inside(sl$G, sl$H, rnd$back(X)))
})

test_that("pol.round preserves the uniform distribution", {
  sl <- .sliver()
  rnd <- pol.round(sl$G, sl$H)
  rej <- function(G, H, lo, hi, n) {
    out <- matrix(0, 0, 2)
    while (nrow(out) < n) {
      Z <- cbind(runif(6 * n, lo[1], hi[1]), runif(6 * n, lo[2], hi[2]))
      out <- rbind(out, Z[apply(G %*% t(Z) >= H, 2, all), , drop = FALSE])
    }
    out[seq_len(n), , drop = FALSE]
  }
  set.seed(3)
  direct <- rej(sl$G, sl$H, c(-40, -40), c(40, 40), 8000)
  via    <- rnd$back(rej(rnd$G, rnd$H, c(-2, -2), c(2, 2), 8000))
  # the Jacobian being constant, both samples follow the same law
  gap <- abs(colMeans(direct) - colMeans(via)) / apply(direct, 2, stats::sd)
  expect_lt(max(gap), 0.08)
})

test_that("pol.round rounds better at the analytic center than at the Chebyshev center", {
  # an asymmetric triangle, where the two centers differ markedly
  G <- rbind(c(1, 0), c(0, 1), c(-1 / 40, -1 / 2)); H <- c(0, 0, -1)
  ra <- pol.round(G, H, center = "analytic")
  rt <- pol.round(G, H, center = "chebyshev")
  eta <- pol.ranges(G = ra$G, H = ra$H)[, 3]
  etc <- pol.ranges(G = rt$G, H = rt$H)[, 3]
  expect_lt(max(eta) / min(eta), max(etc) / min(etc))
  expect_gt(ra$log.volume, rt$log.volume)
})

test_that("pol.round collapses the anisotropy of BOWF-short", {
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  exf <- lim.exfoliate(lim.redpol(df2lim(DF)))
  rnd <- pol.round(G = exf$G, H = exf$H)
  rg0 <- pol.ranges(G = exf$G, H = exf$H)[, 3]
  rg1 <- pol.ranges(G = rnd$G, H = rnd$H)[, 3]
  expect_gt(max(rg0) / min(rg0), 100)          # before: heavily anisotropic
  expect_lt(max(rg1) / min(rg1), 6)            # after: comfortable
})

test_that("pol.round rejects malformed input", {
  sl <- .sliver()
  expect_error(pol.round(sl$G, sl$H[-1]), "incompatible")
  expect_error(pol.round(NULL, sl$H), "0 dimensions")
})
