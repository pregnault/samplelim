# Tests on the behaviour of rlim ----

#### Using rlim() for sampling LIM ####

# A function to check if points are within range
is_within_range <- function(point, ranges) {
  all(ranges[, 1] <= point & point <= ranges[, 2])
}


test_that("rlim() works as expected for LIM object", {
  # lim is mandatory
  expect_error(rlim())
  # Create a LIM object from a Declaration File
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  BOWF <- df2lim(DF)
  nsamp <- 20
  burn <- 0
  thin <- 1
  samp <- rlim(lim = BOWF, nsamp = nsamp, burn = burn, thin = thin,
               seed = 456)
  # Check dimensions of sample are correct
  expect_equal(dim(samp), c(nsamp,ncol(BOWF$A)))
  
  # Check type is correct
  expect_type(samp, "double")
  
  # Check if every point lies within the polytope i.e satistifies both equalities and inequalities
  
  ## Check if every point satisfies all the equalities
  eq_lin_combs <- BOWF$A %*% t(samp)
  apply(eq_lin_combs, 2, function(col) {
    expect_equal(col, BOWF$B)
  })
  
  ## Check if every point satisfies all the inequalities
  ineq_lin_combs <- BOWF$G %*% t(samp)
  apply(ineq_lin_combs, 2, function(col) {
    expect_true(all(col >= BOWF$H))
  })
  
  # Check if every points are within the ranges
  ranges <- lim.ranges(BOWF)
  result <- apply(samp, 1, is_within_range, ranges=ranges)
  expect_true(all(result))
  
  # Check if seed work as expected
  # same seed returns the same values 
  samp2 <- rlim(lim = BOWF, nsamp = nsamp, burn = burn, thin = thin,
             seed = 456)
  expect_equal(samp2,samp)
  #different seed returns a different value
  samp3 <- rlim(lim = BOWF, nsamp = nsamp, burn = burn, thin = thin,
              seed = 123)
  expect_true(any(samp3 != samp))

  # Check if all type of walk work
  # Billard Walk
  samp_BiW <- rlim(lim = BOWF, type = "BiW", nsamp = nsamp, seed = 123)
  expect_equal(dim(samp_BiW), c(nsamp,ncol(BOWF$A)))
  expect_type(samp_BiW, "double")
  
  # Mirror Walk
  samp_MiW <- rlim(lim = BOWF, type = "MiW", nsamp = nsamp, seed = 123)
  expect_equal(dim(samp_MiW), c(nsamp,ncol(BOWF$A)))
  expect_type(samp_MiW, "double")
  
  # Nothing specified (expect to be equal to Mirror Walk)
  samp_None <- rlim(lim = BOWF, nsamp = nsamp, seed = 123)
  expect_equal(samp_None, samp_MiW)

  # Non-existing walk
  expect_error(rlim(lim = BOWF, type = "This_walk_type_does_not_exist", nsamp = nsamp, seed = 123), "walk type")
  
  
  #Test that thinning works
  no_thin <- rlim(lim = BOWF, nsamp = nsamp, thin=1, seed = 123)
  with_thin <-  rlim(lim = BOWF, nsamp = nsamp, thin = 2, seed = 123)
    
  keeped_rows <- seq(2,dim(no_thin)[1]/2,2)
  expect_equal(with_thin[1:length(keeped_rows),],no_thin[keeped_rows,])
  
})


test_that("rlim() works as expected for Hpol object", {

  A <- matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
  b <- c(0,0,1)
  P <- Hpolytope(A = A, b = b)
  nsamp <- 20
  burn <- 0
  thin <- 1
  
  #if there is Hpolytope object, lim must be NULL and an error is expected
  expect_error(rlim(Hpol = P, nsamp = nsamp, burn = burn, thin = thin, seed = 456))
  
  #if there is Hpolytope object and lim is set to NULL an warning must appear
  expect_warning(samp <- rlim(lim = NULL, Hpol = P, nsamp = nsamp, burn = burn, thin = thin, 
                              seed = 456))

  # Check dimensions of sample are correct
  expect_equal(dim(samp), c(nsamp,ncol(P@A)))
  
  # Check type is correct
  expect_type(samp, "double")
  
  # Check if seed work as expected
  # same seed returns the same values 
  expect_warning(samp2 <- rlim(lim = NULL, Hpol = P, nsamp = nsamp, burn = burn, thin = thin, 
                              seed = 456))
  expect_equal(samp2,samp)
  #different seed returns a different value
  expect_warning(samp3 <- rlim(lim = NULL, Hpol = P, nsamp = nsamp, burn = burn, thin = thin, 
                               seed = 123))
  expect_true(any(samp3 != samp))
  
  
  
})














