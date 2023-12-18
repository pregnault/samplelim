# Tests on the behaviour of rlim ----

#### Using rlim() for sampling LIM ####

test_that("rlim() works as expected for LIM object", {
  # lim is mandatory
  expect_error(rlim())
  # Create a LIM object from a Declaration File
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  BOWF <- df2lim(DF)
  nsamp = 20
  burn = 0
  thin = 1
  samp <- rlim(lim = BOWF, nsamp = nsamp, burn = burn, thin = thin,
               seed = 456)
  # Check dimensions of sample are correct
  expect_equal(dim(samp), c(nsamp,ncol(BOWF$A)))
})
