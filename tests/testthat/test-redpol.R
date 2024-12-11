
# Tests on the behaviour of the lim.redpol function ----

test_that("lim.redpol returns error when lim object doesn't have equalities",{
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  model <- df2lim(DF)
  model$A <- NULL
  
  #Check for error for the 2 possible states of test
  expect_error(lim.redpol(model,test=TRUE),"no equalities found")
  expect_error(lim.redpol(model,test=FALSE),"no equalities found")
})

test_that("lim.redpol returns expected projection",{
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  model <- df2lim(DF)
  red <- lim.redpol(model)
  
  #Tests on the type of the outputs
  expect_type(red$G,"double")
  expect_type(red$H,"double")
  expect_type(red$x0,"double")
  expect_type(red$Z,"double")
  
  #Tests on the dimensions of G' (reduced)
  rank=qr(model$A)$rank
  expect_equal(dim(red$G)[1],dim(model$G)[1])
  expect_equal(dim(red$G)[2],dim(model$G)[2]-rank)
  
  #Test on the length of H' (reduced)
  expect_equal(length(red$H),length(model$H))
  
  #Test on the length of x0
  expect_equal(length(red$x0),dim(model$G)[2])
  
  #Test on the dimensions of Z
  expect_equal(dim(red$Z)[1],dim(model$G)[2])
  expect_equal(dim(red$Z)[2],dim(red$G)[2])
  
  #Tests that matrices G' and H' respect decomposition properties
  expect_equal(red$G,model$G %*% red$Z)
  expect_equal(red$H,as.numeric(model$H-model$G %*% red$x0))
})

# test_that("lim.redpol can find inequalities when specified",{
#   
#   A <- matrix(c(1,1), nrow = 1, ncol = 2)
#   B <- 1
#   G <- -matrix(c(1, 0,
#                  0, 1,
#                  -1, 0),
#                byrow = TRUE,
#                nrow = 3, ncol = 2)
#   H <- -matrix(c(2, 12, -2), nrow = 3)
# 
#   lim_exm <- list(A = A, B = B,G = G,H = H)
#   lim.redpol(lim_exm)
# 
# 
# })

# Tests on the behaviour of the red2full function ----

test_that("red2full returns expected values and is reversible",{
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  model <- df2lim(DF)
  red <- lim.redpol(model)
  x0 <- red$x0
  Z <- red$Z
  
  x_red1 <- rlim(lim=red, nsamp = 500)
  x_full <- red2full(x_red1, x0, Z)
  
  expect_type(x_red1,"double")
  expect_equal(dim(x_full)[1], dim(x_red1)[1])
  expect_gte(dim(x_full)[2], dim(x_red1)[2])
  
  #revert
  x_red2 <- full2red(x_full, x0, Z)
  expect_equal(dim(x_red2)[1], dim(x_full)[1])
  expect_lte(dim(x_red2)[2], dim(x_full)[2])
  
  expect_equal(x_red1, x_red2)

})


# Tests on the behaviour of the full2red function ----

test_that("full2red returns expected values and is reversible",{
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  model <- df2lim(DF)
  
  x_full1 <- rlim(model, nsamp = 500)
  red <- lim.redpol(model)
  x_red <- full2red(x_full1, red$x0, red$Z)
  
  expect_type(x_red,"double")
  expect_equal(dim(x_red)[1], dim(x_full1)[1])
  expect_lte(dim(x_red)[2], dim(x_full1)[2])
  
  #revert 
  x_full2 <- red2full(x_red, red$x0, red$Z)
  expect_equal(dim(x_full2)[1], dim(x_red)[1])
  expect_gte(dim(x_full2)[2], dim(x_red)[2])
  
  #expect_equal also checks colnames
  colnames(x_full2) <- colnames(x_full1)
  
  expect_equal(x_full1, x_full2)
})