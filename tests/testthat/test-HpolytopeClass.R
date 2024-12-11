# Tests on the behaviour of the Hpolytope class ----

test_that("The Hpolytope class is instantiated correctly with default values", {
  p_default <- Hpolytope()
  
  # Check that the object is of the correct class
  expect_s4_class(p_default, "Hpolytope")
  
  # Check the default values of the slots
  expect_true(is.matrix(p_default@A))
  expect_equal(p_default@A, as.matrix(0))
  
  expect_true(is.numeric(p_default@b))
  expect_equal(p_default@b, numeric(0))
  
  expect_true(is.numeric(p_default@volume))
  expect_equal(p_default@volume, as.numeric(NaN))
  
  expect_true(is.character(p_default@type))
  expect_equal(p_default@type, "Hpolytope")
})

test_that("The Hpolytope class is instantiated correctly with specified correct values", {
  A <- matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
  b <- c(0,0,1)
  volume <- 123.45
  type <- "CustomPolytope"
  
  p_custom <- Hpolytope(A = A, b = b, volume = volume, type = type)
  
  # Check that the object is of the correct type
  expect_s4_class(p_custom, "Hpolytope")
  
  # Check the specified values of the slots
  expect_equal(p_custom@A, A)
  expect_equal(p_custom@b, b)
  expect_equal(p_custom@volume, volume)
  expect_equal(p_custom@type, type)
})

test_that("Slot values can be changed and acessed", {
  
  A <- matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
  b <- c(0,0,1)
  volume <- 123.45
  type <- "CustomPolytope"
  
  p_custom <- Hpolytope(A = A, b = b, volume = volume, type = type)
  
  new_A <- rbind(p_custom@A, c(-1,1))
  new_b <- append(p_custom@b, c(2))
  new_volume <- 100
  new_type <- "MoreCustomPolytope"
    
  p_custom@A <- new_A
  p_custom@b <- new_b
  p_custom@volume <- new_volume
  p_custom@type <- new_type
  
  # Check the modified values of the slots
  expect_equal(p_custom@A, new_A)
  expect_equal(p_custom@b, new_b)
  expect_equal(p_custom@volume, new_volume)
  expect_equal(p_custom@type, new_type)
  
})

test_that("Errors are raised for incorrect argument types", {
  expect_error(Hpolytope(A = "not a matrix"), 'invalid object for slot "A"')
  expect_error(Hpolytope(b = "not numeric"), 'invalid object for slot "b"')
  expect_error(Hpolytope(volume = "not numeric"), 'invalid object for slot "volume"')
  expect_error(Hpolytope(type = 123), 'invalid object for slot "type"')
})

test_that("Errors are raised when slot are changed with incorrect argument types", {
  p_default <- Hpolytope()
  
  expect_error(p_default@A <- "not a matrix", "matrix")
  expect_error(p_default@b <- "not numeric", "numeric")
  expect_error(p_default@volume <- "not numeric", "numeric")
  expect_error(p_default@type <- 123, "character")
})
