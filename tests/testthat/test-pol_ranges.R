# Tests on the behaviour of the pol.ranges function ----

test_that("pol.ranges returns expected outputs where each variables are well bounded",{
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  model <- df2lim(DF)
  ranges <- pol.ranges(model$A,model$B,model$G,model$H)
  
  #Check if ranges is of the double type
  expect_type(ranges,"double")
  
  #Check if number of lines is equal to the numbers of unknowns
  expect_equal(dim(ranges)[1],model$NUnknowns)
  
  #Check if number of lines is equal to the numbers of unknowns
  expect_equal(dim(ranges)[2],3)
  
  #Check if columns have the right name
  expect_equal(colnames(ranges),c("min","max","range"))
  
  #Check min and max are correctly sorted i.e. for every variable min is less or equal than max
  expect_equal(sum(ranges[,1]<ranges[,2]),model$NUnknowns)
  
  #Check that range is equal to the difference between min and max 
  expect_equal(ranges[,2]-ranges[,1],ranges[,3])
  
  #Check if values are correctly rounded (digits=8)
  rounded_ranges <- round(ranges, digits = 8)
  #if ranges is equal to rounded_ranges that means it is already rounded correctly
  expect_equal(sum(ranges==rounded_ranges),3*model$NUnknowns)
  
})

test_that("pol.ranges properly checks arguments (dim and types)",{
  
  #Ideal arguments for comparison 
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  model <- df2lim(DF)
  expected_ranges <- pol.ranges(model$A,model$B,model$G,model$H)
  
  #When A is a dataframe
  ranges <- pol.ranges(as.data.frame(model$A),model$B,model$G,model$H)
  expect_equal(ranges,expected_ranges)
  
  #When G is a dataframe
  ranges <- pol.ranges(model$A,model$B,data.frame(model$G),model$H)
  expect_equal(ranges,expected_ranges)
  
  #When B is a vector
  ranges <- pol.ranges(model$A,as.vector(model$B),model$G,model$H)
  expect_equal(ranges,expected_ranges)
  
  #When H is a vector
  ranges <- pol.ranges(model$A,model$B,model$G,as.vector(model$H))
  expect_equal(ranges,expected_ranges)
  
  
  #Check if code stops if G is null
  expect_error(pol.ranges(model$A,model$B,NULL,model$H),"G is NULL")
  
  #Check if code stops if G and H have incompatible dimensions
  expect_error(pol.ranges(model$A,model$B,model$G,model$H[2:length(model$H)]),"G and H have incompatible dimensions.")
  
  #Check if code stops if A is not NULL but B is
  expect_error(pol.ranges(model$A,NULL,model$G,model$H),"A and B have incompatible dimensions.")
  
  #Check if code stops if A and B have incompatible dimensions
  expect_error(pol.ranges(model$A,model$B[2:length(model$B)],model$G,model$H),"A and B have incompatible dimensions.")
  
  #Check if code stops if A and G have different number of variables
  expect_error(pol.ranges(model$A,model$B,model$G[,2:ncol(model$G)],model$H),"A and G have incompatible dimensions.")
  
})

test_that("pol.ranges function handles variables with negative values correctly",{
  #A model where a variable has at least one non finite bounds
  G <- -matrix(c(1,0,
                 0, -0.5,
                 0, 1,
                 -1, 0),
               byrow = TRUE,
               nrow = 4, ncol = 2)
  H <- -matrix(c(-5, 3, 12, 30), nrow = 4)
  ranges <- pol.ranges(G = G, H = H)
  
  #Check if ranges is of the double type
  expect_type(ranges,"double")
  
  #Check if number of lines is equal to the numbers of unknowns
  expect_equal(dim(ranges)[1],ncol(G))
  
  #Check if number of lines is equal to the numbers of unknowns
  expect_equal(dim(ranges)[2],3)
  
  #Check if columns have the right name
  expect_equal(colnames(ranges),c("min","max","range"))
  
  #Check min and max are correctly sorted i.e. for every variable min is less or equal than max
  expect_equal(sum(ranges[,1]<ranges[,2]),ncol(G))
  
  #Check that range is equal to the difference between min and max 
  expect_equal(ranges[,2]-ranges[,1],ranges[,3])
  
  #Check if values are correctly rounded (digits=8)
  rounded_ranges <- round(ranges, digits = 8)
  #if ranges is equal to rounded_ranges that means it is already rounded correctly
  expect_equal(sum(ranges==rounded_ranges),3*ncol(G))
  
})

test_that("pol.ranges warns if a variable has an empty or a non-finite domain",{
  #A model where a variable has at least one non finite bounds
  G <- -matrix(c(-1,0,
                 0, -0.5),
               byrow = TRUE,
               nrow = 2, ncol = 2)
  H <- -matrix(c(0, 3), nrow = 2)
  expect_warning(pol.ranges(G=G,H=H),"non-finite")

  #A model where a variable has an empty domain
  G <- -matrix(c(1, 2,
                 -1, -2),
               byrow = TRUE,
               nrow = 2, ncol = 2)
  H <- -matrix(c(4, -10), nrow = 2)
  expect_warning(pol.ranges(G=G,H=H),"empty")
  
})


# Tests on the behaviour of the lim.ranges function ----

test_that("pol.ranges returns expected outputs",{
  DF <- system.file("extdata", "DeclarationFileBOWF-short.txt", package = "samplelim")
  model <- df2lim(DF)
  expected_ranges <- pol.ranges(model$A,model$B,model$G,model$H)
  rownames(expected_ranges) <-model$Unknowns

  ranges <- lim.ranges(model)
  expect_equal(ranges,expected_ranges)
  
})