context("multivariate_difference")


# prepare data ----------------------------------------


# Load our example dataset
data("pplots", package = "codyn") 

#make a dataset with no time
dat1 <- subset(pplots, year == 2003)

#make missing abundance
bdat <- dat1
bdat$relative_cover[1] <- NA

#repeat a species
x <- c("N1P0", 25, 1, 2002, "senecio plattensis", 0.002123142)
bdat2 <- rbind(dat1, x)

#make species name missing
bdat3 <- dat1
bdat3$species[1] <- NA
# run tests -------------------------------------------


test_that("multivariate_difference function returns correct result", {
  
  #test the returned result with default setting, no time
  myresults1 <- multivariate_difference(dat1, abundance.var = "relative_cover",
                           species.var = "species",
                           replicate.var = "plot",
                           treatment.var = "treatment")
  
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 3)
  expect_equal(ncol(myresults1), 5)
  expect_equal(myresults1$compositon_diff[1], 0.1068927, tolerance = 0.00001)
  expect_equal(myresults1$dispersion_diff[1], 0.008186288, tolerance = 0.000001)

  #test that it works with time
  myresults2 <- multivariate_difference(pplots, abundance.var = "relative_cover",
                                    replicate.var = "plot",
                                    species.var = "species",
                                    time.var = "year",
                                    treatment.var = "treatment")
  
  expect_equal(nrow(myresults2), 12)
  expect_equal(ncol(myresults2), 6) 
  
  #test that is doesn't work with missing abundance
  expect_error(multivariate_difference(bdat, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(multivariate_difference(bdat2, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")
  
  #test that is doesn't work with missing species name
  expect_error(multivariate_difference(bdat3, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year", 
                          treatment.var = "treatment"), "Species names are missing")
  
})
