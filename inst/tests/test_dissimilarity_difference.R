context("dissimilarity_difference")


# prepare data ----------------------------------------


# Load our example dataset
data("pplots", package = "codyn") 

#make a dataset with no time
dat1 <- subset(pplots, year == 2002)

#make missing abundance
bdat <- dat1
bdat$relative_cover[1] <- NA

#repeat a species
x <- c("N1P0", 25, 1, 2002, "senecio plattensis", 0.002123142)
bdat2 <- rbind(dat1, x)

# run tests -------------------------------------------

test_that("dissimilarity_difference function returns correct result", {
  
  #test the returned result with default setting, no time
  myresults1 <- dissimilarity_difference(dat1, abundance.var = "relative_cover",
                                         species.var = "species",
                                         replicate.var = "plot",
                                         treatment.var = "treatment")
  
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 3)
  expect_equal(ncol(myresults1), 4)
  expect_equal(myresults1$BC_between_diff[1], 0.2807829, tolerance = 0.00001)
  expect_equal(myresults1$BC_within_diff[1], 0.01396088, tolerance = 0.000001)

  #test that it works with time
  myresults2 <- dissimilarity_difference(pplots, abundance.var = "relative_cover",
                                         replicate.var = "plot",
                                         species.var = "species",
                                         time.var = "year",
                                         treatment.var = "treatment")
  
  expect_equal(nrow(myresults2), 12)
  expect_equal(ncol(myresults2), 5) 
  
  #test that is doesn't work with missing abundance
  expect_error(dissimilarity_difference(bdat, abundance.var = "relative_cover",
                                        replicate.var = "plot",
                                        species.var = "species",
                                        time.var = "year",
                                        treatment.var = "treatment"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(dissimilarity_difference(bdat2, abundance.var = "relative_cover",
                                        replicate.var = "plot",
                                        species.var = "species",
                                        time.var = "year",
                                        treatment.var = "treatment"), "In replicate 25 there is more than one record for species at the time point 2002")
  
})
