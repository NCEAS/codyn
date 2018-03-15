context("dissimilarity_change")


# prepare data ----------------------------------------


# Load our example dataset
data("pplots", package = "codyn") 

#make a dataset with no treatments
dat1 <- subset(pplots, treatment == "N1P0")

#make missing abundance
bdat <- dat1
bdat$relative_cover[1] <- NA

#repeat a species
x <- c("N1P0", 25, 1, 2002, "senecio plattensis", 0.002123142)
bdat2 <- rbind(x, dat1)

# run tests -------------------------------------------


test_that("dissimilarity_change function returns correct result", {
  
  #test the returned result with default setting, no treatment
  myresults1 <- dissimilarity_change(dat1, time.var = "year",
                                     abundance.var = "relative_cover",
                                     species.var = "species",
                                     replicate.var = "plot")
  
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 3)
  expect_equal(ncol(myresults1), 4)
  expect_equal(myresults1$BC_between_change[1], 0.3163601, tolerance = 0.00001)
  expect_equal(myresults1$BC_within_change[1], 0.060276258, tolerance = 0.000001)

  #test that it works with treatment
  myresults2 <- dissimilarity_change(pplots, abundance.var = "relative_cover",
                                      replicate.var = "plot",
                                      species.var = "species",
                                      time.var = "year",
                                      treatment.var = "treatment")
  
  expect_equal(nrow(myresults2), 9)
  expect_equal(ncol(myresults2), 5) 
  
  #test that is doesn't work with missing abundance
  expect_error(dissimilarity_change(bdat, abundance.var = "relative_cover",
                                    replicate.var = "plot",
                                    species.var = "species",
                                    time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(dissimilarity_change(bdat2, abundance.var = "relative_cover",
                                    replicate.var = "plot",
                                    species.var = "species",
                                    time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")

})
