context("centroid_change")


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


test_that("centroid_change function returns correct result", {
  
  #test the returned result with default setting, no treatment
  myresults1 <- centroid_change(dat1, time.var = "year",
                           abundance.var = "relative_cover",
                           species.var = "species",
                           replicate.var = "plot")
  
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 3)
  expect_equal(ncol(myresults1), 3)##this might need to be changed
  expect_equal(myresults1$centroid_distance_change[1], 0.1325553, tolerance = 0.00001)
  expect_equal(myresults1$dispersion_change[1], 0.03804344, tolerance = 0.000001)

  #test that it works with treatment
  myresults2 <- centroid_change(pplots, abundance.var = "relative_cover",
                                    replicate.var = "plot",
                                    species.var = "species",
                                    time.var = "year",
                                    treatment.var = "treatment")
  
  expect_equal(nrow(myresults2), 9)
  expect_equal(ncol(myresults2), 4) ##this might need to be changed
  
  #test that is doesn't work with missing abundance
  expect_error(centroid_change(bdat, abundance.var = "relative_cover",
                               replicate.var = "plot",
                               species.var = "species",
                               time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(centroid_change(bdat2, abundance.var = "relative_cover",
                               replicate.var = "plot",
                               species.var = "species",
                               time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")

})
