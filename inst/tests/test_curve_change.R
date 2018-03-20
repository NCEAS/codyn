context("curve_change")


# prepare data ----------------------------------------


# Load our example dataset
data("pplots", package = "codyn") 

#take a subset without replicates 
dat1 <- subset(pplots, plot ==25 & year %in% c(2002, 2003))

#make missing abundance
bdat <- dat1
bdat$relative_cover[1] <- NA

#repeat a species
x <- c("N1P0", 25, 1, 2002, "senecio plattensis", 0.002123142)
bdat2 <- rbind(dat1, x)

# run tests -------------------------------------------


test_that("curve_change function returns correct result", {
  
  #test the returned result with default setting
  myresults1 <- curve_change(dat1, time.var = "year",
                           abundance.var = "relative_cover",
                           species.var = "species")
  
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 1)
  expect_equal(ncol(myresults1), 3)
  expect_equal(myresults1$curve_change, 0.03259107, tolerance = 0.00001)


  #test that it works with replicates
  myresults2 <- curve_change(pplots, abundance.var = "relative_cover",
                                    replicate.var = "plot",
                                    species.var = "species",
                                    time.var = "year")
  
  expect_equal(nrow(myresults2), 54)
  expect_equal(ncol(myresults2), 4)
  
  #test that is doesn't work with missing abundance
  expect_error(curve_change(bdat, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(curve_change(bdat2, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")
  
})
