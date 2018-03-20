context("curve_difference")


# prepare data ----------------------------------------


# Load our example dataset
data("pplots", package = "codyn") 

#take a subset without time 
dat1 <- subset(pplots, plot %in% c(6, 25) & year == 2002)

#make a subset within a block and no time
dat2 <- subset(pplots, plot %in% c(25, 29) & year == 2002)

#make a subset for 2002 with just N1P0 and N2P0
dat3 <- subset(pplots, treatment %in% c("N2P0", "N1P0") & year == 2002)

#make missing abundance
bdat <- dat1
bdat$relative_cover[1] <- NA

#repeat a species
x <- c("N1P0", 25, 1, 2002, "senecio plattensis", 0.002123142)
bdat2 <- rbind(dat1, x)

# run tests -------------------------------------------


test_that("curve_difference function returns correct result", {
  
  #test the returned result with default setting and no blocking, pooling or time
  myresults1 <- curve_difference(dat1, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species")
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 1)
  expect_equal(ncol(myresults1), 3)
  expect_equal(myresults1$curve_diff, 0.02158507, tolerance = 0.00001)
 
  #test that it works with time
  myresults2 <- curve_difference(pplots, abundance.var = "relative_cover",
                               replicate.var = "plot",
                               species.var = "species",
                               time.var = "year")
  
  expect_equal(nrow(myresults2), 612)
  expect_equal(ncol(myresults2), 4)
  
  #test that it works with time and treatment specified
  myresults2.2 <- curve_difference(pplots, abundance.var = "relative_cover",
                               replicate.var = "plot",
                               species.var = "species",
                               time.var = "year",
                               treatment.var = "treatment")
  
  expect_equal(nrow(myresults2.2), 612)
  expect_equal(ncol(myresults2.2), 6)
  
  #test the returned result with blocking and no time
  myresults3 <- curve_difference(dat2, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species",
                               block.var = "block",
                               treatment.var = "treatment")
  
  expect_is(myresults3, "data.frame")
  expect_equal(nrow(myresults3), 1)
  expect_equal(ncol(myresults3), 6)
  expect_equal(myresults3$curve_diff, 0.008164927, tolerance = 0.00001)

  #test that returned results with blocking and time
  myresults3.2 <- curve_difference(pplots, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species",
                               block.var = "block",
                               treatment.var = "treatment",
                               time.var = "year")
  
  expect_equal(nrow(myresults3.2), 72)
  expect_equal(ncol(myresults3.2), 7)
  
  #test the returned result with pooling and no time
  myresults4 <- curve_difference(dat3, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species", 
                               pool = TRUE,
                               treatment.var = "treatment")
  
  expect_is(myresults4, "data.frame")
  expect_equal(nrow(myresults4), 1)
  expect_equal(ncol(myresults4), 3)
  expect_equal(myresults4$curve_diff, 0.01284473, tolerance = 0.00001)
 
  
  #test the returned result with pooling and time
  myresults4.2 <- curve_difference(pplots, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species", 
                               pool = TRUE,
                               treatment.var = "treatment",
                               time.var = "year")
  
  expect_equal(nrow(myresults4.2), 12)
  expect_equal(ncol(myresults4.2), 4)
  
  #test that is doesn't work with missing abundance
  expect_error(curve_difference(bdat, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(curve_difference(bdat2, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")
  
})
