context("RAC_difference")


# prepare data ----------------------------------------


# Load our example dataset
data("pplots", package = "codyn") 

#take a subset without time 
dat1 <- subset(pplots, plot %in% c(6,25) & year ==2002)

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


test_that("RAC_difference function returns correct result", {
  
  #test the returned result with default setting and no blocking, pooling or time
  myresults1 <- RAC_difference(dat1, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species")
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 1)
  expect_equal(ncol(myresults1), 6)##this might need to be changes
  expect_equal(myresults1$richness_diff, 0.1923077, tolerance = 0.00001)
  expect_equal(myresults1$evenness_diff, 0.0002433149, tolerance = 0.000000001)
  expect_equal(myresults1$rank_diff, 0.1449704, tolerance = 0.00001)
  expect_equal(myresults1$species_diff, 0.4230769, tolerance = 0.00001)

  #test that it works with replicates
  myresults2 <- RAC_difference(pplots, abundance.var = "relative_cover",
                               replicate.var = "plot",
                               species.var = "species",
                               time.var = "year")
  
  expect_equal(nrow(myresults2), 54)
  expect_equal(ncol(myresults2), 8) ##this might need to be changed
  
  #test the returned result with blocking and time
  myresults2 <- RAC_difference(pplots, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species",
                               block.var = "block",
                               treatment.var = "treatment",
                               time.var = "year")
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults2), 1)
  expect_equal(ncol(myresults2), 6)##this might need to be changes
  expect_equal(myresults2$richness_diff, 0.1923077, tolerance = 0.00001)
  expect_equal(myresults2$evenness_diff, 0.0002433149, tolerance = 0.000000001)
  expect_equal(myresults2$rank_diff, 0.1449704, tolerance = 0.00001)
  expect_equal(myresults2$species_diff, 0.4230769, tolerance = 0.00001)
  
  #test the returned result with pooling and time
  myresults1 <- RAC_difference(dat1, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species")
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 1)
  expect_equal(ncol(myresults1), 6)##this might need to be changes
  expect_equal(myresults1$richness_diff, 0.1923077, tolerance = 0.00001)
  expect_equal(myresults1$evenness_diff, 0.0002433149, tolerance = 0.000000001)
  expect_equal(myresults1$rank_diff, 0.1449704, tolerance = 0.00001)
  expect_equal(myresults1$species_diff, 0.4230769, tolerance = 0.00001)
  
  #test that is doesn't work with missing abundance
  expect_error(RAC_difference(bdat, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(RAC_difference(bdat2, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")
  
  #test that is doesn't work with missing species name
  expect_error(RAC_difference(bdat3, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "Species names are missing")
  
})
