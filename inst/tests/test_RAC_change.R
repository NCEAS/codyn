context("RAC_change")


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

#make species name missing
bdat3 <- dat1
bdat3$species[1] <- NA

#make a plot with only 1 species
bdat4 <- subset(dat1, species %in% c("amorpha canescens", "andropogon gerardii"))
# run tests -------------------------------------------


test_that("RAC_change function returns correct result", {
  
  #test the returned result with default setting
  myresults1 <- RAC_change(dat1, time.var = "year",
                           abundance.var = "relative_cover",
                           species.var = "species")
  
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 1)
  expect_equal(ncol(myresults1), 7)
  expect_equal(myresults1$richness_change, -0.15, tolerance = 0.00001)
  expect_equal(myresults1$evenness_change, -0.05537606, tolerance = 0.00001)
  expect_equal(myresults1$rank_change, 0.145, tolerance = 0.00001)
  expect_equal(myresults1$gains, 0.1, tolerance = 0.00001)
  expect_equal(myresults1$losses, 0.25, tolerance = 0.00001)

  #test that it works with replicates
  myresults2 <- RAC_change(pplots, abundance.var = "relative_cover",
                                    replicate.var = "plot",
                                    species.var = "species",
                                    time.var = "year")
  
  expect_equal(nrow(myresults2), 54)
  expect_equal(ncol(myresults2), 8) 
  
  #test that is doesn't work with missing abundance
  expect_error(RAC_change(bdat, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(RAC_change(bdat2, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")
  
  #test that is doesn't work with missing species name
  expect_error(RAC_change(bdat3, abundance.var = "relative_cover",
                          replicate.var = "plot",
                          species.var = "species",
                          time.var = "year"), "Species names are missing")
  
  #test that give warning for evenness NA
  expect_warning(RAC_change(bdat4, abundance.var = "relative_cover",
                            replicate.var = "plot",
                            species.var = "species",
                            time.var = "year"), "Evenness_change values contain NAs because there are plots with only one species")
  
})
