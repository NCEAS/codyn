context("abundance_difference")


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

#make species name missing
bdat3 <- dat1
bdat3$species[1] <- NA
# run tests -------------------------------------------


test_that("abundance_difference function returns correct result", {
  
  #test the returned result with default setting and no blocking, pooling or time
  myresults1 <- abundance_difference(dat1, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species")
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 26)
  expect_equal(ncol(myresults1), 4)
  expect_equal(myresults1$difference[1], 0.003898635, tolerance = 0.00001)

  #test that it works with time
  myresults2 <- abundance_difference(pplots, abundance.var = "relative_cover",
                               replicate.var = "plot",
                               species.var = "species",
                               time.var = "year")
  
  expect_equal(nrow(myresults2), 612)
  expect_equal(ncol(myresults2), 7)
  
  #test that it works with time and treatment specified
  myresults2.2 <- abundance_difference(pplots, abundance.var = "relative_cover",
                                     replicate.var = "plot",
                                     species.var = "species",
                                     time.var = "year", 
                                     treatment.var = "treatment")
  
  expect_equal(nrow(myresults2.2), 612)
  expect_equal(ncol(myresults2.2), 9)
  
  #test the returned result with blocking and no time
  myresults3 <- abundance_difference(dat2, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species",
                               block.var = "block",
                               treatment.var = "treatment")
  
  expect_is(myresults3, "data.frame")
  expect_equal(nrow(myresults3), 22)
  expect_equal(ncol(myresults3), 8)
  expect_equal(myresults3$difference[1], 0.009784736, tolerance = 0.00001)
  
  #test that is works with blocking and time
  myresults3.5 <- abundance_difference(pplots, replicate.var = "plot",
                                     abundance.var = "relative_cover",
                                     species.var = "species",
                                     block.var = "block",
                                     treatment.var = "treatment",
                                     time.var = "year")
  expect_equal(nrow(myresults3.5), 22)
  expect_equal(ncol(myresults3.5), 9)

  #test the returned result with pooling and no time
  myresults4 <- abundance_difference(dat3, replicate.var = "plot",
                               abundance.var = "relative_cover",
                               species.var = "species", 
                               pool = TRUE,
                               treatment.var = "treatment")
  
  expect_is(myresults4, "data.frame")
  expect_equal(nrow(myresults4), 41)
  expect_equal(ncol(myresults4), 4)
  expect_equal(myresults4$difference[1], 0.0010101010, tolerance = 0.00001)

  #test the returned result with pooling and time
  myresults4.5 <- abundance_difference(pplots, replicate.var = "plot",
                                     abundance.var = "relative_cover",
                                     species.var = "species", 
                                     pool = TRUE,
                                     treatment.var = "treatment",
                                     time.var = "year")
  
  expect_is(myresults4.5, "data.frame")
  expect_equal(nrow(myresults4.5), 41)
  expect_equal(ncol(myresults4.5), 4)

  #test that is doesn't work with missing abundance
  expect_error(abundance_difference(bdat, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "Abundance column contains missing values")
  
  #test that is doesn't work with a repeated species
  expect_error(abundance_difference(bdat2, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "In replicate 25 there is more than one record for species at the time point 2002")
  
  #test that is doesn't work with missing species name
  expect_error(abundance_difference(bdat3, abundance.var = "relative_cover",
                              replicate.var = "plot",
                              species.var = "species",
                              time.var = "year"), "Species names are missing")
  
})
