context("community_diversity")


# prepare data ----------------------------------------


# Load our example dataset
data("pplots", package = "codyn") 

#take a subset
#no plots or years
dat1 <- subset(pplots, plot == 25 & year == 2002)

#with time
dat2 <- subset(pplots, plot ==25 & year %in% c(2002, 2003))

#with replicates
dat3 <- subset(pplots, plot %in% c(6, 25) & year ==2002)

#make missing abundance
bdat <- dat1

bdat$relative_cover[1] <- NA

# run tests -------------------------------------------


test_that("community_diversity function returns correct result", {
  
  #test the EQ returned with default setting
  myresults1 <- community_diversity(dat1, abundance.var = "relative_cover")
  
  
  expect_is(myresults1, "data.frame")
  expect_equal(nrow(myresults1), 1)
  expect_equal(ncol(myresults1), 1)
  expect_equal(myresults1$Shannon, 1.842618, tolerance = 0.00001)

  
  #test the Evar returned when specified
  myresults2 <- community_diversity(dat1, abundance.var = "relative_cover",
                                    metric = "InverseSimpson")
  
  expect_equal(myresults2$InverseSimpson, 4.072345, tolerance = 0.00001)
  
  #test that it works with time
  myresults3 <- community_diversity(dat2, abundance.var = "relative_cover",
                                    time.var = "year")
    
  expect_equal(nrow(myresults3), 2)
  expect_equal(ncol(myresults3), 2)
  
  #test that it works with replicates
  myresults4 <- community_diversity(dat3, abundance.var = "relative_cover",
                                    replicate.var = "plot")
  
  expect_equal(nrow(myresults4), 2)
  expect_equal(ncol(myresults4), 2)
  
  #test that it works with both, time and replicates
  myresults5 <- community_diversity(pplots, abundance.var = "relative_cover",
                                    time.var = "year",
                                    replicate.var = "plot")
  
  expect_equal(nrow(myresults5), 72)
  expect_equal(ncol(myresults5), 3)

  #test that is doesn't work with missing abundance
  expect_error(community_diversity(bdat, abundance.var = "relative_cover"), "Abundance column contains missing values")
  
})
