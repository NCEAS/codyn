context("rate_change")

test_that("rate_change loads and returns correct result", {
	# Ensure that trivial tests work correctly
  expect_that(length("a"), equals(1))

	library(codyn)

	# Load our example data set
	# data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
  knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
  expect_that(names(knz_001d)[4], equals("abundance"))

	#give new column names
  knz_001d2 <- knz_001d
  names(knz_001d2)=c("sp", "yr", "sub", "abund")

  #add a random character and factor column
  knz_001d2$randcharacter<-"rchar"
  knz_001d2$randfactor<-as.factor(knz_001d2$randcharacter)

  #take a subset
  dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")

  #rename the subset
  dat2<-dat1
  names(dat2)=c("sp", "yr", "sub", "abund")

  #make subplot a character
  dat3<-dat1
  dat3$subplot<-as.character(dat3$subplot)

  #test the get_slope function
  myresults<-get_slope(dat1, "year", "species", "abundance")
  expect_that(class(myresults), equals("numeric"))
  expect_that(length(myresults), equals(1))
  #test that works with different column names
  myresults2<-get_slope(dat2,  "yr", "sp", "abund")
  expect_that(myresults, equals(myresults2))
  #myresultsNA<-get_slope(dat2, "yr", "sp", "subplot")
  #test that gives a warning if running on factor instead of numeric
  expect_error(get_slope(dat2, "yr", "sp", "subplot"))



	#test the rate_change function
  #test that works on a single replicate
  myresults3<-rate_change(dat1, replicate.var=NA, time.var="year", species.var="species", abundance.var="abundance")
  expect_that(myresults3, equals(myresults2))

  #test that will still run if there are missing levels in a factor "replicate"; deleting levels that are NaN
  myresults4<-rate_change(dat1, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
  #this will give a warning because replicate is a factor without all values present in dat1 - the warning is a good thing
  myresults5<-as.numeric(myresults4[2])
  expect_that(myresults5, equals(myresults3))

  #test that works whether replicate is a character or factor
  myresults6<-rate_change(dat3, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
  expect_that((myresults6[1,2]), equals(myresults3))

  #test that works with multiple replicates
  myresults7<-rate_change(knz_001d, replicate.var="subplot", time.var="year", species.var="species", abundance.var="abundance")
  expect_that(myresults6[1,2], equals(myresults7[1,2]))

    #test that works with different column names
  myresults8<-rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="abund")
  expect_that(myresults7[1,2], equals(myresults8[1,2]))

	#test that works regardless of whether parameter is specified or just ordered
  myresults8<-rate_change(knz_001d, replicate.var="subplot", time.var = "year", species.var = "species", abundance.var = "abundance")
  expect_that(myresults8, is_identical_to(myresults7))

	#test that works with different column orders if names specified
  myresults9<-rate_change(knz_001d, abundance.var="abundance", replicate.var="subplot", species.var="species", time.var="year")
  expect_that(myresults9, is_identical_to(myresults7))

  #test that works with different column names
  myresults10<-rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="abund")
  expect_that(myresults10[1,2], equals(myresults7[1,2]))

	#test that it works even if there are additional unused columns
  knz_001d3<-knz_001d
  knz_001d3$site<-"KNZ"
  myresults11<-rate_change(knz_001d3, replicate.var = "subplot", time.var = "year", species.var = "species", abundance.var = "abundance")
  expect_that(myresults11, is_identical_to(myresults7))

  #test that gives error when abundance column is a character or factor
  expect_error(rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="randcharacter"))
  expect_error(rate_change(knz_001d2, replicate.var="sub", time.var="yr", species.var="sp", abundance.var="randfactor"))

})