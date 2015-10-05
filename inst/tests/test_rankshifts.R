context("rankshifts")

test_that("rankshifts loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))

    library(codyn)

    # Load our example data set
    knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))

    # Basic test if mean rank produces data frame with right structure and values
    result <- meanrank(knz_001d)
    expect_that(length(names(result)), equals(2))
    expect_that(names(result)[[1]], matches("year_pair"))
    expect_that(class(result[,1]), matches("factor"))
    expect_that(names(result)[[2]], matches("MRS"))
    expect_that(class(result[,2]), matches("numeric"))
    expect_that(nrow(result), equals(23))
    expect_that(result[[1,2]], equals(332.9921, tolerance=.00001))
    expect_that(result[[23,2]], equals(670.6388, tolerance=.00001))

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
  	myresults<-meanrank(dat1, "year", "species", "abundance")
  	expect_that(class(myresults[,2]), equals("numeric"))
  	expect_that(length(myresults), equals(2))
  	#test that meanrank function works with different column names
  	myresults2<-meanrank(dat2,  "yr", "sp", "abund")
  	expect_that(myresults2, equals(myresults))
  	#test that gives a warning if running on factor instead of numeric
  	expect_error(meanrank(dat2, "yr", "sp", "subplot"))

		#test the meanrankshift function
  	#test that works on a single replicate
  	myresults3<-meanrankshift(dat1, replicate=NA, time.var="year", species="species", abundance="abundance")
  	expect_that(myresults3, equals(myresults2))

  	#test that works whether replicate is a character or factor
  	myresults4<-meanrankshift(dat3, replicate="subplot", time.var="year", species="species", abundance="abundance")
  	expect_that((myresults4[2,2]), equals(myresults3[2,2]))

  	#test that works with multiple replicates
  	myresults5<-meanrankshift(knz_001d, replicate="subplot", time.var="year", species="species", abundance="abundance")
  	expect_that(myresults4[2,2], equals(myresults5[2,2]))

    #test that works with different column names
  	myresults6<-meanrankshift(knz_001d2, replicate="sub", time.var="yr", species="sp", abundance="abund")
  	expect_that(myresults6[2,2], equals(myresults5[2,2]))

		#test that works regardless of whether parameter is specified or just ordered
  	myresults7<-meanrankshift(knz_001d, "subplot", "year", "species", "abundance")
  	expect_that(myresults7, is_identical_to(myresults5))

		#test that works with different column orders if names specified
  	myresults8<-meanrankshift(knz_001d, abundance="abundance", replicate="subplot", species="species", time.var="year")
  	expect_that(myresults8, is_identical_to(myresults5))

  	#test that works with different column names
  	myresults9<-meanrankshift(knz_001d2, replicate="sub", time.var="yr", species="sp", abundance="abund")
  	expect_that(myresults9[2,2], equals(myresults5[2,2]))

		#test that it works even if there are additional unused columns
  	knz_001d3<-knz_001d
  	knz_001d3$site<-"KNZ"
  	myresults10<-meanrankshift(knz_001d3, "subplot", "year", "species", "abundance")
  	expect_that(myresults10, is_identical_to(myresults5))

  	#test that gives error when abundance column is a character or factor
#   	expect_error(meanrankshift(knz_001d2, replicate="sub", time.var="yr", species="sp", abundance="randcharacter"))
#   	expect_error(meanrankshift(knz_001d2, replicate="sub", time.var="yr", species="sp", abundance="randfactor"))

})