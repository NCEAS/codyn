context("community_stability")

test_that("community_stability loads and returns correct result", {
  # Ensure that trivial tests work correctly
  expect_that(length("a"), equals(1))
  
  library(codyn)
  
  # Load our example data set
  data(knz_001d)
  
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
  
  #test the stability_onerep function
  dat1agg<-aggregate(abundance~year + subplot, data=dat1, sum)
  myresults<-stability_onerep(dat1agg,  "abundance")
  expect_that(class(myresults), equals("numeric"))
  expect_that(length(myresults), equals(1))
  expect_that(myresults, equals(4.1233, tolerance=.00001))
  
  
  #test that works with different column names
  dat2agg<-aggregate(abund~yr + sub, data=dat2, sum)
  myresults2<-stability_onerep(dat2agg, "abund")
  expect_that(myresults, equals(myresults2))
  myresultsNA<-stability_onerep(dat2agg, "subplot")
  #test that gives a warning if running on factor instead of numeric
  expect_that(stability_onerep(dat2agg, "subplot"), gives_warning())
  
  
  #test the community_stability function
  #test that works on a single replicate
  myresults3<-community_stability(dat1, replicate.var=NA, time.var="year", abundance.var="abundance")
  expect_that(myresults3, equals(myresults2))
  
  #test that will still run if there are missing levels in a factor "replicate"; deleting levels that are NaN
  myresults4<-community_stability(dat1, replicate.var="subplot", time.var="year", abundance.var="abundance")
  #this will give a warning because replicate is a factor without all values present in dat1 - the warning is a good thing
  myresults5<-as.numeric(myresults4[2])
  expect_that(myresults5, equals(myresults3)) 
  
  #test that works whether replicate is a character or factor
  myresults6<-community_stability(dat3, replicate.var="subplot", time.var="year", abundance.var="abundance")
  expect_that((myresults6[1,2]), equals(myresults3))  
  
  #test that works with multiple replicates
  myresults7<-community_stability(knz_001d, replicate.var="subplot", time.var="year", abundance.var="abundance")
  expect_that(myresults6[1,2], equals(myresults7[1,2]))  
  
  #test that works with different column names
  myresults8<-community_stability(knz_001d2, replicate.var="sub", time.var="yr", abundance.var="abund")
  expect_that(myresults7[1,2], equals(myresults8[1,2]))
  
  #test that gives error when abundance column is a character or factor
  expect_error(community_stability(knz_001d2, replicate.var="sub", time.var="yr", abundance.var="randcharacter"))
  expect_error(community_stability(knz_001d2, replicate.var="sub", time.var="yr", abundance.var="randfactor"))

  
  # test that works regardless of order of the input replicates
  knz_001dreorder <-knz_001d[order(knz_001d$abundance, knz_001d$year, knz_001d$species),]
  myresults<-community_stability(knz_001d, time.var="year", abundance.var="abundance", 
                             replicate.var="subplot")
  
  myresults_reorder<-community_stability(knz_001dreorder, time.var="year", abundance.var="abundance", 
                                     replicate.var="subplot")
  expect_equal(myresults, myresults_reorder)
  
  })