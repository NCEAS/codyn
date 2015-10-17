context("temporal_torus_translation")

test_that("temporal_torus_translation loads and returns correct result", {
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
  knz_001d2$randcharacter<-sample(letters, size = nrow(knz_001d2), replace = T)
  knz_001d2$randfactor<-as.factor(knz_001d2$randcharacter)
  
  #take a subset
  dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")
  
  #rename the subset
  dat2<-dat1
  names(dat2)=c("sp", "yr", "sub", "abund")
  
  #make subplot a character
  dat3<-dat1
  dat3$subplot<-as.character(dat3$subplot)
  
  #take a two replicate subset
  dat4<-subset(knz_001d, subplot=="A_1"|subplot=="A_2")
  
  #make a species matrix
  
  datmat<-transpose_community(dat1, "year","species",  "abundance")
  
  #Test the genRand function
  myresults<-genRand(datmat)
  expect_true(is.matrix(myresults))
  
  #test that does not generate the same matrix each time
  myresults2<-genRand(datmat)
  myresults3<-sort(as.vector(myresults==myresults2))
  
  expect_true(any(as.vector(myresults==myresults2) == FALSE))
  
  #test that random matrix has same dimensions as datmat
  expect_equal(dim(myresults), dim(datmat))
  
  #test that myresults has the same values as datmat, if both are sorted 
  rand.col.1 <-sort(as.matrix(myresults[,1]))
  dat.col.1 <-sort(as.matrix(datmat[,1]))
  expect_equal(rand.col.1, dat.col.1) 
  
  
  # Test the temporal_torus_translation function
  myresults <- temporal_torus_translation(dat1, "year", "species", "abundance", variance_ratio_matrixdata)
  #test that returns a single numeric value
  expect_equal(length(myresults), 1)
  expect_true(is.numeric(myresults))  
  #test that does not generate the same value every time
  myresults2<-temporal_torus_translation(dat1, "year","species",  "abundance", variance_ratio_matrixdata)
  expect_that(myresults==myresults2, equals(FALSE))
  
  #test that is not sensitive to different column names
  myresults3 <- temporal_torus_translation(dat2, "yr",  "sp", "abund", variance_ratio_matrixdata)
  expect_equal(length(myresults3), 1)
  expect_true(is.numeric(myresults3)) 
  
  
  #Test the temporal_torus_translation_CI function
  myresults<-temporal_torus_translation_CI(dat1, replicate.var=NA, species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  
  #Test that returns an error when abundance is a character or factor column
  expect_error(temporal_torus_translation_CI(knz_001d2, replicate.var="sub", species.var="sp", time.var="yr", abundance.var="randcharacter", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE))

    expect_error(temporal_torus_translation_CI(knz_001d2, replicate.var="sub", species.var="sp", time.var="yr", abundance.var="randfactor", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE))
  
  
  
  ##Test that missing levels in a factor replicate are dropped (so that the same number of subplots if run as a character
  #or as a factor with additional, missing levels)
  
  #For example :
  myresults2<-temporal_torus_translation_CI(dat1, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  expect_true(is.factor(myresults2$subplot))
  
  myresults2b<-temporal_torus_translation_CI(dat3, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  expect_true(is.character(myresults2b$subplot))
  
  expect_equal(as.character(myresults2$subplot), myresults2b$subplot)
  
  myresults2b<-temporal_torus_translation_CI(knz_001d, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  
  #Test that is correct for whether "averagereps" is true or false
  dat5<-dat4
  dat5$subplot<-as.character(dat5$subplot)
  myresults3<-temporal_torus_translation_CI(dat5, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=TRUE)
  myresults4<-temporal_torus_translation_CI(dat5, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
 
  #expect that subplot names are the same between input and output
  expect_equal(unique(dat5$subplot), myresults4$subplot)
  
  #Test both options for "averagereps"
  myresults4<-temporal_torus_translation_CI(knz_001d, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=TRUE)
  expect_that(nrow(myresults4), equals(1))
  myresults5<-temporal_torus_translation_CI(knz_001d, replicate.var="subplot", species.var="species", time.var="year", abundance.var="abundance", FUN=variance_ratio_matrixdata, bootnumber=2, li=0.025, ui=0.975, average.replicates=FALSE)
  expect_that(nrow(myresults5), equals(length(unique(knz_001d$subplot))))
 
})
