context("temporal_torus_translation")

test_that("temporal_torus_translation loads and returns correct result", {
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
  datmat<-calComDat(dat1, "species", "year", "abundance")
  
  

  
  #Test the genRand function
  myresults<-genRand(datmat)
  expect_true(is.matrix(myresults))
  #test that does not generate the same matrix each time
  myresults2<-genRand(datmat)
  myresults3<-sort(as.vector(myresults==myresults2))
  expect_that(unique(myresults3)[1], equals(FALSE))

  #test that random matrix has same dimensions as datmat
  expect_that(dim(myresults), equals(dim(datmat)))
  
  #test that myresults has the same values as datmat 
  myresults4<-sort(as.matrix(myresults[,1]))
  datmat2<-sort(as.matrix(datmat[,1]))
  expect_that(myresults4, equals(datmat2)) 
  
  
  # Test the temporal_torus_translation function
  myresults<-temporal_torus_translation(dat1, "species", "year", "abundance", calVR)
  #test that returns a single numeric value
  expect_that(length(myresults), equals(1))
  expect_true(is.numeric(myresults))  
  #test that does not generate the same value every time
  myresults2<-temporal_torus_translation(dat1, "species", "year", "abundance", calVR)
  expect_that(myresults==myresults2, equals(FALSE))
  #test that is not sensitive to different column names
  myresults3<-temporal_torus_translation(dat2, "sp", "yr", "abund", calVR)
  expect_that(length(myresults3), equals(1))
  expect_true(is.numeric(myresults3)) 
  
  
  
  #Test the temporal_torus_translation_CI function
  myresults<-temporal_torus_translation_CI(dat1, replicate=NA, species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)
  
  ##TO DO: Decide how to handle replicates that are factors with missing levels
  #For example, this:
  myresults2<-temporal_torus_translation_CI(dat1, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)
  #Versus this:
  dat5<-dat4
  dat5$subplot<-as.character(dat5$subplot)
  myresults3<-temporal_torus_translation_CI(dat5, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=TRUE)
  
  #TO DO: Test both options for "averagereps"
  myresults4<-temporal_torus_translation_CI(knz_001d, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=TRUE)
  myresults5<-temporal_torus_translation_CI(knz_001d, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)

  
})