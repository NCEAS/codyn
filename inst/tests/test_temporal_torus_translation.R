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
  
  #Test that returns an error when abundance is a character or factor column
  expect_error(temporal_torus_translation_CI(knz_001d2, replicate="sub", species="sp", year="yr", abundance="randcharacter", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE))
  expect_error(temporal_torus_translation_CI(knz_001d2, replicate="sub", species="sp", year="yr", abundance="randfactor", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE))
  
  
  
  ##Test that missing levels in a factor replicate are dropped (so that the same number of subplots if run as a character
  #or as a factor with additional, missing levels)
  
  #For example that this:
  myresults2<-temporal_torus_translation_CI(dat1, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)
  expect_that(is.factor(myresults2$subplot), equals(TRUE))
  
  myresults2b<-temporal_torus_translation_CI(dat3, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)
  expect_that(is.character(myresults2b$subplot), equals(TRUE))
  
  expect_that(as.character(myresults2$subplot), equals(myresults2b$subplot))
  
  myresults2b<-temporal_torus_translation_CI(knz_001d, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)
  
  #Test that is correct for whether "averagereps" is true or false
  dat5<-dat4
  dat5$subplot<-as.character(dat5$subplot)
  myresults3<-temporal_torus_translation_CI(dat5, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=TRUE)
  myresults4<-temporal_torus_translation_CI(dat5, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)
 
  #expect that subplot names are the same between input and output
  expect_that(unique(dat5$subplot), equals(myresults4$subplot))
  

  
  #Test both options for "averagereps"
  myresults4<-temporal_torus_translation_CI(knz_001d, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=TRUE)
  expect_that(nrow(myresults4), equals(1))
  myresults5<-temporal_torus_translation_CI(knz_001d, replicate="subplot", species="species", year="year", abundance="abundance", FUN=calVR, bootnumber=2, li=0.025, ui=0.975, averagereps=FALSE)
  expect_that(nrow(myresults5), equals(length(unique(knz_001d$subplot))))
 
})
