setwd("C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\For R package")

library(vegan)
library(tidyverse)

pdata<-read.csv('pplots_example.csv')

wide<-pdata%>%
  spread(species, abundance, fill=0)

div<-diversity(wide[,5:102], index = "invsimpson")
rich<-specnumber(wide[,5:102])
info<-wide[,1:4]
info$simpeven<-div/rich


x<-c(80,40,20,10,1)
x1<-c(80,40,20,10,0.5)

Evar(x)
Evar(x1)


EQ(x)
EQ(x1)


SimpsonEvenness(x)
SimpsonEvenness(x1)

y<-subset(pplots, plot==25&year==2002)
y1<-y$relative_cover

SimpsonEvenness(y1)
EQ(y1)
Evar(y1)

evar<-community_structure(pplots, 
                     time.var="year", 
                     abundance.var = "relative_cover",
                     metric = "Evar", 
                     replicate.var = "plot")
eq<-community_structure(pplots, 
                          time.var="year", 
                          abundance.var = "relative_cover",
                          metric = "EQ", 
                        replicate.var = "plot")
se<-community_structure(pplots, 
                          time.var="year", 
                          abundance.var = "relative_cover",
                          metric = "SimpsonEvenness", 
                        replicate.var = "plot")

merge1<-merge(evar, eq, by=c("year", "plot", "richness"))
merge<-merge(merge1, se, by=c('year', 'plot', 'richness'))
