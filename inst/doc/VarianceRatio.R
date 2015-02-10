## ----message=FALSE-------------------------------------------------------
library(codyn)
library(knitr)
library(dplyr)

## ----results='asis'------------------------------------------------------
knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
kable(head(knz_001d))

## ----results='asis'------------------------------------------------------
myresults<-varianceratio(knz_001d, "subplot", "species", "year", "abundance", 1)
kable(myresults)

## ----results='asis'------------------------------------------------------
one_matrix <- knz_001d %>%
  filter(subplot == "A_1") %>%
  select(-subplot) %>%
  codyn:::calComDat("species", "year", "abundance")

kable(one_matrix[1:5, 1:5])

## ------------------------------------------------------------------------
one_ts <- knz_001d %>%
  filter(subplot == "A_1") %>%
  select(-subplot) %>%
  codyn:::calComTS("species", "year", "abundance")

kable(one_ts[1:5, 1:5])

