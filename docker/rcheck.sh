#!/bin/bash

# install R package dependencies 
Rscript --vanilla -e " \
  options(repos = c(CRAN = 'https://cran.rstudio.com')); \
  devtools::install_dev_deps(\"$1\")"

#  prepare and execute an R CMD check
xvfb-run R CMD build $1
pkg=$(ls $1*.gz)
xvfb-run R CMD check --as-cran $pkg
rm $pkg
echo "DONE"
