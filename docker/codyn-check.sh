#/bin/sh
# A shell wrapper to prepare and execute an R CMD check
cd /src/dev
R CMD build codyn
R CMD check --as-cran codyn_0.9.0.9000.tar.gz
echo "DONE"
