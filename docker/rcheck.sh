#/bin/sh
# A shell wrapper to prepare and execute an R CMD check
srcdir=$(ls)
xvfb-run R CMD build $srcdir
pkg=$(ls *.gz)
xvfb-run R CMD check --as-cran $pkg
rm $pkg
echo "DONE"
