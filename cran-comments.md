# Test environments

* Debian GNU/Linux stretch/sid (via Docker), R 3.2.2, R 3.2.3, and R-devel (2016-01-12 r69936)
* Mac OS X 10.10.3, R 3.2.3
* Windows (via win-builder): x86_64-w64-mingw32 (64-bit), R 3.2.3, and R-devel

# R CMD check results

* There were no ERRORs, WARNINGs, or NOTEs on Debian and Mac OS X
* Under Windows, there was 1 NOTE, indicating possibly misspelled words in the DESCRIPTION;
all of these words have been checked and are spelled correctly, albeit technical science terms
    * covariance (20:75)
    * indices (16:5, 16:60, 17:38)
    * synchrony (21:9)

# Downstream dependencies

* There are currently no downstream dependencies for this package.
