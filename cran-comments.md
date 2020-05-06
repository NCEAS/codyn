# Test environments

R CMD check and all tests passed on Linux, MacOS, and Windows, using the following environments:

* Mac OS X 10.14.6
  * R 3.6.2 (x86_64-apple-darwin15.6.0 (64-bit))
  * R 4.0.0 (x86_64-apple-darwin17.0 (64-bit))
* Via R-hub
  * R 4.0.0 debian-gcc-release (r-release)
  * R Under development (unstable) (2020-05-01 r78341) (x86_64-pc-linux-gnu (64-bit))
* Via win_builder
  * R 4.0.0 (2020-04-24) (x86_64-w64-mingw32 (64-bit))
  * R Under development (unstable) (2020-05-05 r78369) (x86_64-w64-mingw32 (64-bit))

# R CMD check results

* There was one NOTE, that this is a new submission, and that words were potentially
misspelled in the DESCRIPTION.  This is not a new submission, as it is an 
update to fix tests that broke after the R 4.0.0 factor changes.  The flagged spelling word
is correctly spelled: 'indices'.

# Downstream dependencies

None, since `devtools::revdep()` reports no downnstream dependencies.
