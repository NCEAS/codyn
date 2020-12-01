# Test environments

R CMD check and all tests passed on Linux, MacOS, and Windows, using the following environments:

* Mac OS X 10.14.6
  * R 4.0.3 (x86_64-apple-darwin17.0 (64-bit))
* Ubuntu
  * R 4.0.2 (2020-06-22) (x86_64-pc-linux-gnu (64-bit))
* Via R-hub
  * R Under development (unstable) (2020-11-27 r79522) debian-gcc-devel
  * R Under development (unstable) (2020-10-24 r79367) fedora-clang-devel
  * R Under development (unstable) (2020-10-24 r79367) fedora-gcc-devel
  * R Under development (unstable) (2020-10-24 r79367) ubuntu-gcc-devel
* Via win_builder
  * R 4.0.3 (2020-10-10) (x86_64-w64-mingw32 (64-bit))
  * R 3.6.3 (2020-02-29) (x86_64-w64-mingw32 (64-bit))
  * R Under development (unstable) (2020-11-27 r79522) (x86_64-w64-mingw32 (64-bit))

# R CMD check results

* There was one NOTE, that words were potentially
misspelled in the DESCRIPTION. The flagged words are
proper names (Avolio, Hallett) or parts of citation text (the latin phrase "et al.", 
which is typically abbreviated as written).

# Downstream dependencies

None, since `devtools::revdep()` reports no downstream dependencies.
