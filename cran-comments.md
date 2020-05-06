# Test environments

R CMD check and all tests passed on Linux, MacOS, and Windows, using the following environments:

* Mac OS X 10.14.6
  * R 3.6.2 (x86_64-apple-darwin15.6.0 (64-bit))
* Via R-hub
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Windows Server 2008 R2 SP1, R-release, 32/64 bit
  * macOS 10.11 El Capitan, R-release (experimental)
  * Debian Linux, R-devel, GCC (2019-09-14 r77190) (x86_64-pc-linux-gnu (64-bit))
* Via win_builder
  * R 3.6.1 (2019-07-05) (x86_64-w64-mingw32 (64-bit))
  * R Under development (unstable) (2019-09-15 r77192) (x86_64-w64-mingw32 (64-bit))
  
# R CMD check results

* There were no NOTES, ERRORs or WARNINGs

# Downstream dependencies

None, since `revdepcheck::revdep_check()` reports no downnstream dependencies.
