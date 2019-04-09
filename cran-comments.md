# Test environments

R CMD check and all tests passed on Linux, MacOS, and Windows, using the following environments:

* Debian GNU/Linux
  * R-devel r76333 (x86_64-pc-linux-gnu (64-bit))
  * R 3.5.3 (x86_64-pc-linux-gnu (64-bit))
* Mac OS X 10.12.6
  * R 3.5.2 (x86_64-apple-darwin15.6.0 (64-bit))
  * R 3.5.3 (x86_64-apple-darwin15.6.0 (64-bit))
* Via R-hub
  * Windows Server 2008 R2 SP1, R-release, 32/64 bit
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
* Via win_builder
  * R 3.6.0 alpha (2019-04-08 r76348) (x86_64-w64-mingw32 (64-bit))

# R CMD check results

* There were no NOTES, ERRORs or WARNINGs

# Downstream dependencies

None, since `revdepcheck::revdep_check()` reports no downnstream dependencies.
