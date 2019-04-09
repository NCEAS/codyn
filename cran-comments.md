# Test environments

R CMD check and all tests passed on Linux, MacOS, and Windows, using the following environments:

* Debian GNU/Linux
  * R-devel r76333 (x86_64-pc-linux-gnu (64-bit))
  * R 3.5.3 (x86_64-pc-linux-gnu (64-bit))
* Mac OS X 10.14.3
  * R 3.5.2 (x86_64-apple-darwin15.6.0 (64-bit))
* Windows (via win-builder)

# R CMD check results

* There were no NOTES, ERRORs or WARNINGs

# Downstream dependencies

Uncertain, since `devtools::revdep_check()` is now defunct.
