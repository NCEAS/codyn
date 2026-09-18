# Test environments

R CMD check and all tests passed on Linux, MacOS, and Windows, using the following environments:

* Mac OS X 15.7.3
  * R 4.6.1 (aarch64-apple-darwin24.6.0 (64-bit))
* Via R-hub
  * Ubuntu 13.3.0-6ubuntu2~24.04.1 (R 4.6.1) (2020-11-27 r79522) x86_64-pc-linux-gnu
  * Ubuntu 24.04.5 LTS (R-devel (unstable) (2026-09-16 r90549)) x86_64-pc-linux-gnu
  * Windows Server 2022 x64 (build 26100) (R-devel (unstable) (2026-09-16 r90549 ucrt)) x86_64-w64-mingw32
  * macOS Sequoia 15.7.9 (R-devel (unstable) (2026-09-16 r90549)) x86_64-apple-darwin20
* Via win_builder
  * Windows Server 2022 x64 (build 20348) (R 4.6.1 (2026-06-24 ucrt)) x86_64-w64-mingw32
  * Windows Server 2022 x64 (build 20348) (R-devel (unstable) (2026-09-16 r90549 ucrt)) x86_64-w64-mingw32

# R CMD check results

* There was one NOTE, that words were potentially
misspelled in the DESCRIPTION. The flagged words are
proper names (Avolio, Hallett) or parts of citation text (the latin phrase "et al.", 
which is typically abbreviated as written).

# Downstream dependencies

None, since `devtools::revdep()` reports no downstream dependencies.
