codyn - Community Dynamics Metrics
=====

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/codyn)](http://cran.r-project.org/package=codyn)
[![Build Status](https://travis-ci.org/laurenmh/codyn.png?branch=master)](https://travis-ci.org/laurenmh/codyn)



- **Authors**:Lauren Hallett [lauren.m.hallett@gmail.com], Sydney K. Jones [syd@sevilleta.unm.edu], Andrew A. MacDonald [aammacdonald@gmail.com], Matthew B. Jones [jones@nceas.ucsb.edu], Dan F. B. Flynn [flynn@fas.harvard.edu], Peter Slaughter [slaughter@nceas.ucsb.edu], Corinna Gries [cgries@wisc.edu], Scott L. Collins [scollins@sevilleta.unm.edu]
- **License**: [Apache 2](http://opensource.org/licenses/Apache-2.0)
- [Package source code on Github](https://github.com/laurenmh/codyn)
- [**Submit Bugs and feature requests**](https://github.com/laurenmh/codyn/issues)

A package to analyze long-term ecological community datasets.

The functions in `codyn` implement metrics that are explicitly temporal, and include the option to calculate them over multiple replicates. Functions fall into two categories: temporal diversity indices and community stability metrics. The diversity indices in `codyn` are temporal analogs to traditional diversity indices such as richness and rank-abundance curves. Specifically, `codyn` includes functions to calculate species turnover, mean rank shifts and lags in community similarity between time points. The community stability metrics in `codyn` calculate overall stability and patterns of species covariance and synchrony over time. Finally, `codyn` contains vignettes that describe methods and reproduce figures from published papers to help users contextualize and apply functions to their own data.

## Installation
From CRAN, the package can be installed using standard tools:
```R
install.packages("codyn")
```

Releases and pre-releases of the software are also available from the NCEAS drat repository, and
can be installed after drat has been installed using:
```R
drat::addRepo("NCEAS")
install.packages("codyn")
```

## Acknowledgements
Work on this package was supported by NSF-ABI grant #1262458 to C Gries, M Jones, and S Collins. Additional support
was provided for working group collaboration by the National Center for Ecological Analysis and Synthesis, a Center funded by the University of California, Santa Barbara, and the State of California.

