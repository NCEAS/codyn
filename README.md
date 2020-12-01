# codyn - Community Dynamics Metrics

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/codyn)](https://cran.r-project.org/package=codyn)
[![Build Status](https://travis-ci.com/NCEAS/codyn.png?branch=master)](https://travis-ci.com/NCEAS/codyn)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/codyn)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

- **Authors**: Lauren Hallett [lauren.m.hallett@gmail.com], Meghan Avolio [meghan.avolio@jhu.edu], Ian T. Carroll [carroll.ian@gmail.com], Sydney K. Jones [syd@sevilleta.unm.edu], Andrew A. MacDonald [a.a.m.macdonald@gmail.com],  Dan F. B. Flynn [flynn@fas.harvard.edu], Peter Slaughter [slaughter@nceas.ucsb.edu], Julie Ripplinger [julie.ripplinger@asu.edu], Scott L. Collins [scollins@sevilleta.unm.edu], Corinna Gries [cgries@wisc.edu], Matthew B. Jones [jones@nceas.ucsb.edu]
- Version 1.x: [doi:10.5063/F1542KJB](https://doi.org/10.5063/F1542KJB)
- Version 2.x: [doi:10.5063/F1N877Z6](https://doi.org/10.5063/F1N877Z6)
- **License**: [Apache 2](https://opensource.org/licenses/Apache-2.0)
- [Package source code on Github](https://github.com/NCEAS/codyn)
- [**Submit Bugs and feature requests**](https://github.com/NCEAS/codyn/issues)

A package to analyze long-term ecological community datasets.

Univariate and multivariate temporal and spatial diversity indices, 
rank abundance curves, and community stability metrics. The functions 
implement metrics that are either explicitly temporal and include the 
option to  calculate them over multiple replicates, or spatial and include 
the option to calculate them over multiple time points. Functions fall into 
five categories: static diversity indices, temporal diversity indices, 
spatial diversity indices, rank abundance curves, and community stability 
metrics. The diversity indices are temporal and spatial analogs to 
traditional diversity indices. Specifically, the package includes functions 
to calculate community richness, evenness and diversity at a given point in 
space and time. In addition, it contains functions to calculate species 
turnover, mean rank shifts, and lags in community similarity between two 
time points.

For an overview of __codyn__, see:
    
- Hallett LM, Jones SK, MacDonald AAA, Jones MB, Flynn DFB, Ripplinger J, Slaughter P, Gries C, Collins SL (2016) *codyn: An R package of community dynamics metrics.* Methods in Ecology and Evolution, 7(10):1146–1151. https://doi.org/10.1111/2041-210X.12569

For a description of the newer spatial methods in __codyn__ v2.x:

- Avolio ML, Carroll IT, Collins SL, Houseman GR, Hallett LM, Isbell F, Koerner SE, Komatsu KJ, Smith MD, Wilcox KR (2019) *A comprehensive approach to analyzing community dynamics using rank abundance curves.* Ecosphere, 10(10):e02881. https://doi.org/10.1002/ecs2.2881


## Installation
From CRAN, the package can be installed using standard tools:
```R
install.packages("codyn")
```

## Automated R CMD check with Docker via rhub

To simplify the process of running `R CMD check` on the package, one can easily run the build and tests using
the [`rhub`](https://github.com/r-hub/rhub) package. Use `rhub::platforms()` to get a list of platforms that can be used to build and test.

```r
library(rhub)
chks <- check(platform = c("debian-gcc-devel", "fedora-gcc-devel"), show_status = FALSE)
```

and the checks can be run locally using rhub as well using a docker container:

```r
library(rhub)
local_check_linux(image="rhub/fedora-gcc-devel")
```

## Acknowledgments

Work on this package was supported by NSF-ABI grant #1262458 to C. Gries, M. Jones, and S. Collins. Additional support was provided for working group collaboration by the National Center for Ecological Analysis and Synthesis, a Center funded by the University of California, Santa Barbara, and the State of California, and a SESYNC Synthesis Postdoctoral Fellowship to MLA.

