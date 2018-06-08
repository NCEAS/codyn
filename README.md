# codyn - Community Dynamics Metrics

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/codyn)](https://cran.r-project.org/package=codyn)
[![Build Status](https://travis-ci.org/NCEAS/codyn.png?branch=master)](https://travis-ci.org/NCEAS/codyn)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/codyn)


- **Authors**: Lauren Hallett [lauren.m.hallett@gmail.com], Meghan Avolio [meghan.avolio@jhu.edu], Ian T. Carroll [icarroll@sesync.org], Sydney K. Jones [syd@sevilleta.unm.edu], Andrew A. MacDonald [a.a.m.macdonald@gmail.com],  Dan F. B. Flynn [flynn@fas.harvard.edu], Peter Slaughter [slaughter@nceas.ucsb.edu], Julie Ripplinger [julie.ripplinger@asu.edu], Scott L. Collins [scollins@sevilleta.unm.edu], Corinna Gries [cgries@wisc.edu], Matthew B. Jones [jones@nceas.ucsb.edu]
- Version 1.x: [doi:10.5063/F1542KJB](http://doi.org/10.5063/F1542KJB)
- Version 2.x: [doi:10.5063/F1N877Z6](http://doi.org/10.5063/F1N877Z6)
- **License**: [Apache 2](http://opensource.org/licenses/Apache-2.0)
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
    
- Hallett et al. (2016) *codyn: An R package of community dynamics metrics*. Methods in Ecology and Evolution. http://doi.org/10.1111/2041-210X.12569

## Installation
From CRAN, the package can be installed using standard tools:
```R
install.packages("codyn")
```

## Automated R CMD check with Docker

To simplify the process of running `R CMD check` on the package, the source distribution on GitHub includes configuration
files to use [Docker](https://www.docker.com/) to download and build standard Debian-based images for the current release of 
R and the current development branch of R. Assuming you already have docker and docker-compose installed, these Docker 
configuration files allow a clean environment to be built and tested with a single command.  Checks can be run against the 
current stable release of R using:

```bash
$ docker-compose run --rm r-check-stable
```

and the checks can be run against the current unstable development version of R using:

```bash
$ docker-compose run --rm r-check-devel
```

## Acknowledgements
Work on this package was supported by NSF-ABI grant #1262458 to C. Gries, M. Jones, and S. Collins. Additional support was provided for working group collaboration by the National Center for Ecological Analysis and Synthesis, a Center funded by the University of California, Santa Barbara, and the State of California, and a SESYNC Synthesis Postdoctoral Fellowship to MLA.

