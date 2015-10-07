codyn - Community Dynamics Metrics
=====

[![Build Status](https://travis-ci.org/laurenmh/codyn.png?branch=master)](https://travis-ci.org/laurenmh/codyn)



- **Authors**:Lauren Hallett [lauren.m.hallett@gmail.com], Sydney Jones [syd@sevilleta.unm.edu], Andrew MacDonald [aammacdonald@gmail.com], Matthew Jones [jones@nceas.ucsb.edu], Dan Flynn [flynn@fas.harvard.edu], Peter Slaughter [slaughter@nceas.ucsb.edu], Corinna Gries [cgries@wisc.edu], Scott Collins [scollins@sevilleta.unm.edu]

- **License**: [Apache 2](http://opensource.org/licenses/Apache-2.0)
- [Package source code on Github](https://github.com/laurenmh/codyn)
- [**Submit Bugs and feature requests**](https://github.com/laurenmh/codyn/issues)

A package to analyze long-term ecological community datasets.

The functions in `codyn` implement metrics that are explicitly temporal, and include the option to calculate them over multiple replicates. Functions fall into two categories: temporal diversity indices and community stability metrics. The diversity indices in `codyn` are temporal analogs to traditional diversity indices such as richness and rank-abundance curves. Specifically, `codyn` includes functions to calculate species turnover, mean rank shifts and lags in community similarity between time points. The community stability metrics in `codyn` calculate overall stability and patterns of species covariance and synchrony over time. Finally, `codyn` contains vignettes that describe methods and reproduce figures from published papers to help users contextualize and apply functions to their own data.

## Installation
```R
install.packages("codyn")
```

Work on this package was supported by NSF-ABI grant #1262458

