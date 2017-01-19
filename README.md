# nl_mr
This package is used to assess non-linear exposure-outcome relationships using instrumental variable (IV) analysis in the context of Mendelian randomisation (MR). In this package, there are two IV methods for investigating the shape of the exposure-outcome relationship: a fractional polynomial method (frac_poly_mr) and a piecewise linear method (piecewise_mr). The population (i.e. one-sample) is divided into strata using the exposure distribution, and a causal effect is estimated, referred to as a localized average causal effect (LACE), in each stratum. The fractional polynomial method fits across these LACE using meta-regression. The piecewise linear method estimates a continuous piecewise linear function by consecutively adding the LACE together. 

# Functions
* frac_poly_mr - this method performs meta-regression across the LACE 
* piecewise_mr - this method

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("/jrs95/nl_mr")

# Example
