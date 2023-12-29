# nlmr <img src='man/figures/logo.png' align="right" height="139"/>

This package is used to assess non-linear exposure-outcome relationships using instrumental variable (IV) analysis in the context of Mendelian randomisation (MR). In this package, there are two IV methods for investigating the shape of the exposure-outcome relationship: a fractional polynomial method (`fracpoly_mr`) and a piecewise linear method (`piecewise_mr`). The population (i.e. one-sample) is divided into strata using the exposure distribution, and a causal effect is estimated, referred to as a localized average causal effect (LACE), in each stratum. The fractional polynomial method fits across these LACE using meta-regression. The piecewise linear method estimates a continuous piecewise linear function by consecutively adding the LACE together. 

## Functions
* `fracpoly_mr`: this method performs an MR analysis using fractional polynomials.  
* `piecewise_mr`: this method performs an MR analysis using a piecewise linear function.  

## Installation
```
install.packages("remotes")
remotes::install_github("jrs95/nlmr")
```

## Example
```
# Libraries
library(nlmr)

# IV (g), exposure (x) & outcome (y)
epsx <- rexp(10000)
u <- runif(10000, 0, 1)
g <- rbinom(10000, 2, 0.3)
epsy <- rnorm(10000)
ag <- 0.25
x <- 1 + ag * g + u + epsx
y <- 0.15 * x^2 + 0.8 * u + epsy

# Covariates (covar)
c1 <- rnorm(10000)
c2 <- rnorm(10000)
c3 <- rbinom(10000, 2, 0.33)
covar <- data.frame(c1 = c1, c2 = c2, c3 = as.factor(c3))

# Analyses
fp <- fracpoly_mr(
  y = y, x = x, g = g, covar = covar,
  family = "gaussian", q = 10, d = 1, ci = "model_se",
  fig = TRUE
)
summary(fp)
plm <- piecewise_mr(
  y = y, x = x, g = g, covar = covar,
  family = "gaussian", q = 10, nboot = 100,
  fig = TRUE
)
summary(plm)
```

## Citation 
Staley JR and Burgess S. Semiparametric methods for estimation of a non-linear exposure-outcome relationship using instrumental variables with application to Mendelian randomization. [Genet Epidemiol](https://pubmed.ncbi.nlm.nih.gov/28317167/) 2017;41(4):341-352.
