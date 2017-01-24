# nlmr
This package is used to assess non-linear exposure-outcome relationships using instrumental variable (IV) analysis in the context of Mendelian randomisation (MR). In this package, there are two IV methods for investigating the shape of the exposure-outcome relationship: a fractional polynomial method (frac_poly_mr) and a piecewise linear method (piecewise_mr). The population (i.e. one-sample) is divided into strata using the exposure distribution, and a causal effect is estimated, referred to as a localized average causal effect (LACE), in each stratum. The fractional polynomial method fits across these LACE using meta-regression. The piecewise linear method estimates a continuous piecewise linear function by consecutively adding the LACE together. 

# Functions
* frac_poly_mr - this method performs IV analysis using fractional polynomials 
* piecewise_mr - this method performs IV analysis using piecewise linear function

# Installation
1. install.packages("devtools")
2. library(devtools) 
3. install_github("jrs95/nlmr")
4. library(nlmr)

# Example
\#\#\# Instrumental variable (g), exposure (x) & outcome (y)  
epsx = rexp(10000)  
u    = runif(10000, 0, 1)  
g    = rbinom(10000, 2, 0.3)  
epsy = rnorm(10000)  
ag = 0.25  
x = 1 + ag\*g + u + epsx  
y = 0.1\*x^2 + 0.8\*u + epsy 

\#\#\# Covariates (c) & covariate types (c_type)  
c1 = rnorm(10000)  
c2 = rnorm(10000)  
c3 = rbinom(10000,2,0.33)  
c = data.frame(c1=c1, c2=c2, c3=as.factor(c3))  
c_type = c("numeric", "numeric", "factor")

\#\#\# Analyses  
fp = frac_poly_mr(y, x, g, c, c_type, family="gaussian", q=10, d=1, ci="model_se", fig=T)  
summary(fp)  
plm = piecewise_mr(y, x, g, c, c_type, family="gaussian", q=10, nboot=5, fig=T)  
summary(plm)

# Reference 
James R Staley & Stephen Burgess, Semiparametric methods for estimation of a non-linear exposure-outcome relationship using instrumental variables with application to Mendelian randomization. 
