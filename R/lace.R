#' @title Repeat rows
#'
#' @description `reprow` creates a `matrix` of a repeated vector where
#'   each row is the same.
#'
#' @name reprow
#'
#' @param x vector to be repeated
#'
#' @param n number of repeats
#'
#' @return `reprow` returns a `matrix` containing the repeated vector.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
reprow <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

#' @title Hamardman product
#'
#' @description `hamardman.prod` computes the Hamardman product of a vector
#'   of regression coefficients and a matrix of covariates.
#'
#' @param coef vector of regression coefficients
#'
#' @param covar a matrix of covariates
#'
#' @return `hamardman.prod` returns a `matrix` containg the Hamardman product
#'   of a vector of regression coefficients and a matrix of covariates.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
hamardman.prod <- function(coef, covar) {
  if (length(coef) != ncol(covar)) {
    stop("the number of coefficients is greater than the number of covariates")
  }
  output <- reprow(coef, nrow(covar)) * covar
  return(output)
}

#' @title IV-free exposure
#'
#' @description `iv_free` computes the IV-free exposure.
#'
#' @name iv_free
#'
#' @import matrixStats
#'
#' @param y vector of outcome values
#'
#' @param x vector of exposure values
#'
#' @param g the instrumental variable
#'
#' @param covar `data.frame` of covariates
#'
#' @param q the number of quantiles the exposure distribution is to be split
#'   into within which a causal effect will be fitted, known as localised
#'   average causal effects (LACE) (default: `10`)
#'
#' @param family a description of the error distribution and link function
#'   to be used in the model and is a `character` string naming either the
#'   gaussian (i.e. `"gaussian"` for continuous data) or binomial
#'   (i.e. `"binomial"` for binary data) family function (default: `"gaussian"`)
#'
#' @return `iv_free` returns a `list` of IV free exposure results:
#'
#' @return \item{xcoef}{the association between the exposure and the instrument}
#'
#' @return \item{x0}{the IV-free exposure}
#'
#' @return \item{x0q}{the quantiles of x0}
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
iv_free <- function(y, x, g, covar = NULL, q, family = "gaussian") {
  if (family == "gaussian") {
    if (!is.null(covar)) {
      model <- lm(x ~ g + covar)
    } else {
      model <- lm(x ~ g)
    }
    if (any(is.na(model$coef))) {
      stop(
        "there are missing regression coefficients in the regression of the ",
        "exposure on the instrument and covariates"
      )
    }
    x0 <- resid(model)
    xcoef <- model$coef[2]
  }
  if (family == "binomial") {
    if (!is.null(covar)) {
      model <- lm(x[y == 0] ~ g[y == 0] + covar[y == 0, , drop = FALSE])
      if (any(is.na(model$coef))) {
        stop(
          "there are missing regression coefficients in the regression of the ",
          "exposure on the instrument and covariates in the controls"
        )
      }
      x0 <- x - (model$coef[1] + model$coef[2] * g + rowSums(hamardman.prod(model$coef[3:length(model$coef)], covar)))
    } else {
      model <- lm(x[y == 0] ~ g[y == 0])
      if (any(is.na(model$coef))) {
        stop(
          "there are missing regression coefficients in the regression of the ",
          "exposure on the instrument and covariates in the controls"
        )
      }
      x0 <- x - (model$coef[1] + model$coef[2] * g)
    }
    xcoef <- model$coef[2]
  }
  quantiles <- quantile(x0, probs = seq(0, 1, 1 / q))
  x0q <- cut(x0, quantiles, include.lowest = TRUE, labels = FALSE)
  output <- list(xcoef = xcoef, x0 = x0, x0q = x0q)
  return(output)
}

#' @title Localised average causal effects
#'
#' @description `lace` computes the localised average causal effect
#'   for each quantile.
#'
#' @name lace
#'
#' @param y vector of outcome values
#'
#' @param x vector of exposure values
#'
#' @param g the instrumental variable
#'
#' @param covar `data.frame` of covariates
#'
#' @param q the number of quantiles the exposure distribution is to be split
#'   into within which a causal effect will be fitted, known as localised
#'   average causal effects (LACE) (default: `10`)
#'
#' @param x0q quantiles of x0 (the IV-free exposure)
#'
#' @param xc_sub compute the association between the exposure and the insturment
#'   in each quantile of x0
#'
#' @param family a description of the error distribution and link function
#'   to be used in the model and is a `character` string naming either the
#'   gaussian (i.e. `"gaussian"` for continuous data) or binomial
#'   (i.e. `"binomial"` for binary data) family function (default: `"gaussian"`)
#'
#' @param xpos the position used to relate `x` to the localised average causal
#'   effect, this can either be the mean of the x-values within each quantile
#'   or a percentile (e.g. 0.5 corresponds to the median value)
#'   (default: `"mean"`)
#'
#' @return `lace` returns a `data.frame` of the LACE for quantile:
#'
#' @return \item{coef}{the LACE in each quantile}
#'
#' @return \item{coef_se}{the standard error of the LACE in each quantile}
#'
#' @return \item{xmean}{the mean of the exposure in each quantile}
#'
#' @return \item{xcoef_sub}{
#'   the association between the exposure and the instrument in each quantile
#' }
#'
#' @return \item{xcoef_sub_se}{
#'   the standard error of the association between the exposure and the
#'   instrument in each quantile
#' }
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
lace <- function(y, x, g, covar = NULL,
  q, x0q, xc_sub = TRUE,
  family = "gaussian", xpos = "mean") {

  coef <- NULL
  coef_se <- NULL
  xmean <- NULL
  xcoef_sub <- NULL
  xcoef_sub_se <- NULL
  for (i in 1:q) {
    if (family == "gaussian") {
      if (is.null(covar)) {
        model <- lm(y[x0q == i] ~ g[x0q == i])
      } else {
        model <- lm(y[x0q == i] ~ g[x0q == i] + covar[x0q == i, , drop = FALSE])
      }
      if (is.na(model$coef[2])) {
        stop(
          "the regression coefficient of the outcome on the instrument in one ",
          "of the quantiles is missing"
        )
      }
      coef[i] <- model$coef[2]
      coef_se[i] <- summary(model)$coef[2, 2]
    }
    if (family == "binomial") {
      if (is.null(covar)) {
        model <- glm(y[x0q == i] ~ g[x0q == i], family = binomial)
      } else {
        model <- glm(
          y[x0q == i] ~ g[x0q == i] + covar[x0q == i, , drop = FALSE],
          family = binomial
        )
      }
      if (is.na(model$coef[2])) {
        stop(
          "the regression coefficient of the outcome on the instrument in one ",
          "of the quantiles is missing"
        )
      }
      coef[i] <- model$coef[2]
      coef_se[i] <- summary(model)$coef[2, 2]
    }
    if (xpos == "mean") {
      xmean[i] <- mean(x[x0q == i])
    }
    if (xpos != "mean") {
      xmean[i] <- quantile(x[x0q == i], probs = xpos)
    }
    if (xc_sub) {
      if (family == "gaussian") {
        if (is.null(covar)) {
          xcoefs <- lm(x[x0q == i] ~ g[x0q == i])
        } else {
          xcoefs <- lm(
            x[x0q == i] ~ g[x0q == i] + covar[x0q == i, , drop = FALSE]
          )
        }
        if (is.na(xcoefs$coef[2])) {
          stop(
            "the regression coefficient of the exposure on the instrument in ",
            "one of the quantiles is missing"
          )
        }
        xcoef_sub[i] <- xcoefs$coef[2]
        xcoef_sub_se[i] <- summary(xcoefs)$coef[2, 2]
      }
      if (family == "binomial") {
        if (is.null(covar)) {
          xcoefs <- lm(x[(x0q == i & y == 0)] ~ g[(x0q == i & y == 0)])
        } else {
          xcoefs <- lm(
            x[(x0q == i & y == 0)] ~
              g[(x0q == i & y == 0)] +
                covar[(x0q == i & y == 0), , drop = FALSE]
          )
        }
        if (is.na(xcoefs$coef[2])) {
          stop(
            "the regression coefficient of the exposure on the instrument in ",
            "one of the quantiles is missing"
          )
        }
        xcoef_sub[i] <- xcoefs$coef[2]
        xcoef_sub_se[i] <- summary(xcoefs)$coef[2, 2]
      }
    }
  }
  if (xc_sub == TRUE) {
    output <- data.frame(
      coef = coef, coef_se = coef_se, xmean = xmean,
      xcoef_sub = xcoef_sub, xcoef_sub_se = xcoef_sub_se
    )
  } else {
    output <- data.frame(coef = coef, coef_se = coef_se, xmean = xmean)
  }
  return(output)

}
