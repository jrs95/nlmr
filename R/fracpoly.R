#' @title Fractional polynomial Mendelian randomization
#'
#' @description `fracpoly_mr` performs a Mendelian randomization (MR)
#'   analysis by fitting fractional polynomial models to localised
#'   average causal effects using meta-regression.
#'
#' @name fracpoly_mr
#'
#' @import metafor
#'
#' @param y vector of outcome values
#'
#' @param x vector of exposure values
#'
#' @param g the instrumental variable
#'
#' @param covar `data.frame` of covariates
#'
#' @param family a description of the error distribution and link function
#'   to be used in the model and is a `character` string naming either the
#'   gaussian (i.e. `"gaussian"` for continuous data) or binomial
#'   (i.e. `"binomial"` for binary data) family function (default: `"gaussian"`)
#'
#' @param q the number of quantiles the exposure distribution is to be split
#'   into within which a causal effect will be fitted, known as localised
#'   average causal effects (LACE) (default: `10`)
#'
#' @param xpos the position used to relate `x` to the localised average causal
#'   effect, this can either be the mean of the x-values within each quantile
#'   or a percentile (e.g. 0.5 corresponds to the median value)
#'   (default: `"mean"`)
#'
#' @param method meta-regression method parsed to the `metafor` package,
#'   e.g. the fixed-effects model (default: `"FE"`)
#'
#' @param d fractional polynomial degree, the options are: `1`, `2`,
#'   or `"both"` (default: `1`)
#'
#' @param pd p-value cut-off for choosing the best-fitting fractional polynomial
#'   of degree 2 over the best-fitting fractional polynomial degree 1, used only
#'   if `d = "both"` (default: `0.05`)
#'
#' @param ci the type of 95\% confidence interval, there are three options:
#'   (i) using the model standard errors (`"model_se"`),
#'   (ii) using bootstrap standard errors (`"bootstrap_se"`), and
#'   (iii) using bootstrap percentile confidence intervals (`"bootstrap_per"`)
#'   (default: `"model_se"`)
#'
#' @param nboot the number of bootstrap replications (default: `100`)
#'
#' @param fig a `logical` statement as to whether the user wants the results
#'   displayed in a figure (default: `FALSE`)
#'
#' @param ref the reference point for the figure, this is the value of the
#'   function that represents the expected difference in the outcome compared
#'   with this reference value when the exposure is set to different values
#'   (default: `mean(x)`)
#'
#' @param pref_x the prefix/label for the x-axis (default: `"x"`)
#'
#' @param pref_x_ref the prefix for the reference value displayed on the y-axis
#'   (default: `"x"`)
#'
#' @param pref_y the prefix/label for the y-axis (default: `"y"`)
#'
#' @param ci_type the type of confidence interval to be displayed on the
#'   graph, where confidence intervals are either presented as bands across the
#'   range of x (option: `"overall"`) or as error bars at the mean in each
#'   quantile of x (option: `"quantile"`) (default:)
#'
#' @param ci_quantiles the number of quantiles at which confidence intervals
#'   are to be displayed (default: `10`)
#'
#' @param breaks breaks on the y-axis of the figure
#'
#' @return `fracpoly_mr` returns a `list` of non-linear MR results from the
#'   fractional polynomial MR approach:
#'
#' @return \item{n}{number of individuals}
#'
#' @return \item{model}{
#'   the model specifications:
#'     number of quantiles (`q`),
#'     position used to relate x to the LACE in each quantiles (`xpos`),
#'     type of confidence interval constructed (`ci`),
#'     number of bootstrap replications performed (`nboot`)
#' }
#'
#' @return \item{powers}{the powers of the chosen polynomial}
#'
#' @return \item{coefficients}{
#'   the regression estimates:
#'     regression coefficients (`beta`),
#'     standard errors of regression coefficients (`se`),
#'     lower 95% confidence interval (`lci`),
#'     upper 95% confidence interval (`uci`),
#'     p-value (`pval`)
#' }
#'
#' @return \item{lace}{
#'   the localised average causal effect estimate in each quantile:
#'     regression coefficients (`beta`),
#'     standard errors of regression coefficients (`se`),
#'     lower 95% confidence interval (`lci`),
#'     upper 95% confidence interval (`uci`),
#'     p-value (`pval`)
#' }
#'
#' @return \item{xcoef}{
#'   the association between the instrument and the exposure in each quantile:
#'     regression coefficients (`beta`),
#'     standard errors of regression coefficients (`se`)
#' }
#'
#' @return \item{p_tests}{
#'   the p-value of the non-linearity tests:
#'     p-value of the test between the fractional polynomial degrees (`fp_d1_d2`),
#'     p-value from the fractional polynomial non-linearity test (`fp`),
#'     p-value from the quadratic test (`quad`),
#'     p-value from the Cochran Q test (`Q`)
#' }
#'
#' @return \item{p_heterogeneity}{
#'   the p-value of heterogeneity:
#'     p-value of the Cochran Q heterogeneity test (`Q`),
#'     p-value from the trend test (`trend`).
#' }
#'
#' @examples
#' # IV (g), exposure (x) & outcome (y)
#' epsx <- rexp(10000)
#' u <- runif(10000, 0, 1)
#' g <- rbinom(10000, 2, 0.3)
#' epsy <- rnorm(10000)
#' ag <- 0.25
#' x <- 1 + ag * g + u + epsx
#' y <- 0.15 * x^2 + 0.8 * u + epsy
#'
#' # Covariates (covar)
#' c1 <- rnorm(10000)
#' c2 <- rnorm(10000)
#' c3 <- rbinom(10000, 2, 0.33)
#' covar <- data.frame(c1 = c1, c2 = c2, c3 = as.factor(c3))
#'
#' # Analyses
#' fp <- fracpoly_mr(
#'   y = y, x = x, g = g, covar = covar,
#'   family = "gaussian", q = 10, d = 1, ci = "model_se",
#'   fig = TRUE
#' )
#' summary(fp)
#' plm <- piecewise_mr(
#'   y = y, x = x, g = g, covar = covar,
#'   family = "gaussian", q = 10, nboot = 100,
#'   fig = TRUE
#' )
#' summary(plm)
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
fracpoly_mr <- function(y, x, g, covar = NULL,
  family = "gaussian", q = 10, xpos = "mean",
  method = "FE", d = 1, pd = 0.05, ci = "model_se", nboot = 100,
  fig = FALSE, ref = mean(x), pref_x = "x", pref_x_ref = "x", pref_y = "y",
  ci_type = "overall", ci_quantiles = 10, breaks = NULL) {

  # Errors
  if (!(is.vector(y) && is.vector(x) && is.vector(g))) {
    stop("either the outcome, exposure or instrument is not a vector")
  }
  if (!is.null(covar) && !is.data.frame(covar)) {
    stop("covar has to be a data.frame")
  }
  if (
    !(
      (is.numeric(y) || is.integer(y)) &&
        (is.numeric(x) || is.integer(x)) &&
        (is.numeric(g) || is.integer(g))
    )
  ) {
    stop("either the outcome, exposure or instrument is not numeric")
  }
  if (any(x <= 1)) {
    stop("fractional polynomial models require the exposure to be >>1")
  }
  if (length(y) <= 1) {
    stop("the outcome is less than or equal to a single value")
  }
  if (length(y) != length(x) || length(y) != length(g)) {
    stop(
      "the number of observations for the outcome, exposure and instrument ",
      "are not all the same"
    )
  }
  if (!is.null(covar) && nrow(covar) != length(y)) {
    stop(
      "the number of observations for the outcome and covariates are not ",
      "the same"
    )
  }
  if (any(is.na(y)) || any(is.na(x)) || any(is.na(g))) {
    stop(
      "there are missing values in either the outcome, exposure or ",
      "instrument"
    )
  }
  if (!is.null(covar) && any(is.na(covar))) {
    stop("there are missing values in the covariates")
  }
  if (!(family == "gaussian" || family == "binomial")) {
    stop("family has to be equal to either \"gaussian\" or \"binomial\"")
  }
  if (family == "binomial" && any(!(y == 1 | y == 0))) {
    stop("y has to be 0 or 1 if family is equal to \"binomial\"")
  }
  if ((length(y) / 10) < q) {
    stop("the quantiles should contain at least 10 observations")
  }
  if (!(xpos == "mean" || (xpos > 0 && xpos < 1))) {
    stop("the position used to relate x to the localised average causal effect")
  }
  if (!(d == 1 || d == 2 || d == "both")) {
    stop("the degree has to be equal to 1, 2 or \"both\"")
  }
  if (!(ci == "model_se" || ci == "bootstrap_se" || ci == "bootstrap_per")) {
    stop(
      "the confidence intervals must be one of \"model_se\", \"bootstrap_se\" ",
      "or \"bootstrap_per\""
    )
  }

  # Covariates
  if (!is.null(covar)) {
    covar <- model.matrix(as.formula(~.), data = covar)[, -1, drop = FALSE]
    if (any(is.na(covar))) {
      stop("there are missing values in the covariates")
    }
  }

  # x0 (IV-Free)
  ivf <- iv_free(y = y, x = x, g = g, covar = covar, q = q, family = family)
  x0 <- ivf$x0
  xcoef <- ivf$xcoef
  x0q <- ivf$x0q

  # LACE
  loc <- lace(
    y = y, x = x, g = g, covar = covar,
    q = q, x0q = x0q, xc_sub = TRUE,
    family = family, xpos = xpos
  )
  coef <- loc$coef / xcoef
  coef_se <- loc$coef_se / xcoef
  xmean <- loc$xmean
  xcoef_sub <- loc$xcoef_sub
  xcoef_sub_se <- loc$xcoef_sub_se

  # Test of IV-exposure assumption
  p_het <- 1 - pchisq(rma(xcoef_sub, vi = (xcoef_sub_se)^2)$QE, df = (q - 1))
  p_het_trend <- rma.uni(
    xcoef_sub ~ xmean,
    vi = xcoef_sub_se^2,
    method = method
  )$pval[2]

  # Best-fitting fractional polynomial
  fracpb <- fracpoly_best(
    coef = coef, coef_se = coef_se, xmean = xmean,
    d = d, pd = pd, method = method
  )
  model <- fracpb$model
  p_ML <- fracpb$p_ML
  p1_ML <- fracpb$p1_ML
  p2_ML <- fracpb$p2_ML
  fp_p <- fracpb$fp_p
  fp_d12_p <- fracpb$fp_d12_p
  d <- fracpb$d

  # Non-linearity tests
  p_quadratic <- rma(coef ~ xmean, (coef_se)^2, method = "FE")$pval[2]
  p_Q <- 1 - pchisq(rma(coef, vi = (coef_se)^2)$QE, df = (q - 1))

  # Bootstrap
  if (ci == "bootstrap_per" || ci == "bootstrap_se") {
    frac_coef_boot <- fracpoly_boot(
      y = y, x = x, g = g, covar = covar,
      q = q, x0q = x0q, xcoef = xcoef,
      family = family, xpos = xpos, method = method, nboot = nboot,
      d = d, p_ML = p_ML, p1_ML = p1_ML, p2_ML = p2_ML
    )
  } else {
    frac_coef_boot <- NULL
  }

  # Output
  if (d == 1) {
    powers <- p_ML + 1
  } else {
    powers <- c((p1_ML + 1), (p2_ML + 1))
  }
  beta <- as.numeric(model$b)
  if (ci == "model_se") {
    cov <- model$vb
    se <- model$se
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * pnorm(-abs(beta / se))
  }
  if (ci == "bootstrap_se") {
    cov <- var(frac_coef_boot)
    se <- sqrt(diag(cov))
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * pnorm(-abs(beta / se))
  }
  if (ci == "bootstrap_per") {
    if (d == 1) {
      se <- NA
      lci <- quantile(frac_coef_boot, probs = 0.025)
      uci <- quantile(frac_coef_boot, probs = 0.975)
      pval <- NA
    }
    if (d == 2) {
      se <- rep(NA, 2)
      lci <- NULL
      uci <- NULL
      pval <- NULL
      lci[1] <- quantile(frac_coef_boot[, 1], probs = 0.025)
      lci[2] <- quantile(frac_coef_boot[, 2], probs = 0.025)
      uci[1] <- quantile(frac_coef_boot[, 1], probs = 0.975)
      uci[2] <- quantile(frac_coef_boot[, 2], probs = 0.975)
      pval <- rep(NA, 2)
    }
  }
  lci <- as.numeric(lci)
  uci <- as.numeric(uci)
  if (ci == "model_se") {
    nboot <- NA
  }

  # Figure
  if (fig == TRUE) {
    if (d == 1) {
      figure <- fracpoly_figure(
        beta = beta, cov = cov, x.min = min(x), x.max = max(x),
        family = family, d = d, p_ML = p_ML, ci = ci,
        frac_coef_boot = frac_coef_boot, ref = ref,
        pref_x = pref_x, pref_x_ref = pref_x_ref, pref_y = pref_y,
        ci_type = ci_type, ci_quantile = ci_quantile, breaks = breaks
      )
    }
    if (d == 2) {
      figure <- fracpoly_figure(
        beta = beta, cov = cov, x.min = min(x), x.max = max(x),
        family = family, d = d, p1_ML = p1_ML, p2_ML = p2_ML, ci = ci,
        frac_coef_boot = frac_coef_boot, ref = ref,
        pref_x = pref_x, pref_x_ref = pref_x_ref, pref_y = pref_y,
        ci_type = ci_type, ci_quantile = ci_quantile, breaks = breaks
      )
    }
  }

  # Output
  model <- as.matrix(
    data.frame(q = q, xpos = xpos, ci_type = ci, nboot = nboot)
  )
  coefficients <- as.matrix(
    data.frame(beta = beta, se = se, lci = lci, uci = uci, pval = pval)
  )
  rownames(coefficients) <- powers
  if (nrow(coefficients) == 2) {
    if (powers[1] == powers[2]) {
      rownames(coefficients) <- c(powers[1], paste0("log ", powers[2]))
    }
  }
  loc <- as.matrix(
    data.frame(
      beta = coef, se = abs(coef_se),
      lci = coef - 1.96 * abs(coef_se),
      uci = coef + 1.96 * abs(coef_se),
      pval = 2 * pnorm(-abs(coef / coef_se))
    )
  )
  rownames(loc) <- seq_len(nrow(loc))
  xcoef_quant <- as.matrix(data.frame(beta = xcoef_sub, se = xcoef_sub_se))
  rownames(xcoef_quant) <- seq_len(nrow(xcoef_quant))
  p_tests <- as.matrix(
    data.frame(fp_d1_d2 = fp_d12_p, fp = fp_p, quad = p_quadratic, Q = p_Q)
  )
  p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
  if (fig == FALSE) {
    output <- list(
      n = length(y), model = model, powers = powers,
      coefficients = coefficients, lace = loc, xcoef = xcoef_quant,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity
    )
  } else {
    output <- list(
      n = length(y), model = model, powers = powers,
      coefficients = coefficients, lace = loc, xcoef = xcoef_quant,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity,
      figure = figure
    )
  }
  class(output) <- "fracpoly_mr"

  # Return
  return(output)

}

#' @title Print fracpoly_mr
#'
#' @description print method for class `"fracpoly_mr"`.
#'
#' @param x an object of class `"fracpoly_mr"`
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
print.fracpoly_mr <- function(x, ...) {
  cat("\nCall: \nfracpoly_mr")
  cat("\n\nPowers:\n")
  cat(x$powers)
  cat("\n\nCoefficients:\n")
  cat(x$coefficients[, 1])
  cat("\n\n")
  if (!is.null(x$figure)) {
    plot(x$figure)
  }
}

#' @title Summary of fracpoly_mr
#'
#' @description summary method for class `"fracpoly_mr"`.
#'
#' @param x an object of class `"fracpoly_mr"`
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
summary.fracpoly_mr <- function(x, ...) {
  model <- as.data.frame(x$model)
  powers <- x$powers
  n <- x$n
  coefficients <- as.data.frame(x$coefficients)
  if (model$ci_type == "bootstrap_per") {
    coefficients <- coefficients[, c(1, 3, 4)]
  }
  p_tests <- as.data.frame(x$p_tests)
  p_heterogeneity <- as.data.frame(x$p_heterogeneity)
  if (is.null(x$figure)) {
    output <- list(
      model = model, powers = powers, n = n,
      coefficients = coefficients,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity
    )
  }
  if (!is.null(x$figure)) {
    output <- list(
      model = model, powers = powers, n = n,
      coefficients = coefficients,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity,
      figure = x$figure
    )
  }
  class(output) <- "summary.fracpoly_mr"
  return(output)
}

#' @title Print summary of fracpoly_mr
#'
#' @description print.summary method for class `"fracpoly_mr"`.
#'
#' @param x an object of class `"fracpoly_mr"`
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
print.summary.fracpoly_mr <- function(x, ...) {
  cat("Call: fracpoly_mr")
  ci_type <- "Model based SEs"
  if (x$model$ci_type == "bootstrap_se") {
    ci_type <- "Bootstrap based SEs"
  }
  if (x$model$ci_type == "bootstrap_per") {
    ci_type <- "Percentile bootstrap"
  }
  if (ci_type == "Model based SEs") {
    cat(
      "\n\nNumber of individuals: ", x$n, "; ",
      "Quantiles: ", as.character(x$model$q), "; ",
      "95%CI: ", ci_type,
      sep = ""
    )
  }
  if (ci_type != "Model based SEs") {
    cat(
      "\n\nNumber of individuals: ", x$n, "; ",
      "Quantiles: ", as.character(x$model$q), "; ",
      "95%CI: ", ci_type, "; ",
      "Number of bootstrap replications: ", as.character(x$model$nboot),
      sep = ""
    )
  }
  cat("\n\nPowers:", x$powers)
  cat("\n\nCoefficients:\n")
  if (ci_type == "Percentile bootstrap") {
    names(x$coefficients) <- c("Estimate", "95%CI Lower", "95%CI Upper")
    printCoefmat(x$coefficients)
  }
  if (ci_type != "Percentile bootstrap") {
    names(x$coefficients) <- c(
      "Estimate", "Std. Error", "95%CI Lower", "95%CI Upper", "p.value"
    )
    printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
  }
  cat("\nNon-linearity tests")
  cat(
    "\nFractional polynomial degree p-value:",
    signif(x$p_tests$fp_d1_d2, digits = 3)
  )
  cat(
    "\nFractional polynomial non-linearity p-value:",
    signif(x$p_tests$fp, digits = 3)
  )
  cat("\nQuadratic p-value:", signif(x$p_tests$quad, digits = 3))
  cat("\nCochran Q p-value:", signif(x$p_tests$Q, digits = 3))
  cat("\n\nHeterogeneity tests")
  cat("\nCochran Q p-value:", signif(x$p_heterogeneity$Q, digits = 3))
  cat("\nTrend p-value:", signif(x$p_heterogeneity$trend, digits = 3))
  cat("\n")
  if (!is.null(x$figure)) {
    plot(x$figure)
  }
}

#' @title Best-fitting fractional polynomial
#'
#' @description `fracpoly_best` computes the best-fitting fractional polynomial
#'   of degrees 1 and 2.
#'
#' @import metafor
#'
#' @param coef coefficients of the localised average causal effects
#'
#' @param coef_se standard errors of the localised average causal effects
#'
#' @param xmean the mean of the exposure in each quantile
#'
#' @param d fractional polynomial degree, the options are: `1`, `2`,
#'   or `"both"` (default: `1`)
#'
#' @param pd p-value cut-off for choosing the best-fitting fractional polynomial
#'   of degree 2 over the best-fitting fractional polynomial degree 1, used only
#'   if `d = "both"` (default: `0.05`)
#'
#' @param method meta-regression method parsed to the `metafor` package,
#'   e.g. the fixed-effects model (default: `"FE"`)
#'
#' @return List of non-linear MR results from the fractional polynomial MR
#'   approach:
#'
#' @return \item{model}{
#'   the best-fitting fractional polynomial of either 1 or 2
#' }
#'
#' @return \item{p_ML}{
#'   the power of the best-fitting fractional polynomial of degree 1
#' }
#'
#' @return \item{p1_ML}{
#'   the first power of the best-fitting fractional polynomial of degree 2
#' }
#'
#' @return \item{p2_ML}{
#'   the second power of the best-fitting fractional polynomial of degree 2
#' }
#'
#' @return \item{fp_p}{
#'   the p-value testing the best-fitting fractional polynomial of degree 1
#'   against the linear model
#' }
#'
#' @return \item{fp_d12_p}{
#'   the p-value testing the best-fitting fractional polynomial of degree 2
#'   against the best-fitting fractional polnomial of degree 1
#' }
#'
#' @return \item{d}{the degree of the best fitting fractional polynomial}
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fracpoly_best <- function(coef, coef_se, xmean,
  d = 1, pd = 0.05, method = "FE") {

  # FP degree 1
  powers <- c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
  likelihood_d1 <- NULL

  for (p1 in powers) {
    if (p1 == -1) {
      x1 <- xmean^p1
    } else {
      x1 <- (p1 + 1) * xmean^p1
    }
    fp_mod <- try(
      rma(coef ~ -1 + x1, vi = (coef_se)^2, method = method),
      silent = TRUE
    )
    if (is(fp_mod, "try-error") == TRUE) {
      likelihood_d1 <- c(likelihood_d1, NA)
    } else {
      if (p1 == 0) {
        fp1 <- fp_mod
        p_ML <- p1
      } else {
        if (
          fp_mod$fit.stats[1, 1] >=
            suppressWarnings(max(likelihood_d1, na.rm = TRUE))
        ) {
          fp1 <- fp_mod
          p_ML <- p1
        }
      }
      likelihood_d1 <- c(likelihood_d1, fp_mod$fit.stats[1, 1])
    }
  }

  maxlik_d1 <- max(likelihood_d1, na.rm = TRUE)
  fp_p <- 1 - pchisq(((-2 * likelihood_d1[1]) - (-2 * maxlik_d1)), df = 1)

  # FP degree 2
  powers1 <- c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
  powers2 <- c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
  likelihood_d2 <- NULL

  for (p11 in powers1) {
    if (p11 == -1) {
      x1 <- xmean^p11
    } else {
      x1 <- (p11 + 1) * xmean^p11
    }
    for (p21 in powers2) {
      if (p11 == p21) {
        if (p21 == -1) {
          x2 <- 2 * (xmean^p21) * log(xmean)
        } else {
          x2 <- ((p21 + 1) * (xmean^p21) * log(xmean) + xmean^p21)
        }
      } else {
        if (p21 == -1) {
          x2 <- xmean^p21
        } else {
          x2 <- (p21 + 1) * xmean^p21
        }
      }
      fp_mod <- try(
        rma(coef ~ -1 + x1 + x2, vi = (coef_se)^2, method = method),
        silent = TRUE
      )
      if (is(fp_mod, "try-error") == TRUE) {
        likelihood_d2 <- c(likelihood_d2, NA)
      } else {
        if (p11 == 0 && p21 == 0) {
          fp2 <- fp_mod
          p1_ML <- p11
          p2_ML <- p21
        } else {
          if (
            fp_mod$fit.stats[1, 1] >=
              suppressWarnings(max(likelihood_d2, na.rm = TRUE))
          ) {
            fp2 <- fp_mod
            p1_ML <- p11
            p2_ML <- p21
          }
        }
        likelihood_d2 <- c(likelihood_d2, fp_mod$fit.stats[1, 1])
      }
    }
    powers2 <- powers2[-1]
  }

  maxlik_d2 <- max(likelihood_d2, na.rm = TRUE)
  fp_d12_p <- 1 - pchisq(((-2 * maxlik_d1) - (-2 * maxlik_d2)), df = 2)

  # Model
  if (d == "both") {
    if (fp_d12_p > pd) {
      d <- 1
    } else {
      d <- 2
    }
  }
  if (d == 1) {
    model <- fp1
    if (length(model$b) != 1) {
      stop(
        "incorrect number of parameters for best fitting fractional ",
        "polynomial of degree 1"
      )
    }
  }
  if (d == 2) {
    model <- fp2
    if (length(model$b) != 2) {
      stop(
        "incorrect number of parameters for best fitting fractional ",
        "polynomial of degree 2"
      )
    }
  }

  # Ouptut
  output <- list(
    model = model, p_ML = p_ML, p1_ML = p1_ML, p2_ML = p2_ML,
    fp_p = fp_p, fp_d12_p = fp_d12_p, d = d
  )

  # Return
  return(output)

}

#' @title Fractional polynomial bootstrap
#'
#' @description fracpoly_boot computes the best-fitting fractional polynomial of
#'   degrees 1 and 2.
#'
#' @name fracpoly_boot
#'
#' @import metafor
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
#' @param xcoef the association between the exposure and the instrument
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
#' @param method meta-regression method parsed to the `metafor` package,
#'   e.g. the fixed-effects model (default: `"FE"`)
#'
#' @param nboot the number of bootstrap replications (default: `100`)
#'
#' @param d fractional polynomial degree, the options are: `1`, `2`,
#'   or `"both"` (default: `1`)
#'
#' @param p_ML the power of the best-fitting fractional polynomial of degree 1
#'   (default: `NULL`)
#'
#' @param p1_ML the first power of the best-fitting fractional polynomial of
#'   degree 2 (default: `NULL`)
#'
#' @param p2_ML the second power of the best-fitting fractional polynomial of
#'   degree 2 (default: `NULL`)
#'
#' @return \item{frac_coef_boot}{
#'   a matrix of fractional polynomial results for each bootstrap replicates
#' }
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fracpoly_boot <- function(y, x, g, covar, q, x0q, xcoef,
  family = "gaussian", xpos = "mean", method = "FE",
  nboot, d, p_ML = NULL, p1_ML = NULL, p2_ML = NULL) {

  # Set-up
  frac_coef_boot <- matrix(0, nrow = nboot, ncol = d)

  # Bootstrapping
  for (i in 1:nboot) {

    ## Bootstrap sample
    indices <- sample.int(length(y), size = length(y), replace = TRUE)
    y_boot <- y[indices]
    x_boot <- x[indices]
    g_boot <- g[indices]
    if (!is.null(covar)) {
      covar_boot <- as.matrix(covar[indices, , drop = FALSE])
    } else {
      covar_boot <- NULL
    }
    x0qb <- x0q[indices]

    ## LACE
    loc <- lace(
      y = y_boot, x = x_boot, g = g_boot, covar = covar_boot,
      q = q, x0q = x0qb, xc_sub = FALSE,
      family = family, xpos = xpos
    )
    coef_boot <- loc$coef / xcoef
    coef_se_boot <- loc$coef_se / xcoef
    xmean_boot <- loc$xmean
    xcoef_sub_boot <- loc$xcoef_sub
    xcoef_sub_se_boot <- loc$xcoef_sub_se

    ## Bootstrap fractional polynomials
    if (d == 1) {
      if (p_ML == -1) {
        x1_boot <- xmean_boot^p_ML
      } else {
        x1_boot <- (p_ML + 1) * xmean_boot^p_ML
      }
      mod <- rma.uni(
        coef_boot ~ -1 + x1_boot,
        vi = (coef_se_boot)^2,
        method = method
      )
      if (length(mod$b) != 1) {
        stop(
          "incorrect number of parameters for best fitting fractional ",
          "polynomial of degree 1 in bootstrap sample"
        )
      }
      frac_coef_boot[i, 1] <- mod$b[1]
    } else {
      if (p1_ML == -1) {
        x1_boot <- xmean_boot^p1_ML
      } else {
        x1_boot <- (p1_ML + 1) * xmean_boot^p1_ML
      }
      if (p1_ML == p2_ML) {
        if (p2_ML == -1) {
          x2_boot <- 2 * (xmean_boot^p2_ML) * log(xmean_boot)
        } else {
          x2_boot <- (p2_ML + 1) * (xmean_boot^p2_ML) * log(xmean_boot) +
            xmean_boot^p2_ML
        }
      } else {
        if (p2_ML == -1) {
          x2_boot <- xmean_boot^p2_ML
        } else {
          x2_boot <- (p2_ML + 1) * xmean_boot^p2_ML
        }
      }
      mod <- rma.uni(
        coef_boot ~ -1 + x1_boot + x2_boot,
        vi = (coef_se_boot)^2,
        method = method
      )
      if (length(mod$b) != 2) {
        stop(
          "incorrect number of parameters for best fitting fractional ",
          "polynomial of degree 2 in bootstrap sample"
        )
      }
      frac_coef_boot[i, 1] <- mod$b[1]
      frac_coef_boot[i, 2] <- mod$b[2]
    }

  }

  # Return
  return(frac_coef_boot)

}

#' @title Fractional polynomial figure
#'
#' @description `fracpoly_figure` plots the best-fitting fractional polynomial.
#'
#' @name fracpoly_figure
#'
#' @import ggplot2
#'
#' @param beta the coefficients of the best-fitting fractional polynomial
#'
#' @param cov the covariance matrix of the best-fitting fractional polynomial
#'
#' @param x.min the minimum value of the exposure
#'
#' @param x.max the maximum value of the exposure
#'
#' @param family a description of the error distribution and link function
#'   to be used in the model and is a `character` string naming either the
#'   gaussian (i.e. `"gaussian"` for continuous data) or binomial
#'   (i.e. `"binomial"` for binary data) family function (default: `"gaussian"`)
#'
#' @param d fractional polynomial degree, the options are: `1`, `2`,
#'   or `"both"` (default: `1`)
#'
#' @param p_ML the power of the best-fitting fractional polynomial of degree 1
#'   (default: `NULL`)
#'
#' @param p1_ML the first power of the best-fitting fractional polynomial of
#'   degree 2 (default: `NULL`)
#'
#' @param p2_ML the second power of the best-fitting fractional polynomial of
#'   degree 2 (default: `NULL`)
#'
#' @param ci the type of 95\% confidence interval, there are three options:
#'   (i) using the model standard errors (`"model_se"`),
#'   (ii) using bootstrap standard errors (`"bootstrap_se"`), and
#'   (iii) using bootstrap percentile confidence intervals (`"bootstrap_per"`)
#'   (default: `"model_se"`)
#'
#' @param frac_coef_boot a matrix of fractional polynomial results for each
#'   bootstrap replicates (default: `NULL`)
#'
#' @param ref the reference point for the figure, this is the value of the
#'   function that represents the expected difference in the outcome compared
#'   with this reference value when the exposure is set to different values
#'
#' @param pref_x the prefix/label for the x-axis (default: `"x"`)
#'
#' @param pref_x_ref the prefix for the reference value displayed on the y-axis
#'   (default: `"x"`)
#'
#' @param pref_y the prefix/label for the y-axis (default: `"y"`)
#'
#' @param ci_type the type of confidence interval to be displayed on the
#'   graph, where confidence intervals are either presented as bands across the
#'   range of x (option: `"overall"`) or as error bars at the mean in each
#'   quantile of x (option: `"quantile"`) (default:)
#'
#' @param ci_quantiles the number of quantiles at which confidence intervals
#'   are to be displayed (default: `10`)
#'
#' @param breaks breaks on the y-axis of the figure
#'
#' @return `fracpoly_figure` returns a plot of the best-fitting fractional
#'   polynomial.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fracpoly_figure <- function(beta, cov, x.min, x.max, family = "gaussian",
  d = 1, p_ML = NULL, p1_ML = NULL, p2_ML = NULL,
  ci = "model_se", frac_coef_boot = NULL, ref,
  pref_x = "x", pref_x_ref = "x", pref_y = "y",
  ci_type = "overall", ci_quantiles = 10, breaks = NULL) {

  if (ci_type == "overall") {
    plot.data <- data.frame(x = runif(10000, x.min, x.max))
    plot.data.ref <- data.frame(x = ref, y = 0)
    if (d == 1) {
      if (p_ML == -1) {
        plot.data$yest <- beta * log(plot.data$x) - (beta * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((log(plot.data$x) - log(ref))^2 * cov[1, 1])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot) -
            reprow(log(ref) %*% t(frac_coef_boot), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p_ML != -1) {
        plot.data$yest <- beta * plot.data$x^(p_ML + 1) - beta * ref^(p_ML + 1)
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p_ML + 1) - ref^(p_ML + 1))^2 * cov[1, 1]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p_ML + 1) %*% t(frac_coef_boot[, 1]) -
            reprow(
              ref^(p_ML + 1) %*% t(frac_coef_boot[, 1]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
    }
    if (d == 2) {
      if (p1_ML == -1 && p2_ML == -1) {
        plot.data$yest <- beta[1] * log(plot.data$x) +
          beta[2] * log(plot.data$x) * log(plot.data$x) -
          (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (log(plot.data$x) - log(ref))^2 * cov[1, 1] +
              2 * (log(plot.data$x) - log(ref)) * (log(plot.data$x) * log(plot.data$x) - log(ref) * log(ref)) * cov[1, 2] +
              (log(plot.data$x) * log(plot.data$x) - log(ref) * log(ref))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot[, 1]) +
            log(plot.data$x) * log(plot.data$x) %*% t(frac_coef_boot[, 2]) -
            reprow(
              log(ref) %*% t(frac_coef_boot[, 1]) +
                log(ref) * log(ref) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML == -1 && p2_ML != -1 && p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * log(plot.data$x) +
          beta[2] * plot.data$x^(p2_ML + 1) -
          (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (log(plot.data$x) - log(ref))^2 * cov[1, 1] +
              2 * (log(plot.data$x) - log(ref)) * (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1)) * cov[1, 2] +
              (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot[, 1]) +
            plot.data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
            reprow(
              log(ref) %*% t(frac_coef_boot[, 1]) +
                ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 && p2_ML == -1 && p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) +
          beta[2] * log(plot.data$x) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] +
              2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (log(plot.data) - log(ref)) * cov[1, 2] +
              (log(plot.data$x) - log(ref))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
            log(plot.data$x) %*% t(frac_coef_boot[, 2]) -
            reprow(
              ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                log(ref) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 && p2_ML != -1 && p1_ML == p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) +
          beta[2] * plot.data$x^(p2_ML + 1) * log(plot.data$x) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1) * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] +
              2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (plot.data$x^(p2_ML + 1) * log(plot.data$x) - ref^(p2_ML + 1) * log(ref)) * cov[1, 2] +
              (plot.data$x^(p2_ML + 1) * log(plot.data$x) - ref^(p2_ML + 1) * log(ref))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
            (plot.data$x^(p2_ML + 1) * log(plot.data$x)) %*% t(frac_coef_boot[, 2]) -
            reprow(
              ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                (ref^(p2_ML + 1) * log(ref)) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 && p2_ML != -1 && p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) +
          beta[2] * plot.data$x^(p2_ML + 1) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] +
              2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1)) * cov[1, 2] +
              (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
            plot.data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
            reprow(
              ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
    }
    if (family != "binomial") {
      figure <- ggplot(data = plot.data, mapping = aes(x = x)) +
        geom_hline(aes(yintercept = 0), colour = "grey") +
        geom_line(aes(y = yest), color = "black") +
        geom_line(aes(y = lci), color = "grey") +
        geom_line(aes(y = uci), color = "grey") +
        geom_point(
          mapping = aes(x = x, y = y), data = plot.data.ref,
          colour = "red", size = 4
        ) +
        theme_bw() +
        labs(
          x = pref_x,
          y = bquote(
            .(pref_y) ~
              " [" ~ .(pref_x_ref)["ref"] ~ "=" ~ .(round(ref, 2)) ~ "]"
          )
        ) +
        theme(
          axis.title = element_text(vjust = 0.5, size = 20),
          axis.text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      if (!is.null(breaks)) {
        suppressMessages(figure <- figure + scale_y_continuous(breaks = breaks))
      }
    }
    if (family == "binomial") {
      pref_y <- paste0("Odds ratio of ", pref_y)
      plot.data$yest <- exp(plot.data$yest)
      plot.data$uci <- exp(plot.data$uci)
      plot.data$lci <- exp(plot.data$lci)
      plot.data.ref$y <- exp(0)
      figure <- ggplot(data = plot.data, mapping = aes(x = x)) +
        geom_hline(aes(yintercept = 1), colour = "grey") +
        geom_line(aes(y = yest), color = "black") +
        geom_line(aes(y = lci), color = "grey") +
        geom_line(aes(y = uci), color = "grey") +
        geom_point(
          mapping = aes(x = x, y = y), data = plot.data.ref,
          colour = "red", size = 4
        ) +
        theme_bw() +
        labs(
          x = pref_x,
          y = bquote(
            .(pref_y) ~
              " [" ~ .(pref_x_ref)["ref"] ~ "=" ~ .(round(ref, 2)) ~ "]"
          )
        ) +
        theme(
          axis.title = element_text(vjust = 0.5, size = 20),
          axis.text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      if (!is.null(breaks)) {
        figure <- figure +
          scale_y_continuous(trans = "log", breaks = breaks)
      } else {
        ybreaks <- ggplot_build(figure)$layout$panel_params[[1]]$y$breaks
        figure <- figure +
          scale_y_continuous(trans = "log", breaks = ybreaks)
      }
    }
  }
  if (ci_type == "quantile") {
    xmin <- min(x)
    xmax <- max(x)
    prob <- 1 / ci_quantile
    x_ci <- quantile(x, probs = seq(0, 1, prob))
    x_quantiles_ci <- cut(x, x_ci, include.lowest = TRUE)
    x_quantiles_ci <- as.numeric(x_quantiles_ci)
    xmean_ci <- NULL
    for (i in 1:ci_quantile) {
      xmean_ci[i] <- mean(x[x_quantiles_ci == i])
    }
    plot.data <- data.frame(x = c(ref, xmean_ci))
    plot.data.est <- data.frame(
      x = runif(10000, quantile(x, probs = 0.001), quantile(x, probs = 0.999))
    )
    if (d == 1) {
      if (p_ML == -1) {
        plot.data$yest <- beta * log(plot.data$x) - (beta * log(ref))
        plot.data.est$yest <- beta * log(plot.data.est$x) - (beta * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((log(plot.data$x) - log(ref))^2 * cov)
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot) -
            reprow(log(ref) %*% t(frac_coef_boot), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p_ML != -1) {
        plot.data$yest <- beta * plot.data$x^(p_ML + 1) - beta * ref^(p_ML + 1)
        plot.data.est$yest <- beta * plot.data.est$x^(p_ML + 1) -
          beta * ref^(p_ML + 1)
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p_ML + 1) - ref^(p_ML + 1))^2 * cov
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p_ML + 1) %*% t(frac_coef_boot) -
            reprow(ref^(p_ML + 1) %*% t(frac_coef_boot), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
    }
    if (d == 2) {
      if (p1_ML == -1 && p2_ML == -1) {
        plot.data$yest <- beta[1] * log(plot.data$x) +
          beta[2] * log(plot.data$x) * log(plot.data$x) -
          (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
        plot.data.est$yest <- beta[1] * log(plot.data.est$x) +
          beta[2] * log(plot.data.est$x) * log(plot.data.est$x) -
          (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (log(plot.data$x) - log(ref))^2 * cov[1, 1] +
              2 * (log(plot.data$x) - log(ref)) * (log(plot.data$x) * log(plot.data$x) - log(ref) * log(ref)) * cov[1, 2] +
              (log(plot.data$x) * log(plot.data$x) - log(ref) * log(ref))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot[, 1]) +
            log(plot.data$x) * log(plot.data$x) %*% t(frac_coef_boot[, 2]) -
            reprow(
              log(ref) %*% t(frac_coef_boot[, 1]) +
                log(ref) * log(ref) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML == -1 && p2_ML != -1 && p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * log(plot.data$x) +
          beta[2] * plot.data$x^(p2_ML + 1) -
          (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
        plot.data.est$yest <- beta[1] * log(plot.data.est$x) +
          beta[2] * plot.data.est$x^(p2_ML + 1) -
          (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (log(plot.data$x) - log(ref))^2 * cov[1, 1] +
              2 * (log(plot.data$x) - log(ref)) * (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1)) * cov[1, 2] +
              (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot[, 1]) +
            plot.data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
            reprow(
              log(ref) %*% t(frac_coef_boot[, 1]) +
                ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 && p2_ML == -1 && p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) +
          beta[2] * log(plot.data$x) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] +
              2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (log(plot.data) - log(ref)) * cov[1, 2] +
              (log(plot.data$x) - log(ref))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
            log(plot.data$x) %*% t(frac_coef_boot[, 2]) -
            reprow(
              ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                log(ref) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 && p2_ML != -1 && p1_ML == p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) +
          beta[2] * plot.data$x^(p2_ML + 1) * log(plot.data$x) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1) * log(ref))
        plot.data.est$yest <- beta[1] * plot.data.est$x^(p1_ML + 1) +
          beta[2] * plot.data.est$x^(p2_ML + 1) * log(plot.data.est$x) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1) * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] +
              2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (plot.data$x^(p2_ML + 1) * log(plot.data$x) - ref^(p2_ML + 1) * log(ref)) * cov[1, 2] +
              (plot.data$x^(p2_ML + 1) * log(plot.data$x) - ref^(p2_ML + 1) * log(ref))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
            (plot.data$x^(p2_ML + 1) * log(plot.data$x)) %*% t(frac_coef_boot[, 2]) -
            reprow(
              ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                (ref^(p2_ML + 1) * log(ref)) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 && p2_ML != -1 && p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) +
          beta[2] * plot.data$x^(p2_ML + 1) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
        plot.data.est$yest <- beta[1] * plot.data.est$x^(p1_ML + 1) +
          beta[2] * plot.data.est$x^(p2_ML + 1) -
          (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt(
            (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] +
              2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1)) * cov[1, 2] +
              (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1))^2 * cov[2, 2]
          )
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
            plot.data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
            reprow(
              ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
              n = nrow(plot.data)
            )
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
    }
    highlight <- c("red", rep("black", (nrow(plot.data) - 1)))
    if (family != "binomial") {
      figure <- ggplot(data = plot.data, mapping = aes(x = x)) +
        geom_hline(aes(yintercept = 0), colour = "grey") +
        geom_line(aes(x = x, y = yest), color = "black", data = plot.data.est) +
        geom_errorbar(
          mapping = aes(x = x, ymin = lci, ymax = uci),
          color = "grey", width = 0.025
        ) +
        geom_point(aes(y = yest), color = highlight, size = 4) +
        theme_bw() +
        labs(
          x = pref_x,
          y = bquote(
            .(pref_y) ~
              " [" ~ .(pref_x_ref)["ref"] ~ "=" ~ .(round(ref, 2)) ~ "]"
          )
        ) +
        theme(
          axis.title = element_text(vjust = 0.5, size = 20),
          axis.text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      if (!is.null(breaks)) {
        figure <- figure +
          scale_y_continuous(breaks = breaks)
      }
    }
    if (family == "binomial") {
      pref_y <- paste0("Odds ratio of ", pref_y)
      plot.data$yest <- exp(plot.data$yest)
      plot.data$uci <- exp(plot.data$uci)
      plot.data$lci <- exp(plot.data$lci)
      plot.data.est$yest <- exp(plot.data.est$yest)
      figure <- ggplot(data = plot.data, mapping = aes(x = x)) +
        geom_hline(aes(yintercept = 1), colour = "grey") +
        geom_line(
          mapping = aes(x = x, y = yest), data = plot.data.est,
          color = "black"
        ) +
        geom_errorbar(
          mapping = aes(x = x, ymin = lci, ymax = uci),
          color = "grey", width = 0.025
        ) +
        geom_point(aes(y = yest), color = highlight, size = 4) +
        theme_bw() +
        labs(
          x = pref_x,
          y = bquote(
            .(pref_y) ~
              " [" ~ .(pref_x_ref)["ref"] ~ "=" ~ .(round(ref, 2)) ~ "]"
          )
        ) +
        theme(
          axis.title = element_text(vjust = 0.5, size = 20),
          axis.text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      if (!is.null(breaks)) {
        figure <- figure +
          scale_y_continuous(trans = "log", breaks = breaks)
      } else {
        ybreaks <- ggplot_build(figure)$layout$panel_params[[1]]$y$breaks
        figure <- figure +
          scale_y_continuous(trans = "log", breaks = ybreaks)
      }
    }
  }
  return(figure)

}
