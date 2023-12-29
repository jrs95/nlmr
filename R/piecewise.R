#' @title Piecewise linear Mendelian randomization
#'
#' @description `piecewise_mr` performs a Mendelian randomization (MR) analysis
#'   by fitting a piecewise linear function to localised average causal effects.
#'
#' @name piecewise_mr
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
#' @param nboot the number of bootstrap replications (default: `100`)
#'
#' @param fig a `logical` statement as to whether the user wants the results
#'   displayed in a figure (default: `TRUE`)
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
#' @param ci_quantiles the number of quantiles at which confidence intervals
#'   are to be displayed (default: `10`)
#'
#' @param breaks breaks on the y-axis of the figure
#'
#' @return `piecewise_mr` returns a `list` of non-linear MR results from the
#'   piecewise linear function MR approach:
#'
#' @return \item{n}{number of individuals}
#'
#' @return \item{model}{
#'   the model specifications:
#'     number of quantiles (`q`),
#'     number of bootstrap replications performed (`nboot`)
#' }
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
#'   y = y , x = x, g = g, covar = covar,
#'   family = "gaussian", q = 10, d = 1, ci = "model_se",
#'   fig = TRUE
#' )
#' summary(fp)
#' plm <- piecewise_mr(
#'   y = y , x = x, g = g, covar = covar,
#'   family = "gaussian", q = 10, nboot = 100,
#'   fig = TRUE
#' )
#' summary(plm)
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
piecewise_mr <- function(
  y, x, g, covar = NULL,
  family = "gaussian", q = 10, xpos = "mean", nboot = 100,
  fig = TRUE, ref = mean(x), pref_x = "x", pref_x_ref = "x", pref_y = "y",
  ci_quantiles = 10, breaks = NULL
) {

  # Errors
  if (!(is.vector(y) && is.vector(x) && is.vector(g)))
    stop("either the outcome, exposure or instrument is not a vector")
  if (!is.null(covar) && !is.data.frame(covar))
    stop("covar has to be a data.frame")
  if (
    !(
      (is.numeric(y) || is.integer(y)) &&
        (is.numeric(x) || is.integer(x)) &&
        (is.numeric(g) || is.integer(g))
    )
  )
    stop("either the outcome, exposure or instrument is not numeric")
  if (length(y) <= 1)
    stop("the outcome is less than or equal to a single value")
  if (length(y) != length(x) || length(y) != length(g))
    stop(
      "the number of observations for the outcome, exposure and instrument ",
      "are not all the same"
    )
  if (!is.null(covar) && nrow(covar) != length(y))
    stop(
      "the number of observations for the outcome and covariates are not ",
      "the same"
    )
  if (any(is.na(y)) || any(is.na(x)) || any(is.na(g)))
    stop(
      "there are missing values in either the outcome, exposure or ",
      "instrument"
    )
  if (!is.null(covar) && any(is.na(covar)))
    stop("there are missing values in the covariates")
  if (!(family == "gaussian" || family == "binomial"))
    stop("family has to be equal to either \"gaussian\" or \"binomial\"")
  if (family == "binomial" && any(!(y == 1 | y == 0)))
    stop("y has to be 0 or 1 if family is equal to \"binomial\"")
  if ((length(y) / 10) < q)
    stop("the quantiles should contain at least 10 observations")

  # Covariates
  if (!is.null(covar)) {
    covar <- model.matrix(as.formula(~.), data = covar)[, -1, drop = FALSE]
    if (any(is.na(covar)))
      stop("there are missing values in the covariates")
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
    method = "DL"
  )$pval[2]

  # Non-linearity tests
  p_quadratic <- rma(coef ~ xmean, (coef_se)^2, method = "FE")$pval[2]
  p_Q <- 1 - pchisq(rma(coef, vi = (coef_se)^2)$QE, df = (q - 1))

  # Results
  lci <- coef - 1.96 * coef_se
  uci <- coef + 1.96 * coef_se
  pval <- 2 * pnorm(-abs(coef / coef_se))

  # Figure
  if (fig == TRUE) {
    figure <- piecewise_figure(
      y = y, x = x, g = g, covar = covar,
      q = q, xcoef = xcoef, coef = coef, x0q = x0q,
      family = family, nboot = nboot, ref = ref,
      pref_x = pref_x, pref_x_ref = pref_x_ref, pref_y = pref_y,
      ci_quantiles = ci_quantiles, breaks = breaks
    )
  }

  # Output
  model <- as.matrix(data.frame(q = q, nboot = nboot))
  lace <- as.matrix(
    data.frame(beta = coef, se = coef_se, lci = lci, uci = uci, pval = pval)
  )
  rownames(lace) <- seq_len(nrow(lace))
  xcoef_quant <- as.matrix(data.frame(beta = xcoef_sub, se = xcoef_sub_se))
  rownames(xcoef_quant) <- seq_len(nrow(xcoef_quant))
  p_tests <- as.matrix(data.frame(quad = p_quadratic, Q = p_Q))
  p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
  if (fig == FALSE) {
    output <- list(
      n = length(y), model = model,
      lace = lace, xcoef = xcoef_quant,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity
    )
  } else {
    output <- list(
      n = length(y), model = model,
      lace = lace, xcoef = xcoef_quant,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity,
      figure = figure
    )
  }
  class(output) <- "piecewise_mr"

  # Return
  return(output)

}

#' @title Print piecewise_mr
#'
#' @description print method for class `"piecewise_mr"`.
#'
#' @param x an object of class `"piecewise_mr"`
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
print.piecewise_mr <- function(x, ...) {
  cat("\nCall: \npiecewise_mr")
  cat("\n\nCoefficients:\n")
  cat(x$lace[, 1])
  cat("\n\n")
  if (!is.null(x$figure)) {
    plot(x$figure)
  }
}

#' @title Summary of piecewise_mr
#'
#' @description summary method for class `"piecewise_mr"`.
#'
#' @param x an object of class `"piecewise_mr"`
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
summary.piecewise_mr <- function(x, ...) {
  model <- as.data.frame(x$model)
  n <- x$n
  coefficients <- as.data.frame(x$lace)
  p_tests <- as.data.frame(x$p_tests)
  p_heterogeneity <- as.data.frame(x$p_heterogeneity)
  if (is.null(x$figure)) {
    output <- list(
      model = model, n = n,
      coefficients = coefficients,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity
    )
  }
  if (!is.null(x$figure)) {
    output <- list(
      model = model, n = n,
      coefficients = coefficients,
      p_tests = p_tests, p_heterogeneity = p_heterogeneity,
      figure = x$figure
    )
  }
  class(output) <- "summary.piecewise_mr"
  return(output)
}

#' @title Print summary of piecewise_mr
#'
#' @description print.summary method for class `"piecewise_mr"`.
#'
#' @param x an object of class `"piecewise_mr"`.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
print.summary.piecewise_mr <- function(x, ...) {
  cat("Call: piecewise_mr")
  cat(
    "\n\nNumber of individuals: ", x$n, "; ",
    "Quantiles: ", as.character(x$model$q), "; ",
    "Number of bootstrap replications: ", as.character(x$model$nboot),
    sep = ""
  )
  cat("\n\nLACE:\n")
  names(x$coefficients) <- c(
    "Estimate", "Std. Error", "95%CI Lower", "95%CI Upper", "p.value"
  )
  printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
  cat("\nNon-linearity tests")
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

#' @title Piecewise linear figure
#'
#' @description `piecewise_figure` plots the piecewise linear function.
#'
#' @name piecewise_figure
#'
#' @import matrixStats
#'
#' @import ggplot2
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
#' @param xcoef the association between the exposure and the instrument
#'
#' @param coef the coefficients of the localised causal effects
#'
#' @param x0q quantiles of x0 (the IV-free exposure)
#'
#' @param family a description of the error distribution and link function
#'   to be used in the model and is a `character` string naming either the
#'   gaussian (i.e. `"gaussian"` for continuous data) or binomial
#'   (i.e. `"binomial"` for binary data) family function (default: `"gaussian"`)
#'
#' @param nboot the number of bootstrap replications (default: `100`)
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
#' @param ci_quantiles the number of quantiles at which confidence intervals
#'   are to be displayed (default: `10`)
#'
#' @param breaks breaks on the y-axis of the figure
#'
#' @return `piecewise_figure` returns a plot of the piecewise linear function.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
piecewise_figure <- function(
  y, x, g, covar = NULL,
  q = 10, xcoef, coef, x0q,
  family = "gaussian", nboot = 100, ref = mean(x),
  pref_x = "x", pref_x_ref = "x", pref_y = "y",
  ci_quantiles = 10, breaks = NULL
) {

  # Piecewise function
  m <- NULL
  l <- q + 1
  prob <- 1 / q
  quantiles_x <- quantile(x, probs = seq(0, 1, prob))
  x_quantiles <- cut(x, quantiles_x, include.lowest = TRUE)
  x_quantiles <- as.numeric(x_quantiles)
  for (i in 1:(l - 1)) {
    if (quantiles_x[i] <= ref && quantiles_x[(i + 1)] >= ref) {
      ref_pos <- i + 1
    }
  }
  for (k in 1:l) {
    if (k == 1) {
      m[k] <- quantile(x, probs = 0.0000000001)
    }
    if (k == l) {
      m[k] <- quantile(x, probs = 0.9999999999)
    }
    if (k > 1 && k < l) {
      m[k] <- max(x[x_quantiles == k - 1])
    }
  }
  y_mm <- NULL
  for (k in 1:l) {
    if (k == 1) {
      y_mm[k] <- 0
    }
    if (k >= 2) {
      y_mm[k] <- (coef[k - 1] * m[k] - coef[k - 1] * m[k - 1]) + y_mm[k - 1]
    }
    if (k == ref_pos) {
      y_ref <- (coef[k - 1] * ref - coef[k - 1] * m[k - 1]) + y_mm[k - 1]
    }
  }
  y_mm_ref <- y_mm - y_ref

  prob <- 1 / ci_quantiles
  ci_quant <- as.numeric(quantile(x, probs = seq(0, 1, prob)))
  x_quantiles_ci <- cut(x, ci_quant, include.lowest = TRUE)
  x_quantiles_ci <- as.numeric(x_quantiles_ci)
  y_mm_quant <- NULL
  xmean_quant <- NULL
  for (j in 1:ci_quantiles) {
    x_ci <- mean(x[x_quantiles_ci == j])
    xmean_quant[j] <- mean(x[x_quantiles_ci == j])
    for (i in 1:(l - 1)) {
      if (quantiles_x[i] <= x_ci && quantiles_x[(i + 1)] >= x_ci) {
        ci_pos <- i + 1
      }
    }
    for (k in 1:l) {
      if (k == ci_pos) {
        y_mm_quant[j] <- (coef[k - 1] * x_ci - coef[k - 1] * m[k - 1]) +
          y_mm[k - 1]
      }
    }
  }
  y_mm_quant_ref <- y_mm_quant - y_ref

  # Bootstrap CIs
  non_p_boot <- matrix(NA, nrow = nboot, ncol = q)
  for (i in 1:nboot) {
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
    for (j in 1:q) {
      if (!is.null(covar_boot)) {
        if (family == "gaussian") {
          non_p_boot[i, j] <- (
            lm(y_boot[x0qb == j] ~ g_boot[x0qb == j] + covar_boot[x0qb == j, , drop = FALSE])$coef[2]
          ) / xcoef
        }
        if (family == "binomial") {
          non_p_boot[i, j] <- (
            glm(y_boot[x0qb == j] ~ g_boot[x0qb == j] + covar_boot[x0qb == j, , drop = FALSE], family = binomial)$coef[2]
          ) / xcoef
        }
      } else {
        if (family == "gaussian") {
          non_p_boot[i, j] <- (
            lm(y_boot[x0qb == j] ~ g_boot[x0qb == j])$coef[2]
          ) / xcoef
        }
        if (family == "binomial") {
          non_p_boot[i, j] <- (
            glm(y_boot[x0qb == j] ~ g_boot[x0qb == j], family = binomial)$coef[2]
          ) / xcoef
        }
      }
    }
  }
  y_boot <- matrix(0, nrow = nrow(non_p_boot), ncol = l)
  y_mm_boot_ref <- matrix(0, nrow = nrow(non_p_boot), ncol = l)
  for (i in seq_len(nrow(non_p_boot))) {
    for (k in 1:l) {
      if (k == ref_pos) {
        y_boot_ref <- (non_p_boot[i, k - 1] * ref - non_p_boot[i, k - 1] * m[k - 1]) +
          y_boot[i, k - 1]
      }
    }
    for (k in 1:l) {
      if (k == 1) {
        y_boot[i, k] <- 0
      }
      if (k >= 2) {
        y_boot[i, k] <- (non_p_boot[i, k - 1] * m[k] - non_p_boot[i, k - 1] * m[k - 1]) +
          y_boot[i, k - 1]
      }
      if (k == ref_pos) {
        y_boot_ref <- (non_p_boot[i, k - 1] * ref - non_p_boot[i, k - 1] * m[k - 1]) +
          y_boot[i, k - 1]
      }
    }
    y_mm_boot_ref[i, ] <- y_boot[i, ] - y_boot_ref
  }
  y_quant_boot_ref <- matrix(NA, nrow = nrow(non_p_boot), ncol = ci_quantiles)
  for (j in 1:ci_quantiles) {
    x_ci <- mean(x[x_quantiles_ci == j])
    for (w in 1:(l - 1)) {
      if (quantiles_x[w] <= x_ci && quantiles_x[(w + 1)] >= x_ci) {
        ci_pos <- w + 1
      }
    }
    for (i in seq_len(nrow(non_p_boot))) {
      for (k in 1:l) {
        if (k == ref_pos) {
          y_boot_ref <- (non_p_boot[i, k - 1] * ref - non_p_boot[i, k - 1] * m[k - 1]) +
            y_boot[i, k - 1]
        }
      }
      for (k in 1:l) {
        if (k == ci_pos) {
          y_quant_boot <- (non_p_boot[i, k - 1] * x_ci - non_p_boot[i, k - 1] * m[k - 1]) +
            y_boot[i, k - 1]
        }
      }
      y_quant_boot_ref[i, j] <- y_quant_boot - y_boot_ref
    }
  }
  y_lci <- colQuantiles(y_quant_boot_ref, probs = 0.025)
  y_uci <- colQuantiles(y_quant_boot_ref, probs = 0.975)

  # Figure
  plot.data <- data.frame(x = m, y = y_mm_ref)
  plot.data.1 <- data.frame(
    x = xmean_quant, y = y_mm_quant_ref, y_lci = y_lci, y_uci = y_uci
  )
  plot.data.2 <- data.frame(x = ref, y = 0)
  if (family != "binomial") {
    figure <- ggplot(data = plot.data, mapping = aes(x)) +
      geom_hline(aes(yintercept = 0), colour = "grey") +
      geom_line(aes(x = x, y = y), colour = "black") +
      geom_errorbar(
        mapping = aes(x = x, ymin = y_lci, ymax = y_uci), data = plot.data.1,
        color = "grey", width = 0.025
      ) +
      geom_point(
        mapping = aes(x = x, y = y), data = plot.data.1,
        colour = "black", size = 3
      ) +
      geom_point(
        mapping = aes(x = x, y = y), data = plot.data.2,
        colour = "red", size = 3
      ) + theme_bw() +
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
    plot.data$y <- exp(plot.data$y)
    plot.data.1$y <- exp(plot.data.1$y)
    plot.data.1$y_lci <- exp(plot.data.1$y_lci)
    plot.data.1$y_uci <- exp(plot.data.1$y_uci)
    plot.data.2$y <- exp(plot.data.2$y)
    figure <- ggplot(data = plot.data, mapping = aes(x)) +
      geom_hline(aes(yintercept = 1), colour = "grey") +
      geom_line(aes(x = x, y = y), colour = "black") +
      geom_errorbar(
        mapping = aes(x = x, ymin = y_lci, ymax = y_uci),
        data = plot.data.1,
        color = "grey", width = 0.025
      ) +
      geom_point(
        mapping = aes(x = x, y = y), data = plot.data.1,
        colour = "black", size = 3
      ) +
      geom_point(
        mapping = aes(x = x, y = y), data = plot.data.2,
        colour = "red", size = 3
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

  # Return
  return(figure)

}
