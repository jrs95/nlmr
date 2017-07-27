#' Instrumental variable analysis using a piecewise linear function
#'
#' piecewise_mr performs instumental variable analysis by fitting a piecewise linear function to localised average causal effects.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param c matrix of covariates.
#' @param c_type vector of covariate types. These can either "numeric" or "factor" depending on whether the variables are continuous or categorical.
#' @param family a description of the error distribution and link function to be used in the model. For frac_poly_mr this can be a character string naming either the gaussian (i.e. for continuous data) or binomial (i.e. for binary data) family function.
#' @param q the number of quantiles the exposure distribution is to be split into. Within each quantile a causal effect will be fitted, known as a localised average causal effect (LACE). The default is deciles (i.e. 10 quantiles).
#' @param nboot the number of bootstrap replications (if required). The default is 100 replications.
#' @param fig a logical statement as to whether the user wants the results displayed in a figure. The default is false.
#' @param ref the reference point for the figure. This is the value of the function that represents the expected difference in the outcome compared with this reference value when the exposure is set to different values. The default is the mean of x.
#' @param pref_x the prefix/label for the x-axis. The default is "x".
#' @param pref_x_ref the prefix for the reference value displayed on the y-axis. The default is "x".
#' @param pref_y the prefix/label for the y-axis. The default is "y".
#' @param ci_quantile the number of quantiles at which confidence intervals are to be displayed. The default is deciles (i.e. 10).
#' @param breaks breaks on the y-axis of the figure.
#' @return n number of individuals.
#' @return model the model specifications. The first column is the number of quantiles (q); the second column is the number of bootstrap replications performed (nboot).
#' @return coefficients the regression estimates. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se); the third column is the lower confidence interval (lci); the fourth column is the upper confidence interval (uci); the fifth column is the p-value (pval).
#' @return lace the localised average causal effect estimate in each quantile. The first column is the regression coefficients (beta); the second column is the standard errors of the regression coefficients (se).
#' @return xcoef the association between the instrument and the exposure in each quantile. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se).
#' @return p_tests the p-value of the non-linearity tests. The first column is the p-value of the test between the fractional polynomial degrees (fp_d1_d2); the second column is the p-value from the fractional polynomial non-linearity test (fp); the third column is the p-value from the quadratic test (quad); the fourth column is the p-value from the Cochran Q test (Q).
#' @return p_heterogeneity the p-value of heterogeneity. The first column is the p-value of the Cochran Q heterogeneity test (Q); the second column is the p-value from the trend test (trend).
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
piecewise_mr <- function(y, x, g, c=NULL, c_type=NULL, family="gaussian", q=10, nboot=100, fig=T, ref=mean(x), pref_x="x", pref_x_ref="x", pref_y="y", ci_quantiles=10, breaks=NULL){

  ##### Error messages #####
  if(!(is.vector(y) | is.vector(x) | is.vector(g))) stop('either the outcome, exposure or instrument is not a vector')
  if(length(y)<=1) stop('the outcome is less than or equal to a single value')
  if(!(length(y)==length(x) & length(y)==length(g)) | (if(!is.null(c)){(nrow(c)!=length(y))}else{FALSE})) stop('the number of observations for the outcome, exposure, instrument and covarites are different')
  if(any(is.na(y)) | any(is.na(x)) | any(is.na(g)) | (if(!is.null(c)){any(is.na(c))}else{FALSE})) stop('there are missing values in either the outcome, exposure, instrument or covariates')
  if((!is.null(c) | !is.null(c_type)) & (if(!is.null(c)){ncol(c)!=length(c_type)}else{FALSE})) stop('the number of columns of the covariates matrix does not match the number of covariate types')
  if(!(family=="gaussian" | family=="binomial")) stop('family has to be equal to either "gaussian" or "binomial"')
  if((length(y)/10)<q) stop('the quantiles should contain at least 10 observations')
  
  ##### Covariates #####
  c1 <- c[,c_type!="factor"]
  if(length(c1)>0){c1 <- as.matrix(as.data.frame(c1)); class(c1) <- "numeric"}
  c2.1 <- as.data.frame(c[,c_type=="factor"])
  c2 <- data.frame(id=1:length(y))
  if(length(c2.1)>0){
    for(i in 1:ncol(c2.1)){
      cat <- model.matrix(~cat, data=data.frame(cat=as.factor(c2.1[,i])))
      cat <- cat[,2:ncol(cat)]
      c2 <- cbind(c2, cat)
    }
    c2 <- c2[,2:ncol(c2)]
    c2 <- as.matrix(as.data.frame(c2))
  }
  if(length(c1)==0){c1 <- data.frame(c1=(rep(1,length(y)))); c1 <- as.matrix(c1); class(c1) <- "numeric"}
  if(length(c2.1)==0){c2 <- data.frame(c2=(rep(1,length(y)))); c2 <- as.matrix(c2)}
  
  ##### x0 (IV-Free) #####
  if(family=="gaussian"){x0 <- resid(lm(x~g+c1+c2)); xcoef <- lm(x~g+c1+c2)$coef[2]}
  if(family=="binomial"){dataset <- data.frame(y=y, x=x, g=g); dataset <- cbind(dataset, c1); dataset <- cbind(dataset, c2); model <- lm(x~g+c1+c2, subset=dataset$y==0); x0 <- x - predict(model, newdata=dataset); xcoef <- summary(model)$coefficients[2,1]}
  prob <- (100/q)/100
  quantiles <- quantile(x0, probs=seq(0,1, prob))
  x0_quantiles <- cut(x0, quantiles, include.lowest=T)
  x0_quantiles <- as.numeric(x0_quantiles)
  data <- data.frame(y,x,g,x0,x0_quantiles)
  N <- nrow(data)
  
  ##### LACE in each quantile #####
  xcoef <- lm(x~g+c1+c2)$coef[2]
  coef <- NULL
  se <- NULL
  xmean<-NULL
  xcoef_sub <- NULL
  xcoef_sub_se <- NULL
  for(i in 1:q){
    if(family=="gaussian"){
      coef[i] <- lm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], data=data)$coef[2]/xcoef
      se[i] <- sqrt((summary(lm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], data=data))$coef[2,2]/xcoef)^2)
    }
    if(family=="binomial"){
      coef[i] <- glm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], family=binomial, data=data)$coef[2]/xcoef
      se[i] <- sqrt((summary(glm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], family=binomial, data=data))$coef[2,2]/xcoef)^2)
    }
    xmean[i] <- mean(x[x0_quantiles==i])
    if(family=="gaussian"){
      xcoef_sub[i] <- lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], data=data)$coef[2]
      xcoef_sub_se[i] <- summary(lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], data=data))$coef[2,2]
    }
    if(family=="binomial"){
      xcoef_sub[i] <- lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], subset=dataset$y==0, data=data)$coef[2]
      xcoef_sub_se[i] <- summary(lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], subset=dataset$y==0, data=data))$coef[2,2]
    }
  }
  
  ##### Test of IV-exposure assumption #####
  p_het <- 1- pchisq(rma(xcoef_sub, vi=(xcoef_sub_se)^2)$QE, df=(q-1))
  p_het_trend <- rma.uni(xcoef_sub ~ xmean, vi=xcoef_sub_se^2, method="DL")$pval[2]
  
  ##### Tests #####
  p_quadratic <- rma(coef/xcoef ~ xmean, (se/xcoef)^2, method="FE")$pval[2]
  p_Q <- 1 - pchisq(rma(coef/xcoef, vi=(se/xcoef)^2)$QE, df=(q-1))
  
  ##### Results #####
  lci <- coef - 1.96*se
  uci <- coef + 1.96*se
  pval <- 2*pnorm(-abs(coef/se))
  
  ##### Figure #####
  if(fig==T){
    
    ## Main Effect
    m <- NULL
    l <- q + 1
    quantiles_x <- quantile(x, probs= seq(0,1, prob))
    x_quantiles <- cut(x, quantiles_x, include.lowest=T)
    x_quantiles <- as.numeric(x_quantiles)
    for(i in 1:(l-1)){if(quantiles_x[i]<=ref & quantiles_x[(i+1)]>=ref){ref_pos <- i+1}}
    for(k in 1:l){
      if(k==1){m[k] <- quantile(x, probs=0.0000000001)}
      if(k==l){m[k] <- quantile(x, probs=0.9999999999)}
      if(k>1 & k<l){m[k] <- max(x[x_quantiles==k-1])}
    }
    y_mm<-NULL
    for(k in 1:l){
      if(k==1){y_mm[k] <- 0}
      if(k>=2){y_mm[k] <- (coef[k-1]*m[k] - coef[k-1]*m[k-1]) + y_mm[k-1]}
      if(k==ref_pos){y_ref <- (coef[k-1]*ref - coef[k-1]*m[k-1]) + y_mm[k-1]}
    }
    y_mm_ref <- y_mm - y_ref
    
    prob <- (100/ci_quantiles)/100
    ci_quant <- as.numeric(quantile(x, probs=seq(0,1,prob)))
    x_quantiles_ci <- cut(x, ci_quant, include.lowest=T)
    x_quantiles_ci <- as.numeric(x_quantiles_ci)
    y_mm_quant<-NULL
    xmean_quant <- NULL
    for(j in 1:ci_quantiles){
      x_ci <- mean(x[x_quantiles_ci==j])
      xmean_quant[j] <- mean(x[x_quantiles_ci==j])
      for(i in 1:(l-1)){if(quantiles_x[i]<=x_ci & quantiles_x[(i+1)]>=x_ci){ci_pos <- i+1}}
      for(k in 1:l){
        if(k==ci_pos){y_mm_quant[j] <- (coef[k-1]*x_ci - coef[k-1]*m[k-1]) + y_mm[k-1]}
      }
    }
    y_mm_quant_ref <- y_mm_quant - y_ref
    
    ## Bootstrap CIs
    non_p_boot <- matrix(, nrow = nboot, ncol = q)
    for(i in 1:nboot){
      indices <- sample.int(nrow(data), size = nrow(data), replace = T)
      data1 <- data[indices,]
      c11 <- as.matrix(c1[indices,])
      c21 <- as.matrix(c2[indices,])
      x0_quantiles_boot <- x0_quantiles[indices]
      for(j in 1:q){
        if(family=="gaussian"){non_p_boot[i,j] <- (lm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,])$coef[2])/xcoef}
        if(family=="binomial"){non_p_boot[i,j] <- (glm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,], family="binomial")$coef[2])/xcoef}
      }
    }
    y_boot <- matrix(, nrow = nrow(non_p_boot), ncol = l)
    y_mm_boot_ref <- matrix(, nrow = nrow(non_p_boot), ncol = l)
    for(i in 1:nrow(non_p_boot)){
      for(k in 1:l){
        if(k==ref_pos){y_boot_ref <- (non_p_boot[i,k-1]*ref - non_p_boot[i,k-1]*m[k-1]) + y_boot[i,k-1]}
      }
      for(k in 1:l){
        if(k==1){y_boot[i,k] <- 0}
        if(k>=2){y_boot[i,k] <- (non_p_boot[i,k-1]*m[k] - non_p_boot[i,k-1]*m[k-1]) + y_boot[i,k-1]}
        if(k==ref_pos){y_boot_ref <- (non_p_boot[i,k-1]*ref - non_p_boot[i,k-1]*m[k-1]) + y_boot[i,k-1]}
      }
      y_mm_boot_ref[i,] <- y_boot[i,] - y_boot_ref
    }
    y_quant_boot_ref <- matrix(, nrow = nrow(non_p_boot), ncol = ci_quantiles)
    for(j in 1:ci_quantiles){
      x_ci <- mean(x[x_quantiles_ci==j])
      for(w in 1:(l-1)){if(quantiles_x[w]<=x_ci & quantiles_x[(w+1)]>=x_ci){ci_pos <- w+1}}
      for(i in 1:nrow(non_p_boot)){
        for(k in 1:l){
          if(k==ref_pos){y_boot_ref <- (non_p_boot[i,k-1]*ref - non_p_boot[i,k-1]*m[k-1]) + y_boot[i,k-1]}
        }
        for(k in 1:l){
          if(k==ci_pos){y_quant_boot <- (non_p_boot[i,k-1]*x_ci - non_p_boot[i,k-1]*m[k-1]) + y_boot[i,k-1]}
        }
        y_quant_boot_ref[i,j] <- y_quant_boot - y_boot_ref
      }
    }
    y_lci <- colQuantiles(y_quant_boot_ref, probs=0.025)
    y_uci <- colQuantiles(y_quant_boot_ref, probs=0.975)
    
    ## Figure
    plot.data <- data.frame(x=m,y=y_mm_ref)
    plot.data.1 <- data.frame(x=xmean_quant, y=y_mm_quant_ref, y_lci=y_lci, y_uci=y_uci)
    plot.data.2 <- data.frame(x=ref, y=0)
    if(family!="binomial"){figure <- ggplot(plot.data, aes(x)); figure <- figure + geom_hline(aes(yintercept=0), colour="grey") + geom_line(aes(x=x, y=y), colour="black"); figure <- figure + geom_errorbar(aes(x=x, ymin=y_lci, ymax=y_uci), data=plot.data.1, color="grey", width = 0.025) + geom_point(aes(x=x, y=y), data=plot.data.1, colour="black", size=3) + geom_point(aes(x=x, y=y), data=plot.data.2, colour="red", size=3); figure <- figure + theme_bw() + labs(x=pref_x,y=bquote(.(pref_y)~" ["~.(pref_x_ref)["ref"]~"="~.(round(ref,2))~"]")) + theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); if(!is.null(breaks)){figure <- figure + scale_y_continuous(breaks=breaks)}}
    if(family=="binomial"){plot.data$y <- exp(plot.data$y); plot.data.1$y <- exp(plot.data.1$y); plot.data.1$y_lci <- exp(plot.data.1$y_lci); plot.data.1$y_uci <- exp(plot.data.1$y_uci); plot.data.2$y <- exp(plot.data.2$y); figure <- ggplot(plot.data, aes(x)); figure <- figure + geom_hline(aes(yintercept=1), colour="grey") + geom_line(aes(x=x, y=y), colour="black"); figure <- figure + geom_errorbar(aes(x=x, ymin=y_lci, ymax=y_uci), data=plot.data.1, color="grey", width = 0.025) + geom_point(aes(x=x, y=y), data=plot.data.1, colour="black", size=3) + geom_point(aes(x=x, y=y), data=plot.data.2, colour="red", size=3); figure <- figure + theme_bw() + labs(x=pref_x,y=bquote(.(pref_y)~" ["~.(pref_x_ref)["ref"]~"="~.(round(ref,2))~"]")) + theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); if(!is.null(breaks)){figure <- figure + scale_y_continuous(breaks=breaks)}; figure <- figure  + coord_trans(y="log")}
  }
  
  ##### Return #####
  model <- as.matrix(data.frame(q=q, nboot=nboot))
  lace <- as.matrix(data.frame(beta=coef, se=se, lci=lci, uci=uci, pval=pval))
  rownames(lace) <- 1:nrow(lace)
  xcoef_quant <- as.matrix(data.frame(beta=xcoef_sub, se=xcoef_sub_se))
  rownames(xcoef_quant) <- 1:nrow(xcoef_quant)
  p_tests <- as.matrix(data.frame(quad=p_quadratic, Q=p_Q))
  p_heterogeneity <- as.matrix(data.frame(Q=p_het, trend=p_het_trend))
  if(fig==F){results <- list(n=N, model=model, lace=lace, xcoef=xcoef_quant, p_tests=p_tests, p_heterogeneity=p_heterogeneity)}
  else{results <- list(n=N, model=model, lace=lace, xcoef=xcoef_quant, p_tests=p_tests, p_heterogeneity=p_heterogeneity, figure=figure)}
  class(results) <- "piecewise_mr"
  return(results)
}


#' Print Piecewise Linear Fits
#'
#' print method for class "piecewise_mr".
#' @param x an object of class "piecewise_mr".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
print.piecewise_mr <- function(x, ...){
  cat("\nCall: \npiecewise_mr")
  cat("\n\nCoefficients:\n")
  cat(x$lace[,1])
  cat("\n\n")
  if(!is.null(x$figure)){plot(x$figure)}
}


#' Summarizing Piecewise Linear Fits
#'
#' summary method for class "piecewise_mr".
#' @param x an object of class "piecewise_mr".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
summary.piecewise_mr <- function(x, ...){
  model <- as.data.frame(x$model)
  n <- x$n
  coefficients <- as.data.frame(x$lace)
  #coefficients$ci <- paste0("(",coefficients$lci,", ",coefficients$uci,")")
  p_tests <- as.data.frame(x$p_tests)
  p_heterogeneity <- as.data.frame(x$p_heterogeneity)
  if(is.null(x$figure)){summ <- list(model=model, n=n, coefficients=coefficients, p_tests=p_tests, p_heterogeneity=p_heterogeneity)}
  if(!is.null(x$figure)){summ <- list(model=model, n=n, coefficients=coefficients, p_tests=p_tests, p_heterogeneity=p_heterogeneity, figure=x$figure)}
  class(summ) <- "summary.piecewise_mr"
  return(summ)
}


#' Print Summary Piecewise Linear Fits
#'
#' print summary method for class "piecewise_mr".
#' @param x an object of class "piecewise_mr".
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
print.summary.piecewise_mr <- function(x, ...){
  cat("Call: piecewise_mr")
  cat("\n\nNumber of individuals: ",x$n,"; Quantiles: ", as.character(x$model$q), "; Number of bootstrap replications: ", as.character(x$model$nboot), sep="")
  cat("\n\nLACE:\n")
  names(x$coefficients) <- c("Estimate", "Std. Error", "95%CI Lower", "95%CI Upper", "p.value")
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
  cat("\nNon-linearity tests")
  cat("\nQuadratic p-value:", signif(x$p_tests$quad, digits=3))
  cat("\nCochran Q p-value:", signif(x$p_tests$Q, digits=3))
  cat("\n\nHeterogeneity tests")
  cat("\nCochran Q p-value:", signif(x$p_heterogeneity$Q, digits=3))
  cat("\nTrend p-value:", signif(x$p_heterogeneity$trend, digits=3))
  cat("\n")
  if(!is.null(x$figure)){plot(x$figure)}
}
