#' Instrumental variable analysis using fractional polynomials
#'
#' frac_poly_mr performs instumental variable analysis by fitting fractional polynomial models to localised average causal effects using meta-regression.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param g the instrumental variable.
#' @param c matrix of covariates.
#' @param c_type vector of covariate types. These can either "numeric" or "factor" depending on whether the variables are continuous or categorical.
#' @param family a description of the error distribution and link function to be used in the model. For frac_poly_mr this can be a character string naming either the gaussian (i.e. for continuous data) or binomial (i.e. for binary data) family function.
#' @param q the number of quantiles the exposure distribution is to be split into. Within each quantile a causal effect will be fitted, known as a localised average causal effect (LACE). The default is deciles (i.e. 10 quantiles).
#' @param xpos the position used to relate x to the localised average causal effect. The default is the mean of the x-values within each quantile, otherwise specify a percentile (e.g. 0.5 corresponds to the median value).
#' @param method meta-regression method parsed to the rma package. The default is fixed-effects ("FE").
#' @param d fractional polynomial degree. The default is degree 1. The other options are: 1, 2, or "both".
#' @param pd p-value cut-off for choosing the best-fitting fractional polynomial of degree 2 over the best-fitting fractional polynomial degree 1. This option is only used if both . The default is 0.05.
#' @param ci the type of 95\% confidence interval. There are three options: (i) using the model standard errors ("model_se"), (ii) using bootstrap standard errors ("bootstrap_se"), (iii) using bootstrap percentile confidence intervals ("bootstrap_per"). The default is the model standard errors.
#' @param nboot the number of bootstrap replications (if required). The default is 100 replications.
#' @param fig a logical statement as to whether the user wants the results displayed in a figure. The default is false.
#' @param ref the reference point for the figure. This is the value of the function that represents the expected difference in the outcome compared with this reference value when the exposure is set to different values. The default is the mean of x.
#' @param pref_x the prefix/label for the x-axis. The default is "x".
#' @param pref_x_ref the prefix for the reference value displayed on the y-axis. The default is "x".
#' @param pref_y the prefix/label for the y-axis. The default is "y".
#' @param ci_type the type of confidence interval to be displayed on the graph. The default is "overall" where confidence intervals are presented as bands across the range of x. The alternative option is "quantile" where the confidence intervals are presented as error bars at the mean in each quantile of x.
#' @param ci_quantile the number of quantiles at which confidence intervals are to be displayed. The default is deciles (i.e. 10).
#' @param breaks breaks on the y-axis of the figure.
#' @return n number of individuals.
#' @return model the model specifications. The first column is the number of quantiles (q); the second column is the position used to relate x to the LACE in each quantiles (xpos); the third column is the type of confidence interval constructed (ci); the fourth column is the number of bootstrap replications performed (nboot).
#' @return powers the powers of the chosen polynomial.
#' @return coefficients the regression estimates. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se); the third column is the lower confidence interval (lci); the fourth column is the upper confidence interval (uci); the fifth column is the p-value (pval).
#' @return lace the localised average causal effect estimate in each quantile. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se); the third column is the lower confidence interval (lci); the fourth column is the upper confidence interval (uci); the fifth column is the p-value (pval).
#' @return xcoef the association between the instrument and the exposure in each quantile. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se).
#' @return p_tests the p-value of the non-linearity tests. The first column is the p-value of the test between the fractional polynomial degrees (fp_d1_d2); the second column is the p-value from the fractional polynomial non-linearity test (fp); the third column is the p-value from the quadratic test (quad); the fourth column is the p-value from the Cochran Q test (Q).
#' @return p_heterogeneity the p-value of heterogeneity. The first column is the p-value of the Cochran Q heterogeneity test (Q); the second column is the p-value from the trend test (trend).
#' @author James R Staley (js16174@bristol.ac.uk)
#' @export
frac_poly_mr <- function(y, x, g, c=NULL, c_type=NULL, family="gaussian", q=10, xpos="mean", method="FE", d=1, pd=0.05, ci="model_se", nboot=100, fig=F, ref=mean(x), pref_x="x", pref_x_ref="x", pref_y="y", ci_type="overall", ci_quantile=10, breaks=NULL){
  
  ##### Error messages #####
  if(!(is.vector(y) | is.vector(x) | is.vector(g))) stop('the outcome, exposure, and instrument are not all vectors')
  if(any(is.na(y)) | any(is.na(x)) | any(is.na(g)) | any(is.na(c))) stop('there are missing values in either the outcome, exposure, instrument or covariates')
  if(!(length(y)==length(x) & length(y)==length(g)) | (!is.null(c) & !(nrow(c)==length(y)))) stop('the number of observations for the outcome, exposure, instrument and covarites are different')
  if((!is.null(c) | !is.null(c_type)) & ncol(c)!=length(c_type)) stop('the number of columns of the covariates matrix does not match the number of covariate types')
  if(!(family=="gaussian" | family=="binomial")) stop('family has to be equal to either "gaussian" or "binomial"')
  if((length(y)/10)<q) stop('the quantiles should contain at least 10 observations')
  if(!(xpos=="mean" | (xpos>0 & xpos<1))) stop('the position used to relate x to the localised average causal effect')
  if(!(d==1 | d==2 | d=="both")) stop('the degree has to be equal to 1, 2 or "both"')
  if(!(ci=="model_se" | ci=="bootstrap_se" | ci=="bootstrap_per")) stop('the confidence intervals must be one of "model_se", "bootstrap_se" and "bootstrap_per"')
  
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
  frac_coef <- NULL
  frac_se <- NULL
  xmean<-NULL
  xcoef_sub <- NULL
  xcoef_sub_se <- NULL
  for(i in 1:q){
    if(family=="gaussian"){
      frac_coef[i] <- lm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,])$coef[2]
      frac_se[i] <- summary(lm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,]))$coef[2,2]
    }
    if(family=="binomial"){
      frac_coef[i] <- glm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], family=binomial)$coef[2]
      frac_se[i] <- summary(glm(y[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], family=binomial))$coef[2,2]
    }
    if(xpos=="mean"){xmean[i] <- mean(x[x0_quantiles==i])}
    if(xpos!="mean"){xmean[i] <- quantile(x[x0_quantiles==i], probs=xpos)}
    if(family=="gaussian"){
      xcoef_sub[i] <- lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,])$coef[2]
      xcoef_sub_se[i] <- summary(lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,]))$coef[2,2]
    }
    if(family=="binomial"){
      xcoef_sub[i] <- lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], subset=dataset$y==0)$coef[2]
      xcoef_sub_se[i] <- summary(lm(x[x0_quantiles==i]~g[x0_quantiles==i]+c1[x0_quantiles==i,]+c2[x0_quantiles==i,], subset=dataset$y==0))$coef[2,2]
    }
  }
  
  ##### Test of IV-exposure assumption #####
  p_het <- 1- pchisq(rma(xcoef_sub, vi=(xcoef_sub_se)^2)$QE, df=(q-1))
  p_het_trend <- rma.uni(xcoef_sub ~ xmean, vi=xcoef_sub_se^2, method=method)$pval[2]
  
  ##### Best-fitting fractional polynomial of degree 1 #####
  powers<-c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
  p<-NULL
  ML<-NULL
  j<-1
  for(p1 in powers){
    if(p1==-1){x1 <- xmean^p1}else{x1 <- (p1+1)*xmean^p1}
    p[j]<-p1
    ML[j] <- summary(rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method))$fit.stats[1,1]
    j<-j+1
  }
  fits <- data.frame(p, ML)
  fits$max <- 0
  fits$max[fits$ML==max(fits$ML)] <- 1
  p_ML <- fits$p[fits$max==1]
  
  ##### Best-fitting fractional polynomial of degree 2 #####
  if(d==1 | d==2 | d=="both"){
    powers1 <- c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
    powers2 <- c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
    p1<-NULL
    p2 <-NULL
    ML <- NULL
    j <- 1
    for(p11 in powers1){
      if(p11==-1){x1 <- xmean^p11}else{x1 <- (p11+1)*xmean^p11}
      for(p21 in powers2){
        if(p11==p21){if(p21==-1){x2 <- 2*(xmean^p21)*log(xmean)}else{x2 <- ((p21+1)*(xmean^p21)*log(xmean) + xmean^p21)}}
        else{if(p21==-1){x2 <- xmean^p21}else{x2 <- (p21+1)*xmean^p21}}
        p1[j]<-p11
        p2[j]<-p21
        cc <- try(rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method), silent=TRUE)
        if(is(cc, "try-error")==T){ML[j] <- NA}
        if(is(cc, "try-error")==F){ML[j] <- summary(rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method))$fit.stats[1,1]}
        j<-j+1
      }
      powers2<-powers2[-1]
    }
    fits <- data.frame(p1, p2, ML)
    fits$max <- 0
    fits$max[fits$ML==max(fits$ML, na.rm=T)] <- 1
    p1_ML <- fits$p1[fits$max==1]
    p2_ML <- fits$p2[fits$max==1]
  }
  
  ##### Best-fitting fractional polynomial of either degree 1 or degree 2 #####
  p_d1_d2 <- NA
  if(d==1 | d==2 | d=="both"){
    if(p_ML==-1){x1<-xmean^p_ML}else{x1 <- (p_ML+1)*xmean^p_ML}
    best_fracp_d1 <- rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method)
    dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2,1]
    if(p1_ML==-1){x1 <- xmean^p1_ML}else{x1 <- (p1_ML+1)*xmean^p1_ML}
    if(p1_ML==p2_ML){if(p2_ML==-1){x2 <- 2*(xmean^p2_ML)*log(xmean)}else{x2 <- ((p2_ML+1)*(xmean^p2_ML)*log(xmean) + xmean^p2_ML)}}
    else{if(p2_ML==-1){x2 <- xmean^p2_ML}else{x2 <- (p2_ML+1)*xmean^p2_ML}}
    best_fracp_d2 <- rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method)
    dev_best_fracp_d2 <- best_fracp_d2$fit.stats[2,1]
    p_d1_d2 <- 1 - pchisq((dev_best_fracp_d1 - dev_best_fracp_d2), df=2)
    if(p_d1_d2>=pd){d1 <- 1}else{d1 <- 2}
    if(d=="both"){d <- d1}
  }
  
  ##### Model #####
  if(d==1){if(p_ML==-1){x1<-xmean^p_ML}else{x1 <- (p_ML+1)*xmean^p_ML}; model <- rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method)}
  if(d==2){if(p1_ML==-1){x1<-xmean^p1_ML}else{x1 <- (p1_ML+1)*xmean^p1_ML}; if(p1_ML==p2_ML){if(p2_ML==-1){x2 <- 2*(xmean^p2_ML)*log(xmean)}else{x2 <- ((p2_ML+1)*(xmean^p2_ML)*log(xmean) + xmean^p2_ML)}}else{if(p2_ML==-1){x2 <- xmean^p2_ML}else{x2 <- (p2_ML+1)*xmean^p2_ML}}; model <- rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method)}
  
  ##### Bootstrap #####
  if(ci=="bootstrap_per" | ci=="bootstrap_se"){
    if(d==1){
      frac_coef_boot <- NULL
      for(i in 1:nboot){
        indices <- sample.int(nrow(data), size = nrow(data), replace = T)
        data1 <- data[indices,]
        c11 <- as.matrix(c1[indices,])
        c21 <- as.matrix(c2[indices,])
        x0_quantiles_boot <- x0_quantiles[indices]
        frac_coef1 <- NULL
        frac_se1 <- NULL
        xmean1 <- NULL
        for(j in 1:q){
          if(family=="gaussian"){
            frac_coef1[j] <- lm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,])$coef[2]
            frac_se1[j] <- summary(lm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,]))$coef[2,2]
          }
          if(family=="binomial"){
            frac_coef1[j] <- glm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,], family=binomial)$coef[2]
            frac_se1[j] <- summary(glm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,], family=binomial))$coef[2,2]
          }
          if(xpos=="mean"){xmean1[j] <- mean(data1$x[x0_quantiles_boot==j])} # x0_quantiles_boot==j could also be data1$x0_quantiles==j
          if(xpos!="mean"){xmean1[j] <- quantile(data1$x[x0_quantiles_boot==j], probs=xpos)}
        }
        if(p_ML==-1){x111<-xmean1^p_ML}else{x111 <- (p_ML+1)*xmean1^p_ML}
        mod <- rma.uni(frac_coef1/xcoef ~ -1 + x111, vi=(frac_se1/xcoef)^2, method=method)
        frac_coef_boot[i] <- mod$b[1]
      }
    }
    if(d==2){
      frac_coef_boot <- matrix(, nrow = nboot, ncol = 2)
      for(i in 1:nboot){
        indices <- sample.int(nrow(data), size = nrow(data), replace = T)
        data1 <- data[indices,]
        c11 <- as.matrix(c1[indices,])
        c21 <- as.matrix(c2[indices,])
        x0_quantiles_boot <- x0_quantiles[indices]
        frac_coef1 <- NULL
        frac_se1 <- NULL
        xmean1 <- NULL
        for(j in 1:q){
          if(family=="gaussian"){
            frac_coef1[j] <- lm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,])$coef[2]
            frac_se1[j] <- summary(lm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,]))$coef[2,2]
          }
          if(family=="binomial"){
            frac_coef1[j] <- glm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,], family=binomial)$coef[2]
            frac_se1[j] <- summary(glm(data1$y[x0_quantiles_boot==j]~data1$g[x0_quantiles_boot==j]+c11[x0_quantiles_boot==j,]+c21[x0_quantiles_boot==j,], family=binomial))$coef[2,2]
          }
          if(xpos=="mean"){xmean1[j] <- mean(data1$x[x0_quantiles_boot==j])} # x0_quantiles_boot==j could also be data1$x0_quantiles==j
          if(xpos!="mean"){xmean1[j] <- quantile(data1$x[x0_quantiles_boot==j], probs=xpos)}
        }
        if(p1_ML==-1){x111<-xmean1^p1_ML}else{x111 <- (p1_ML+1)*xmean1^p1_ML}
        if(p1_ML==p2_ML){if(p2_ML==-1){x211 <- 2*(xmean1^p2_ML)*log(xmean1)}else{x211 <- ((p2_ML+1)*(xmean1^p2_ML)*log(xmean1) + xmean1^p2_ML)}}
        else{if(p2_ML==-1){x211 <- xmean1^p2_ML}else{x211 <- (p2_ML+1)*xmean1^p2_ML}}
        mod <- rma.uni(frac_coef1/xcoef ~ -1 + x111 + x211, vi=(frac_se1/xcoef)^2, method=method)
        frac_coef_boot[i,1] <- mod$b[1]
        frac_coef_boot[i,2] <- mod$b[2]
      }
    }
  }
  
  ##### Fractional polynomial degree 1 test against linearity #####
  if(p_ML==-1){x1<-xmean^p_ML}else{x1 <- (p_ML+1)*xmean^p_ML}
  linear <- rma(frac_coef/xcoef ~ 1, vi=(frac_se/xcoef)^2, method=method)
  dev_linear <- linear$fit.stats[2,1]
  best_fracp_d1 <- rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method)
  dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2,1]
  p_fp <- 1 - pchisq((dev_linear - dev_best_fracp_d1), df=1)
  
  ##### Other tests #####
  p_quadratic <- rma(frac_coef/xcoef ~ xmean, vi=(frac_se/xcoef)^2, method=method)$pval[2]
  p_Q <- 1 - pchisq(rma(frac_coef/xcoef, vi=(frac_se/xcoef)^2)$QE, df=(q-1))
  
  ##### Results #####
  beta <- as.numeric(model$b)
  if(ci=="model_se"){if(d==1){powers <- p_ML + 1}; if(d==2){powers <- c(p1_ML, p2_ML); powers <- powers + 1}; cov <- model$vb; se <- model$se; lci <- beta - 1.96*se; uci <- beta + 1.96*se; pval <- 2*pnorm(-abs(beta/se))}
  if(ci=="bootstrap_se"){if(d==1){powers <- p_ML + 1; cov <- var(frac_coef_boot); se <- sqrt(cov)}; if(d==2){powers <- c(p1_ML, p2_ML); powers <- powers + 1; cov <- cov(frac_coef_boot); se <- sqrt(diag(cov))}; lci <- beta - 1.96*se; uci <- beta + 1.96*se; pval <- 2*pnorm(-abs(beta/se))}
  if(ci=="bootstrap_per"){if(d==1){powers <- p_ML + 1; se <- NA; lci <- quantile(frac_coef_boot, probs=0.025); uci <- quantile(frac_coef_boot, probs=0.975); pval <- NA}; if(d==2){powers <- c(p1_ML, p2_ML); powers <- powers + 1; se <- rep(NA, 2); lci <- NULL; uci <- NULL; pval <- NULL; lci[1] <- quantile(frac_coef_boot[,1], probs=0.025); lci[2] <- quantile(frac_coef_boot[,2], probs=0.025); uci[1] <- quantile(frac_coef_boot[,1], probs=0.975); uci[2] <- quantile(frac_coef_boot[,2], probs=0.975); pval <- rep(NA,2)}}
  lci <- as.numeric(lci); uci <- as.numeric(uci)
  if(ci=="model_se"){nboot<-NA}
  
  ##### Figure #####
  if(fig==T){
    if(ci_type=="overall"){
      plot.data <- data.frame(x=runif(10000, min(x), max(x)))
      plot.data.1 <- data.frame(x=ref, y=0)
      if(d==1){
        if(p_ML==-1){plot.data$yest <- beta*log(plot.data$x) - (beta*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((log(plot.data$x)-log(ref))^2*cov); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- log(plot.data$x)%*%t(frac_coef_boot) - reprow(log(ref)%*%t(frac_coef_boot), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p_ML!=-1){plot.data$yest <- beta*plot.data$x^(p_ML+1) - beta*ref^(p_ML+1); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p_ML+1)-ref^(p_ML+1))^2*cov); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p_ML+1)%*%t(frac_coef_boot) - reprow(ref^(p_ML+1)%*%t(frac_coef_boot), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
      }
      if(d==2){
        if(p1_ML==-1 & p2_ML==-1){plot.data$yest <- beta[1]*log(plot.data$x) + beta[2]*log(plot.data$x)*log(plot.data$x) - (beta[1]*log(ref) + beta[2]*log(ref)*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((log(plot.data$x)-log(ref))^2*cov[1,1] + 2*(log(plot.data$x)-log(ref))*(log(plot.data$x)*log(plot.data$x)-log(ref)*log(ref))*cov[1,2] + (log(plot.data$x)*log(plot.data$x)-log(ref)*log(ref))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- log(plot.data$x)%*%t(frac_coef_boot[,1]) + log(plot.data$x)*log(plot.data$x)%*%t(frac_coef_boot[,2]) - reprow(log(ref)%*%t(frac_coef_boot[,1]) + log(ref)*log(ref)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML==-1 & p2_ML!=-1 & p1_ML!=p2_ML){plot.data$yest <- beta[1]*log(plot.data$x) + beta[2]*plot.data$x^(p2_ML+1) - (beta[1]*log(ref) + beta[2]*ref^(p2_ML+1)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((log(plot.data$x)-log(ref))^2*cov[1,1] + 2*(log(plot.data$x)-log(ref))*(plot.data$x^(p2_ML+1)-ref^(p2_ML+1))*cov[1,2] + (plot.data$x^(p2_ML+1)-ref^(p2_ML+1))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- log(plot.data$x)%*%t(frac_coef_boot[,1]) + plot.data$x^(p2_ML+1)%*%t(frac_coef_boot[,2]) - reprow(log(ref)%*%t(frac_coef_boot[,1]) + ref^(p2_ML+1)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML!=-1 & p2_ML==-1 & p1_ML!=p2_ML){plot.data$yest <- beta[1]*plot.data$x^(p1_ML+1) + beta[2]*log(plot.data$x) - (beta[1]*ref^(p1_ML+1) + beta[2]*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p1_ML+1)-ref^(p1_ML+1))^2*cov[1,1] + 2*(plot.data$x^(p1_ML+1)-ref^(p1_ML+1))*(log(plot.data)-log(ref))*cov[1,2] + (log(plot.data$x)-log(ref))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p1_ML+1)%*%t(frac_coef_boot[,1]) + log(plot.data$x)%*%t(frac_coef_boot[,2]) - reprow(ref^(p1_ML+1)%*%t(frac_coef_boot[,1]) + log(ref)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML!=-1 & p2_ML!=-1 & p1_ML==p2_ML){plot.data$yest <- beta[1]*plot.data$x^(p1_ML+1) + beta[2]*plot.data$x^(p2_ML+1)*log(plot.data$x) - (beta[1]*ref^(p1_ML+1) + beta[2]*ref^(p2_ML+1)*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p1_ML+1)-ref^(p1_ML+1))^2*cov[1,1] + 2*(plot.data$x^(p1_ML+1)-ref^(p1_ML+1))*(plot.data$x^(p2_ML+1)*log(plot.data$x)-ref^(p2_ML+1)*log(ref))*cov[1,2] + (plot.data$x^(p2_ML+1)*log(plot.data$x)-ref^(p2_ML+1)*log(ref))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p1_ML+1)%*%t(frac_coef_boot[,1]) + (plot.data$x^(p2_ML+1)*log(plot.data$x))%*%t(frac_coef_boot[,2]) - reprow(ref^(p1_ML+1)%*%t(frac_coef_boot[,1]) + (ref^(p2_ML+1)*log(ref))%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML!=-1 & p2_ML!=-1 & p1_ML!=p2_ML){plot.data$yest <- beta[1]*plot.data$x^(p1_ML+1) + beta[2]*plot.data$x^(p2_ML+1) - (beta[1]*ref^(p1_ML+1) + beta[2]*ref^(p2_ML+1)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p1_ML+1)-ref^(p1_ML+1))^2*cov[1,1] + 2*(plot.data$x^(p1_ML+1)-ref^(p1_ML+1))*(plot.data$x^(p2_ML+1)-ref^(p2_ML+1))*cov[1,2] + (plot.data$x^(p2_ML+1)-ref^(p2_ML+1))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p1_ML+1)%*%t(frac_coef_boot[,1]) + plot.data$x^(p2_ML+1)%*%t(frac_coef_boot[,2]) - reprow(ref^(p1_ML+1)%*%t(frac_coef_boot[,1]) + ref^(p2_ML+1)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
      }
      if(family!="binomial"){figure <- ggplot(plot.data, aes(x=x)); figure <- figure + geom_hline(aes(yintercept=0), colour="grey") + geom_line(aes(y=yest), color="black") + geom_line(aes(y=lci), color="grey") + geom_line(aes(y=uci), color="grey") + theme_bw() + labs(x=pref_x,y=bquote(.(pref_y)~" ["~.(pref_x_ref)["ref"]~"="~.(round(ref,2))~"]")) + theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + geom_point(aes(x=x, y=y), data=plot.data.1, colour="red", size=4)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); if(!is.null(breaks)){suppressMessages(figure <- figure + scale_y_continuous(breaks=breaks))}}
      if(family=="binomial"){plot.data$yest <- exp(plot.data$yest); plot.data$uci <- exp(plot.data$uci); plot.data$lci <- exp(plot.data$lci); plot.data.1$y <- exp(0); figure <- ggplot(plot.data, aes(x=x)); figure <- figure + geom_hline(aes(yintercept=1), colour="grey") + geom_line(aes(y=yest), color="black") + geom_line(aes(y=lci), color="grey") + geom_line(aes(y=uci), color="grey") + theme_bw() + labs(x=pref_x,y=bquote(.(pref_y)~" ["~.(pref_x_ref)["ref"]~"="~.(round(ref,2))~"]")) + theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + geom_point(aes(x=x, y=y), data=plot.data.1, colour="red", size=4) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); if(!is.null(breaks)){figure <- figure + scale_y_continuous(breaks=breaks)}; figure <- figure + coord_trans(y="log")}
    }
    if(ci_type=="quantile"){
      xmin <- min(x)
      xmax <- max(x)
      prob <- (100/ci_quantile)/100
      x_ci <- quantile(x, probs=seq(0,1,prob))
      x_quantiles_ci <- cut(x, x_ci, include.lowest=T)
      x_quantiles_ci <- as.numeric(x_quantiles_ci)
      xmean_ci <- NULL
      for(i in 1:ci_quantile){xmean_ci[i] <- mean(x[x_quantiles_ci==i])}
      plot.data <- data.frame(x=c(ref, xmean_ci))
      plot.data.1 <- data.frame(x=runif(10000, quantile(x, probs=0.001), quantile(x, probs=0.999)))
      if(d==1){
        if(p_ML==-1){plot.data$yest <- beta*log(plot.data$x) - (beta*log(ref)); plot.data.1$yest <- beta*log(plot.data.1$x) - (beta*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((log(plot.data$x)-log(ref))^2*cov); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- log(plot.data$x)%*%t(frac_coef_boot) - reprow(log(ref)%*%t(frac_coef_boot), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p_ML!=-1){plot.data$yest <- beta*plot.data$x^(p_ML+1) - beta*ref^(p_ML+1); plot.data.1$yest <- beta*plot.data.1$x^(p_ML+1) - beta*ref^(p_ML+1); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p_ML+1)-ref^(p_ML+1))^2*cov); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p_ML+1)%*%t(frac_coef_boot) - reprow(ref^(p_ML+1)%*%t(frac_coef_boot), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
      }
      if(d==2){
        if(p1_ML==-1 & p2_ML==-1){plot.data$yest <- beta[1]*log(plot.data$x) + beta[2]*log(plot.data$x)*log(plot.data$x) - (beta[1]*log(ref) + beta[2]*log(ref)*log(ref)); plot.data.1$yest <- beta[1]*log(plot.data.1$x) + beta[2]*log(plot.data.1$x)*log(plot.data.1$x) - (beta[1]*log(ref) + beta[2]*log(ref)*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((log(plot.data$x)-log(ref))^2*cov[1,1] + 2*(log(plot.data$x)-log(ref))*(log(plot.data$x)*log(plot.data$x)-log(ref)*log(ref))*cov[1,2] + (log(plot.data$x)*log(plot.data$x)-log(ref)*log(ref))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- log(plot.data$x)%*%t(frac_coef_boot[,1]) + log(plot.data$x)*log(plot.data$x)%*%t(frac_coef_boot[,2]) - reprow(log(ref)%*%t(frac_coef_boot[,1]) + log(ref)*log(ref)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML==-1 & p2_ML!=-1 & p1_ML!=p2_ML){plot.data$yest <- beta[1]*log(plot.data$x) + beta[2]*plot.data$x^(p2_ML+1) - (beta[1]*log(ref) + beta[2]*ref^(p2_ML+1)); plot.data.1$yest <- beta[1]*log(plot.data.1$x) + beta[2]*plot.data.1$x^(p2_ML+1) - (beta[1]*log(ref) + beta[2]*ref^(p2_ML+1)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((log(plot.data$x)-log(ref))^2*cov[1,1] + 2*(log(plot.data$x)-log(ref))*(plot.data$x^(p2_ML+1)-ref^(p2_ML+1))*cov[1,2] + (plot.data$x^(p2_ML+1)-ref^(p2_ML+1))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- log(plot.data$x)%*%t(frac_coef_boot[,1]) + plot.data$x^(p2_ML+1)%*%t(frac_coef_boot[,2]) - reprow(log(ref)%*%t(frac_coef_boot[,1]) + ref^(p2_ML+1)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML!=-1 & p2_ML==-1 & p1_ML!=p2_ML){plot.data$yest <- beta[1]*plot.data$x^(p1_ML+1) + beta[2]*log(plot.data$x) - (beta[1]*ref^(p1_ML+1) + beta[2]*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p1_ML+1)-ref^(p1_ML+1))^2*cov[1,1] + 2*(plot.data$x^(p1_ML+1)-ref^(p1_ML+1))*(log(plot.data)-log(ref))*cov[1,2] + (log(plot.data$x)-log(ref))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p1_ML+1)%*%t(frac_coef_boot[,1]) + log(plot.data$x)%*%t(frac_coef_boot[,2]) - reprow(ref^(p1_ML+1)%*%t(frac_coef_boot[,1]) + log(ref)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML!=-1 & p2_ML!=-1 & p1_ML==p2_ML){plot.data$yest <- beta[1]*plot.data$x^(p1_ML+1) + beta[2]*plot.data$x^(p2_ML+1)*log(plot.data$x) - (beta[1]*ref^(p1_ML+1) + beta[2]*ref^(p2_ML+1)*log(ref)); plot.data.1$yest <- beta[1]*plot.data.1$x^(p1_ML+1) + beta[2]*plot.data.1$x^(p2_ML+1)*log(plot.data.1$x) - (beta[1]*ref^(p1_ML+1) + beta[2]*ref^(p2_ML+1)*log(ref)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p1_ML+1)-ref^(p1_ML+1))^2*cov[1,1] + 2*(plot.data$x^(p1_ML+1)-ref^(p1_ML+1))*(plot.data$x^(p2_ML+1)*log(plot.data$x)-ref^(p2_ML+1)*log(ref))*cov[1,2] + (plot.data$x^(p2_ML+1)*log(plot.data$x)-ref^(p2_ML+1)*log(ref))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p1_ML+1)%*%t(frac_coef_boot[,1]) + (plot.data$x^(p2_ML+1)*log(plot.data$x))%*%t(frac_coef_boot[,2]) - reprow(ref^(p1_ML+1)%*%t(frac_coef_boot[,1]) + (ref^(p2_ML+1)*log(ref))%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
        if(p1_ML!=-1 & p2_ML!=-1 & p1_ML!=p2_ML){plot.data$yest <- beta[1]*plot.data$x^(p1_ML+1) + beta[2]*plot.data$x^(p2_ML+1) - (beta[1]*ref^(p1_ML+1) + beta[2]*ref^(p2_ML+1)); plot.data.1$yest <- beta[1]*plot.data.1$x^(p1_ML+1) + beta[2]*plot.data.1$x^(p2_ML+1) - (beta[1]*ref^(p1_ML+1) + beta[2]*ref^(p2_ML+1)); if(ci!="bootstrap_per"){plot.data$yse <- sqrt((plot.data$x^(p1_ML+1)-ref^(p1_ML+1))^2*cov[1,1] + 2*(plot.data$x^(p1_ML+1)-ref^(p1_ML+1))*(plot.data$x^(p2_ML+1)-ref^(p2_ML+1))*cov[1,2] + (plot.data$x^(p2_ML+1)-ref^(p2_ML+1))^2*cov[2,2]); plot.data$lci <- plot.data$yest - 1.96*plot.data$yse; plot.data$uci <- plot.data$yest + 1.96*plot.data$yse}else{boot <- plot.data$x^(p1_ML+1)%*%t(frac_coef_boot[,1]) + plot.data$x^(p2_ML+1)%*%t(frac_coef_boot[,2]) - reprow(ref^(p1_ML+1)%*%t(frac_coef_boot[,1]) + ref^(p2_ML+1)%*%t(frac_coef_boot[,2]), n=nrow(plot.data)); plot.data$lci <- rowQuantiles(boot, probs=0.025); plot.data$uci <- rowQuantiles(boot, probs=0.975)}}
      }
      highlight <- c("red", rep("black", (nrow(plot.data)-1)))
      if(family!="binomial"){figure <- ggplot(plot.data, aes(x=x)); figure <- figure + geom_hline(aes(yintercept=0), colour="grey") + geom_line(aes(x=x, y=yest), color="black", data=plot.data.1) + geom_errorbar(mapping=aes(x=x, ymin=lci, ymax=uci), color="grey", width=0.025) + geom_point(aes(y=yest), color=highlight, size=4) + theme_bw() + labs(x=pref_x,y=bquote(.(pref_y)~" ["~.(pref_x_ref)["ref"]~"="~.(round(ref,2))~"]")) + theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); if(!is.null(breaks)){figure <- figure + scale_y_continuous(breaks=breaks)}}
      if(family=="binomial"){plot.data$yest <- exp(plot.data$yest); plot.data$uci <- exp(plot.data$uci); plot.data$lci <- exp(plot.data$lci); figure <- ggplot(plot.data, aes(x=x)); plot.data.1$yest <- exp(plot.data.1$yest); figure <- figure + geom_hline(aes(yintercept=1), colour="grey") + geom_line(aes(x=x, y=yest), color="black", data=plot.data.1) + geom_errorbar(mapping=aes(x=x, ymin=lci, ymax=uci), color="grey", width=0.025) + geom_point(aes(y=yest), color=highlight, size=4) + theme_bw() + labs(x=pref_x,y=bquote(.(pref_y)~" ["~.(pref_x_ref)["ref"]~"="~.(round(ref,2))~"]")) + theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); if(!is.null(breaks)){figure <- figure + scale_y_continuous(breaks=breaks)}; figure <- figure + coord_trans(y="log")}
    }
  }
  
  ##### Return #####
  model <- as.matrix(data.frame(q=q, xpos=xpos, ci_type=ci, nboot=nboot))
  coefficients <- as.matrix(data.frame(beta=beta, se=se, lci=lci, uci=uci, pval=pval))
  rownames(coefficients) <- powers
  lace <- as.matrix(data.frame(beta=(frac_coef/xcoef), se=(abs(frac_se/xcoef)), lci=(frac_coef/xcoef - 1.96*(abs(frac_se/xcoef))), uci=(frac_coef/xcoef + 1.96*(abs(frac_se/xcoef))), pval=(2*pnorm(-abs(frac_coef/frac_se)))))
  rownames(lace) <- 1:nrow(lace)
  xcoef_quant <- as.matrix(data.frame(beta=xcoef_sub, se=xcoef_sub_se))
  rownames(xcoef_quant) <- 1:nrow(xcoef_quant)
  p_tests <- as.matrix(data.frame(fp_d1_d2=p_d1_d2, fp=p_fp, quad=p_quadratic, Q=p_Q))
  p_heterogeneity <- as.matrix(data.frame(Q=p_het, trend=p_het_trend))
  if(fig==F){results <- list(n=N, model=model, powers=powers, coefficients=coefficients, lace=lace, xcoef=xcoef_quant, p_tests=p_tests, p_heterogeneity=p_heterogeneity)}
  else{results <- list(n=N, model=model, powers=powers, coefficients=coefficients, lace=lace, xcoef=xcoef_quant, p_tests=p_tests, p_heterogeneity=p_heterogeneity, figure=figure)}
  class(results) <- "frac_poly_mr"
  return(results)
}


#' Print Fractional Polynomial Fits
#'
#' print method for class "frac_poly_mr".
#' @param x an object of class "frac_poly_mr".
#' @author James R Staley (js16174@bristol.ac.uk)
#' @export
print.frac_poly_mr <- function(x, ...){
  cat("\nCall: \nfrac_poly_mr")
  cat("\n\nPowers:\n")
  cat(x$powers)
  cat("\n\nCoefficients:\n")
  cat(x$coefficients[,1])
  cat("\n\n")
}


#' Summarizing Fractional Polynomial Fits
#'
#' summary method for class "frac_poly_mr".
#' @param x an object of class "frac_poly_mr".
#' @author James R Staley (js16174@bristol.ac.uk)
#' @export
summary.frac_poly_mr <- function(x, ...){
  model <- as.data.frame(x$model)
  powers <- x$powers
  n <- x$n
  coefficients <- as.data.frame(x$coefficients)
  if(model$ci_type=="bootstrap_per"){coefficients <- coefficients[,c(1,3,4)]}
  #coefficients$ci <- paste0("(",format(coefficients$lci, 7),", ",format(coefficients$lci, 7),")")
  p_tests <- as.data.frame(x$p_tests)
  p_heterogeneity <- as.data.frame(x$p_heterogeneity)
  summ <- list(model=model, powers=powers, n=n, coefficients=coefficients, p_tests=p_tests, p_heterogeneity=p_heterogeneity)
  class(summ) <- "summary.frac_poly_mr"
  return(summ)
}


#' Print Summary Fractional Polynomial Fits
#'
#' print.summary method for class "frac_poly_mr".
#' @param x an object of class "frac_poly_mr".
#' @author James R Staley (js16174@bristol.ac.uk)
#' @export
print.summary.frac_poly_mr <- function(x, ...){
  cat("Call: frac_poly_mr")
  ci_type <- "Model based SEs"
  if(x$model$ci_type=="bootstrap_se"){ci_type <- "Bootstrap based SEs"}
  if(x$model$ci_type=="bootstrap_per"){ci_type <- "Percentile bootstrap"}
  if(ci_type=="Model based SEs"){cat("\n\nNumber of individuals: ", x$n,"; Quantiles: ", as.character(x$model$q), "; 95%CI: ", ci_type, sep="")}
  if(ci_type!="Model based SEs"){cat("\n\nNumber of individuals: ", x$n,"; Quantiles: ", as.character(x$model$q), "; 95%CI: ", ci_type, "; Number of bootstrap replications: ", as.character(x$model$nboot), sep="")}
  cat("\n\nPowers:", x$powers)
  cat("\n\nCoefficients:\n")
  if(ci_type=="Percentile bootstrap"){names(x$coefficients) <- c("Estimate", "95%CI Lower", "95%CI Upper"); printCoefmat(x$coefficients)}
  if(ci_type!="Percentile bootstrap"){names(x$coefficients) <- c("Estimate", "Std. Error", "95%CI Lower", "95%CI Upper", "p.value"); printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)}
  cat("\nNon-linearity tests")
  cat("\nFractional polynomial degree p-value:", signif(x$p_tests$fp_d1_d2, digits=3))
  cat("\nFractional polynomial non-linearity p-value:", signif(x$p_tests$fp, digits=3))
  cat("\nQuadratic p-value:", signif(x$p_tests$quad, digits=3))
  cat("\nCochran Q p-value:", signif(x$p_tests$Q, digits=3))
  cat("\n\nHeterogeneity tests")
  cat("\nCochran Q p-value:", signif(x$p_heterogeneity$Q, digits=3))
  cat("\nTrend p-value:", signif(x$p_heterogeneity$trend, digits=3))
}