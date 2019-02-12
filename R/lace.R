#' Repeat rows
#'
#' This function creates a matrix of a repeated vector where each row is the same.
#' @param x vector to be repeated
#' @param n number of repeats
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
reprow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#' Hamardman product 
#'
#' hamardman.prod computes the Hamardman product of a vector of regression coefficients and a matrix of covariates.
#' @param coef vector of regression coefficients.
#' @param covar a matrix of covariates
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
hamardman.prod <- function(coef, covar){
  if(length(coef)!=ncol(covar)) stop("the number of coefficients is greater than the number of covariates")
  results <- reprow(coef,nrow(covar))*covar
  return(results)
}

#' IV-free exposure
#'
#' iv_free computes the IV-free exposure.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param g the instrumental variable.
#' @param q the number of quantiles the exposure distribution is to be split into. Within each quantile a causal effect will be fitted, known as a localised average causal effect (LACE). The default is deciles (i.e. 10 quantiles).
#' @param covar a matrix of covariates.
#' @param family a description of the error distribution and link function to be used in the model (either gaussian or binomial can be specified).
#' @return \item{xcoef}{the association between the exposure and the instrument.} 
#' @return \item{x0}{the IV-free exposure.}
#' @return \item{x0q}{the quantiles of x0.} 
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
iv_free <- function(y, x, g, covar, q, family="gaussian"){
  if(family=="gaussian"){
    if(!is.null(covar)){model <- lm(x~g+covar)}else{model <- lm(x~g)}
    x0 <- resid(model)
    xcoef <- model$coef[2]
  }
  if(family=="binomial"){
    if(!is.null(covar)){
      model <- lm(x[y==0]~g[y==0]+covar[y==0,])
      print(model$coef)
      x0 <- x - (model$coef[1] + model$coef[2]*g + rowSums(hamardman.prod(model$coef[3:length(model$coef)],covar)))
    }else{
      model <- lm(x[y==0]~g[y==0])
      x0 <- x - (model$coef[1] + model$coef[2]*g)
    }
    xcoef <- model$coef[2]
  }
  quantiles <- quantile(x0, probs=seq(0,1,1/q))
  x0q <- cut(x0, quantiles, include.lowest=T, labels=F)
  results <- list(xcoef=xcoef, x0=x0, x0q=x0q)
  return(results)
}

#' Localised avergae causal effects 
#'
#' lace computes the localised average causal effect for quantile.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param g the instrumental variable.
#' @param covar data.frame of covariates.
#' @param q the number of quantiles the exposure distribution is to be split into. Within each quantile a causal effect will be fitted, known as a localised average causal effect (LACE). The default is deciles (i.e. 10 quantiles).
#' @param x0q quantiles of x0 (the IV-free exposure).
#' @param xc_sub compute the association between the exposure and the insturment in each quantile of x0.
#' @param family a description of the error distribution and link function to be used in the model (either gaussian or binomial can be specified).
#' @param xpos the position used to relate x to the localised average causal effect. The default is the mean of the x-values within each quantile, otherwise specify a percentile (e.g. 0.5 corresponds to the median value).
#' @return List of results for the localised average causal effects for quantile.
#' @return \item{coef}{the localised average causal effect in each quantile.} 
#' @return \item{coef_se}{the IV-free exposure.}
#' @return \item{xmean}{the mean of the exposure in each quantile.} 
#' @return \item{xcoef_sub}{the association between the exposure and the instrument in each quantile.} 
#' @return \item{xcoef_sub_se}{the standard error of the association between the exposure and the instrument in each quantile}
#' @author James R Staley <js16174@bristol.ac.uk>
#' @export
lace <- function(y, x, g, covar=NULL, q, x0q, xc_sub=TRUE, family="gaussian", xpos="mean"){
  coef <- NULL
  coef_se <- NULL
  xmean<-NULL
  xcoef_sub <- NULL
  xcoef_sub_se <- NULL
  for(i in 1:q){
    if(family=="gaussian"){
      if(is.null(covar)){model <- lm(y[x0q==i]~g[x0q==i])}else{model<- lm(y[x0q==i]~g[x0q==i]+covar[x0q==i,])}
      coef[i] <- model$coef[2]  
      coef_se[i] <- summary(model)$coef[2,2]
    }
    if(family=="binomial"){
      if(is.null(covar)){model <- glm(y[x0q==i]~g[x0q==i],family="binomial")}else{model <- glm(y[x0q==i]~g[x0q==i]+covar[x0q==i,],family="binomial")}
      coef[i] <- model$coef[2]
      coef_se[i] <- summary(model)$coef[2,2]
    }
    if(xpos=="mean"){xmean[i] <- mean(x[x0q==i])}  
    if(xpos!="mean"){xmean[i] <- quantile(x[x0q==i], probs=xpos)}
    if(xc_sub){
      if(family=="gaussian"){
        if(is.null(covar)){xcoefs <- lm(x[x0q==i]~g[x0q==i])}else{xcoefs <- lm(x[x0q==i]~g[x0q==i]+covar[x0q==i,])}
        xcoef_sub[i] <- xcoefs$coef[2]
        xcoef_sub_se[i] <- summary(xcoefs)$coef[2,2]
      }
      if(family=="binomial"){
        if(is.null(covar)){xcoefs <- lm(x[(x0q==i & y==0)]~g[(x0q==i & y==0)])}else{xcoefs <- lm(x[(x0q==i & y==0)]~g[(x0q==i & y==0)]+covar[(x0q==i & y==0),])}
        xcoef_sub[i] <- xcoefs$coef[2]
        xcoef_sub_se[i] <- summary(xcoefs)$coef[2,2]
      }
    }
  }
  if(xc_sub){
    results <- data.frame(coef=coef, coef_se=coef_se, xmean=xmean, xcoef_sub=xcoef_sub, xcoef_sub_se=xcoef_sub_se)
  }else{
    results <- data.frame(coef=coef, coef_se=coef_se, xmean=xmean)
  }
  return(results)
}
