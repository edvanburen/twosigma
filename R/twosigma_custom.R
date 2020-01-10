##' Fit the TWO-component SInGle cell Model-based Association method of ... with custom user-specified model formulas.
##' @param count Vector of non-negative integer read counts. Users should ensure that this matches matches the LHS of the formula in "mean_form."
##' @param mean_form Custom two-sided model formula for the (conditional) mean model. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the dependent variable matches the argument to the parameter "count."
##' @param zi_form Custom one-sided model formula for the zero-inflation model. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param id Vector of individual-level ID's. Used for random effect prediction.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.
##' @section Details:
##' This function is likely only needed if users wish to include random effect terms beyond random intercepts. Users should be confident in their abilities to specify random effects using the syntax of lme4.
##' @return An object of class \code{glmmTMB}.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm
##' @export twosigma_custom
twosigma_custom<-function(count,mean_form,zi_form,id,
  disp_covar=NULL
  ,weights=rep(1,length(count))
  ,control = glmmTMBControl()){
  #if(!grepl("nbinom2",family$family)){
  #  stop("Only the Negative Binomial Distribution is implemented in TWO-SIGMA")
  #}
  # Check that response is only valid counts if using Negative Binomial
  #if(grepl("nbinom2",family$family) & (sum(!as.matrix(count,ncol=1)%%1==0)>0 | min(count)<0)){

  check_twosigma_custom_input(count,mean_form,zi_form,id,disp_covar)
  count<-as.numeric(count)
  if(is.null(disp_covar)){
    disp_form<- ~1 #Default is intercept only
  }else{
    if(is.atomic(disp_covar)&length(disp_covar)==1){
      if(disp_covar==0){stop("Invalid dispersion covariate option")}else{
        if(disp_covar==1){
          disp_form<-~1
        }else{
          stop("Invalid dispersion covariate option") #No zero-inflation component)
        }
      }
    }
  }

  if(is.null(disp_form)){
    disp_form<- ~ disp_covar}

  formulas<-list(mean_form=mean_form,zi_form=zi_form,disp_form=disp_form)

  fit<-glmmTMB(formula=formulas$mean_form
    ,ziformula=formulas$zi_form
    ,weights=weights
    ,dispformula = formulas$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)

  return(fit)
}


