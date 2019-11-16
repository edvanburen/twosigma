##' lr.twosigma: Conveinent wrapper function for performing joint likelihood ratio tests in TWO-SIGMA.
##' @param count Vector of non-negative integer read counts.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in glmmTMB.
##' @section Details:
##' This function is a wrapper for conducting fixed effect likelihood ratio tests with twosigma.  There is no (absolutely no) checking to make sure that the alt and null model formulas represent a valid likelihood ratio test when fit together.  Users must ensure that inputted formulas represent valid nested models.
##'
##' @return A list containing glmmTMB objects of model fits under the null and alternative, the 2 d.f. Likelihood Ratio statistic, and the p-value.
##' @export lr.twosigma_custom
# Likely will want to remove the ability to input a formula for this fn to work properly

lr.twosigma_custom<-function(count,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null
                      ,disp_covar=NULL,id,adhoc=FALSE
                      ,weights=rep(1,length(count))
                      ,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
                        ,step.max=.00001,step.min=.00001
                        ,rel.tol=1e-5,x.tol=1e-5))){
  check_twosigma_custom_input(count,mean_form_alt,zi_form_alt,id,disp_covar)
  check_twosigma_custom_input(count,mean_form_null,zi_form_null,id,disp_covar)

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

  formulas_alt<-list(mean_form=mean_form_alt,zi_form=zi_form_alt,disp_form=disp_form)
  formulas_null<-list(mean_form=mean_form_null,zi_form=zi_form_null,disp_form=disp_form)

  fit_alt<-glmmTMB(formula=formulas_alt$mean_form
    ,ziformula=formulas_alt$zi_form
    ,weights=weights
    ,dispformula = formulas_alt$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)

fit_null<-glmmTMB(formula=formulas_null$mean_form
    ,ziformula=formulas_null$zi_form
    ,weights=weights
    ,dispformula = formulas_null$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)

LR_stat<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
p.val<-1-pchisq(LR_stat,df=2)
return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,p.val=p.val))
}
