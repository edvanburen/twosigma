##' Conveinent wrapper function for performing joint likelihood ratio tests with the TWO-SIGMA model of ... with custom user-specified formulas.
##' @param count Vector of non-negative integer read counts.
##' @param mean_form_alt Custom two-sided model formula for the (conditional) mean model under the null. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the dependent variable matches the argument to the parameter "count."
##' @param zi_form_alt Custom one-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param mean_form_null Custom two-sided model formula for the (conditional) mean model under the null. Syntax is as in \code{mean_form_alt}.
##' @param zi_form_null Custom one-sided model formula for the zero-inflation model under the null. Syntax is as in \code{zi_form_alt}.
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param lr.df Degrees of Freedom for the constructed likelihood ratio test. Must be a non-negative integer.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in glmmTMB.
##' @section Details:
##' This function is a wrapper for conducting fixed effect likelihood ratio tests with twosigma.  There is no checking to make sure that the alt and null model formulas represent a valid likelihood ratio test when fit together.  Users must ensure that inputted formulas represent valid nested models. If either model fails to converge, or the LR statistic is negative, both the statistic and p-value are assigned as NA.
##'
##' @return A list containing glmmTMB objects of model fits under the null and alternative, the Likelihood Ratio statistic, associated p-value, and all model formulas used.
##' @export lr.twosigma_custom
# Likely will want to remove the ability to input a formula for this fn to work properly

lr.twosigma_custom<-function(count,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null
                      ,id,lr.df,
                      covar_logFC=NULL,
                      disp_covar=NULL
                      ,weights=rep(1,length(count))
                      ,control=glmmTMBControl())
                      #,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
                      #  ,step.max=.00001,step.min=.00001
                      #  ,rel.tol=1e-5,x.tol=1e-5)))
  {
  check_twosigma_custom_input(count,mean_form_alt,zi_form_alt,id,disp_covar)
  check_twosigma_custom_input(count,mean_form_null,zi_form_null,id,disp_covar)
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

sum_null<-summary(fit_null)
sum_alt<-summary(fit_alt)
#browser()
#if(is.character(contrast)){}
names_null<-rownames(sum_null$coefficients$cond)
names_alt<-rownames(sum_alt$coefficients$cond)

index<-which(grepl(covar_logFC,names_alt))
est<-sum_alt$coefficients$cond[index,1]

names_alt<-rownames(sum_alt$coefficients$cond)
index<-which(grepl(covar_logFC,names_alt))
est_zi<-sum_alt$coefficients$zi[index,1]

LR_stat<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
if(LR_stat<0 | (!fit_alt$sdr$pdHess) | (!fit_null$sdr$pdHess)){
  LR_stat<-NA
  message("LR stat set to NA, indicative of model specification or fitting problem")}
p.val<-1-pchisq(LR_stat,df=lr.df)
return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,LR_p.val=p.val,mean_comp_logFC=est,zi_comp_est=est_zi,mean_form_alt=mean_form_alt,zi_form_alt=zi_form_alt,mean_form_null=mean_form_null,zi_form_null=zi_form_null))
}
