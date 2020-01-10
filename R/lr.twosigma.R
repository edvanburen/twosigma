##' Conveinent wrapper function for performing joint likelihood ratio tests using the TWO-SIGMA model of ...
##' @param count Vector of non-negative integer read counts.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param contrast Either a string indicating the column name of the covariate to test or an integer referring to its column position in BOTH the mean_covar and zi_covar matrices (if the two matrices differ using a string name is preferred). Argument is ignored if mean_covar and zi_covar are both a single covariate (that covariate is assumed of interest).
##' @param mean_re Should random intercepts be included in the (conditional) mean model?
##' @param zi_re Should random intercepts be included in the zero-inflation model?
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param adhoc Should the adhoc method be used by default to judge if random effects are needed?
##' @param adhoc_thresh Value below which the adhoc p-value is deemed significant (and thus RE are deemed necessary). Only used if adhoc==TRUE.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in glmmTMB.
##' @section Details:
##' This function assumes that the variable being tested is in both components of the model (and thus that the zero-inflation component exists and contains more than an Intercept). Users wishing to do fixed effect testing in other cases or specify custom model formulas they will need to construct the statistics themselves using either two separate calls to \code{twosigma} or the \code{lr.twosigma_custom} fumction. If adhoc==TRUE, any input in mean_re and zi_re will be ignored. If either model fails to converge, or the LR statistic is negative, both the statistic and p-value are assigned as NA.
##'
##' @return A list containing glmmTMB objects of model fits under the null and alternative, the 2 d.f. Likelihood Ratio statistic, and the p-value.
##' @export lr.twosigma

lr.twosigma<-function(count,mean_covar,zi_covar,contrast
                      ,mean_re=FALSE,zi_re=FALSE,
                       id,adhoc=FALSE,adhoc_thresh=0.1
                      ,disp_covar=NULL
                      ,weights=rep(1,length(count))
                      ,control=glmmTMBControl())
                      #,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
                       # ,step.max=.00001,step.min=.00001
                      #  ,rel.tol=1e-5,x.tol=1e-5)))
  {
  check_twosigma_input(count,mean_covar,zi_covar
    ,mean_re,zi_re
    ,disp_covar,id=id,adhoc=adhoc)
  count<-as.numeric(count)
  if(adhoc==TRUE){
    if(is.atomic(zi_covar)&length(zi_covar)==1){
      if(zi_covar==0){stop("adhoc method only implemented when ZI model contains at minimum an intercept. Please either set adhoc=FALSE or specify at minimum an intercept in the ZI model.")}}
    p.val<-adhoc.twosigma(count=count,mean_covar=mean_covar,zi_covar = zi_covar,id=id,weights=weights)
    if(p.val<adhoc_thresh){
      mean_re=TRUE
      zi_re=TRUE
      message("adhoc method used to set both mean_re and zi_re to TRUE. Set adhoc=FALSE to customize mean_re and zi_re.")
    }else{
      mean_re=FALSE
      zi_re=FALSE
      message("adhoc method used to set both mean_re and zi_re to FALSE. Set adhoc=FALSE to customize user-inputted values for mean_re and zi_re.")
    }
  }

  formulas<-create_model_formulas(mean_covar,zi_covar
    ,mean_form=NULL,zi_form=NULL
    ,mean_re,zi_re
    ,disp_covar)

  if(is.atomic(zi_covar)&length(zi_covar)==1){
    if(zi_covar==1 | zi_covar==0){
      stop("This function is meant for joint testing of a covariate when it is present in both components.")
    }
  }
  if(is.atomic(mean_covar)&length(mean_covar)==1){
    if(mean_covar==1){
      stop("This function is meant for joint testing of a covariate when it is present in both components.")
    }
  }
  if(is.vector(mean_covar)&is.vector(zi_covar)){
    if(!identical(mean_covar,zi_covar)){
      stop("Mean covariate is not idential to zi covariate. This function is only designed for cases in which the same covariate is being designed in both components.")
    }
  }
  if(is.character(contrast)){
      #if(!is.vector(mean_covar)&!is.vector(zi_covar))
      if(is.matrix(mean_covar)&is.matrix(zi_covar))
      {
        if(!(contrast%in%colnames(mean_covar) & contrast%in%colnames(zi_covar))){
          stop("contrast not found in both matrices")
      }
    }
  }
  fit_alt<-glmmTMB(formula=formulas$mean_form
    ,ziformula=formulas$zi_form
    ,weights=weights
    ,dispformula = formulas$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)
  #If numeric we are assuming that the variable is in the same position
  # need also to point out that users have some responsibilities here

  if(is.numeric(contrast)){
    if(is.matrix(mean_covar)&is.matrix(zi_covar)){
      if(contrast>max(ncol(mean_covar),ncol(zi_covar))){
        stop("Contrast seems to be ill-defined")
      }
    }
    if(is.vector(mean_covar) & is.vector(zi_covar)){
      mean_covar<-1
      zi_covar<-1
      formulas<-create_model_formulas(mean_covar,zi_covar
        ,mean_form=NULL,zi_form=NULL
        ,mean_re,zi_re
        ,disp_covar)
      #formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
      #formulas$zi_form<-as.formula(gsub("zi_covar","1",format(formulas$zi_form)))
    }
    if(is.matrix(mean_covar) & is.vector(zi_covar)){
      mean_covar<-mean_covar[,-contrast]
      zi_covar<-1
      formulas<-create_model_formulas(mean_covar,zi_covar
        ,mean_form=NULL,zi_form=NULL
        ,mean_re,zi_re
        ,disp_covar)
      #formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-contrast]",format(formulas$mean_form)))
      #formulas$zi_form<-as.formula(gsub("zi_covar","1",format(formulas$zi_form)))
    }
    if(is.vector(mean_covar) & is.matrix(zi_covar)){
      mean_covar<-1
      zi_covar<-zi_covar[,-contrast]
      formulas<-create_model_formulas(mean_covar,zi_covar
        ,mean_form=NULL,zi_form=NULL
        ,mean_re,zi_re
        ,disp_covar)
      #formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
      #formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-contrast]",format(formulas$zi_form)))
    }
    if(is.matrix(mean_covar) & is.matrix(zi_covar)){
      mean_covar<-mean_covar[,-contrast]
      zi_covar<-zi_covar[,-contrast]
      formulas<-create_model_formulas(mean_covar,zi_covar
        ,mean_form=NULL,zi_form=NULL
        ,mean_re,zi_re
        ,disp_covar)
      #formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-contrast]",format(formulas$mean_form)))
      #formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-contrast]",format(formulas$zi_form)))
    }
  }else{
    if(is.matrix(mean_covar)){
      index<-which(colnames(mean_covar)==contrast)
      if(length(index)==0){
        stop("Contrast Name not found in colnames of mean model covariate data matrix")
      }
      #formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-index]",format(formulas$mean_form)))
      mean_covar<-mean_covar[,-index]
    }
    if(is.vector(mean_covar)){
      #formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
      mean_covar<-1
    }
    if(is.matrix(zi_covar)){
    index<-which(colnames(zi_covar)==contrast)
    if(is.null(index)){
      stop("Contrast Name not found in colnames of zi model covariate data matrix")
    }
    #formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-index]",format(formulas$zi_form)))
    zi_covar<-zi_covar[,-index]
    }
    if(is.vector(zi_covar)){
      #formulas$zi_form<-as.formula(gsub("zi_covar","1",format(formulas$zi_form)))
      zi_covar<-1
    }
    formulas<-create_model_formulas(mean_covar,zi_covar
      ,mean_form=NULL,zi_form=NULL
      ,mean_re,zi_re
      ,disp_covar)
  }

fit_null<-glmmTMB(formula=formulas$mean_form
    ,ziformula=formulas$zi_form
    ,weights=weights
    ,dispformula = formulas$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)
LR_stat<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
if(LR_stat<0 | (!fit_alt$sdr$pdHess) | (!fit_null$sdr$pdHess)){
  LR_stat<-NA
  message("LR stat set to NA, indicative of model specification or fitting problem")}
p.val<-1-pchisq(LR_stat,df=2)
return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,LR_p.val=p.val))
}
