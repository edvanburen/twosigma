##' lr.twosigma: Conveinent wrapper function for performing joint likelihood ratio tests in TWO-SIGMA
##' @param count Vector of non-negative integer read counts.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param contrast Either a string indicating the column name of the covariate to test or an integer referring to its column position in BOTH the mean_covar and zi_covar matrices (if the two matrices differ using a string name is preferred). Argument is ignored if mean_covar and zi_covar are both a single covariate (as that covariate is assumed of interest).
##' @param mean_re Should random intercepts be included in the (conditional) mean model?
##' @param zi_re Should random intercepts be included in the zero-inflation model?
##' @param id Vector of individual-level ID's. Used for random effect prediction.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in glmmTMB.
##' @section Details:
##' This function assumes that the variable being tested is in both components of the model (and thus that the zero-inflation component exists and contains more than an Intercept). Users wishing to do fixed effect testing in other cases will need to construct the statistics themselves using two separate calls to \code{twosigma}. If users wish to specify custom model formulas they should also construct the likelihood ratio statistics themselves to ensure the test is being conducted properly.
##'
##' @return A list containing glmmTMB objects of model fits under the null and alternative, the 2 d.f. Likelihood Ratio statistic, and the p-value.
##' @export lr.twosigma
# Likely will want to remove the ability to input a formula for this fn to work properly

lr.twosigma<-function(count,mean_covar,zi_covar,contrast#,joint=TRUE
                      #,mean_form=NULL,zi_form=NULL
                      ,mean_re=TRUE,zi_re=TRUE
                      ,disp_covar=NULL,id
                      ,weights=rep(1,length(count))
                      ,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
                        ,step.max=.00001,step.min=.00001
                        ,rel.tol=1e-5,x.tol=1e-5))){
  check_twosigma_input(count,mean_covar,zi_covar
    ,mean_re,zi_re
    ,disp_covar)

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

  if(is.character(contrast)){
      if(!is.atomic(mean_covar)&!is.atomic(zi_covar))
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
    if(!is.atomic(mean_covar)&!is.atomic(zi_covar)){
      if(contrast>max(ncol(mean_covar),ncol(zi_covar))){
        stop("Contrast seems to be ill-defined")
      }
    }
    if(is.atomic(mean_covar) & is.atomic(zi_covar)){
      formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
      formulas$zi_form<-as.formula(gsub("zi_covar","1",format(formulas$zi_form)))
    }
    if(!is.atomic(mean_covar) & is.atomic(zi_covar)){
      formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-contrast]",format(formulas$mean_form)))
      formulas$zi_form<-as.formula(gsub("zi_covar","1",format(formulas$zi_form)))
    }
    if(is.atomic(mean_covar) & !is.atomic(zi_covar)){
      formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
      formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-contrast]",format(formulas$zi_form)))
    }
    if(!is.atomic(mean_covar) & !is.atomic(zi_covar)){
      formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-contrast]",format(formulas$mean_form)))
      formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-contrast]",format(formulas$zi_form)))
    }
  }else{
    if(!is.atomic(mean_covar)){
      index<-which(colnames(mean_covar)==contrast)
      if(length(index)==0){
        stop("Contrast Name not found in colnames of mean model covariate data matrix")
      }
      formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-index]",format(formulas$mean_form)))
    }
    if(is.atomic(mean_covar)){
      formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
    }
    if(!is.atomic(zi_covar)){
    index<-which(colnames(zi_covar)==contrast)
    if(is.null(index)){
      stop("Contrast Name not found in colnames of zi model covariate data matrix")
    }
    formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-index]",format(formulas$zi_form)))
    }
    if(is.atomic(zi_covar)){
      formulas$zi_form<-as.formula(gsub("zi_covar","1",format(formulas$zi_form)))
    }
  }

fit_null<-glmmTMB(formula=formulas$mean_form
    ,ziformula=formulas$zi_form
    ,weights=weights
    ,dispformula = formulas$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)
LR_stat<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
p.val<-1-pchisq(LR_stat,df=2)
return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,p.val=p.val))
}
