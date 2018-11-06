# Likely will want to remove the ability to input a formula for this fn to work properly
lr.twosigma<-function(contrast,joint=TRUE,count,mean_covar=NULL,zi_covar=NULL
                      ,mean_form=NULL,zi_form=NULL
                      ,mean_re=TRUE,zi_re=TRUE
                      ,disp_covar=NULL,data=NULL,id
                      ,weights=rep(1,length(count))
                      ,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
                        ,step.max=.00001,step.min=.00001
                        ,rel.tol=1e-5,x.tol=1e-5))){

  check_twosigma_input(count,mean_covar,zi_covar
    ,mean_form,zi_form
    ,mean_re,zi_re
    ,disp_covar)

  formulas<-create_model_formulas(mean_covar,zi_covar
    ,mean_form,zi_form
    ,mean_re,zi_re
    ,disp_covar)

  fit_alt<-glmmTMB(formula=formulas$mean_form
    ,ziformula=formulas$zi_form
    ,weights=weights
    ,dispformula = formulas$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)

  #If numeric we are assuming that the variable is in the same position
  if(is.numeric(contrast)){
    formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-contrast]",format(formulas$mean_form)))
    formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-contrast]",format(formulas$zi_form)))
  }else{
    index<-which(colnames(mean_covar)==contrast)
    formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-index]",format(formulas$mean_form)))
    formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-index]",format(formulas$zi_form)))
  }
fit_null<-glmmTMB(formula=formulas$mean_form
    ,ziformula=formulas$zi_form
    ,weights=weights
    ,dispformula = formulas$disp_form
    ,family=nbinom2,verbose = F
    ,control = control)
LR_stat<-2*(fit_null$fit$objective-fit_alt$fit$objective)
return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat))
}
