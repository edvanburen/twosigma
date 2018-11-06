twosigma<-function(count,mean_covar=NULL,zi_covar=NULL
                   ,mean_form=NULL,zi_form=NULL,
                   mean_re=TRUE,zi_re=TRUE
                   ,disp_covar=NULL,data=NULL,id #need to be able to use data option
                   ,weights=rep(1,length(count))
                   ,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
                                                          ,step.max=.00001,step.min=.00001
                                                          ,rel.tol=1e-5,x.tol=1e-5))){
  #if(!grepl("nbinom2",family$family)){
  #  stop("Only the Negative Binomial Distribution is implemented in TWO-SIGMA")
  #}
  # Check that response is only valid counts if using Negative Binomial
    #if(grepl("nbinom2",family$family) & (sum(!as.matrix(count,ncol=1)%%1==0)>0 | min(count)<0)){

  check_twosigma_input(count,mean_covar,zi_covar
                       ,mean_form,zi_form
                       ,mean_re,zi_re
                       ,disp_covar)

  formulas<-create_model_formulas(mean_covar,zi_covar
    ,mean_form,zi_form
    ,mean_re,zi_re
    ,disp_covar)

  fit<-glmmTMB(formula=formulas$mean_form
            ,ziformula=formulas$zi_form
            ,weights=weights
            ,dispformula = formulas$disp_form
            ,family=nbinom2,verbose = F
            ,control = control)
return(fit)
}


