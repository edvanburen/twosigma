##' Conveinent wrapper function for performing joint likelihood ratio tests using the TWO-SIGMA model.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns correspoding to cells. It is recommended to make the rownames the gene names for better output.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param covar_to_test Either a string indicating the column name of the covariate to test or an integer referring to its column position in BOTH the mean_covar and zi_covar matrices (if the two matrices differ using a string name is preferred). Argument is ignored if mean_covar and zi_covar are both a single covariate (that covariate is assumed of interest).
##' @param mean_re Should random intercepts be included in the (conditional) mean model?
##' @param zi_re Should random intercepts be included in the zero-inflation model?
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param return_full_fits If TRUE, fit objects of class glmmTMB are returned. If FALSE, only objects of class summary.glmmTMB are returned. The latter require a much larger amount of memory to store.
##' @param adhoc Should the adhoc method be used by default to judge if random effects are needed?
##' @param adhoc_thresh Value below which the adhoc p-value is deemed significant (and thus RE are deemed necessary). Only used if adhoc==TRUE.
##' @param silent If TRUE, progress is not printed.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.  See \code{?glmmTMBControl}.
##' @param control Control parameters for optimization in glmmTMB.
##' @param ncores Number of cores used for parallelization. Defaults to 1.
##' @section Details:
##' This function assumes that the variable being tested is in both components of the model (and thus that the zero-inflation component exists and contains more than an Intercept). Users wishing to do fixed effect testing in other cases or specify custom model formulas they will need to construct the statistics themselves using either two separate calls to \code{twosigma} or the \code{lr.twosigma_custom} fumction. If \code{adhoc=TRUE}, any input in mean_re and zi_re will be ignored. If either model fails to converge, or the LR statistic is negative, both the statistic and p-value are assigned as NA.
##'
##' @return A list containing model fit objects of class \code{glmmTMB} (only if \code{return_full_fits}=TRUE), model fit summary objects of class \code{summary.glmmTMB}, the 2 d.f. Likelihood Ratio statistic, and the p-value.
##' @export lr.twosigma

lr.twosigma<-function(count_matrix,mean_covar,zi_covar,covar_to_test
                      ,mean_re=FALSE,zi_re=FALSE,
                       id,return_full_fits=TRUE,adhoc=FALSE,adhoc_thresh=0.1
                      ,silent=FALSE,ncores=1
                      ,disp_covar=NULL
                      ,weights=rep(1,ncol(count_matrix))
                      ,control=glmmTMBControl())
                      #,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
                       # ,step.max=.00001,step.min=.00001
                      #  ,rel.tol=1e-5,x.tol=1e-5)))
  {
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_covar","zi_covar","covar_to_test","id")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}
  ngenes<-nrow(count_matrix)
  LR_stat<-rep(NA,length=ngenes)
  p.val<-rep(NA,length=ngenes)
  sum_fit_alt<-vector('list',length=ngenes)
  sum_fit_null<-vector('list',length=ngenes)

  if(return_full_fits==TRUE){
    fits_all_null<-vector('list',length=ngenes)
    fits_all_alt<-vector('list',length=ngenes)
  }
  mc<-mean_covar
  zc<-zi_covar
  cl <- makeCluster(ncores)
  vars<-c("mean_covar","zi_covar","mc","zc")
  clusterExport(cl,list=vars,envir = environment())
  registerDoSNOW(cl)
  pb <- progress_bar$new(
    format = "num genes complete = :num [:bar] :elapsed | eta: :eta",
    total = ngenes,    # 100
    width = 60)

  progress <- function(n){
    pb$tick(tokens = list(num = n))
  }
  opts <- list(progress = progress)
  print("Running Gene-Level Models")
  #browser()
  a<-foreach(i=1:ngenes,.options.snow = opts)%dopar%{
    mean_covar<-mc
    zi_covar<-zc
    count<-count_matrix[i,,drop=FALSE]
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
        print("adhoc method used to set both mean_re and zi_re to FALSE. Set adhoc=FALSE to customize user-inputted values for mean_re and zi_re.")
      }
    }
    formulas<-create_model_formulas(mean_covar,zi_covar
      ,mean_form=NULL,zi_form=NULL
      ,mean_re,zi_re
      ,disp_covar)
    if(is.atomic(zi_covar)&length(zi_covar)==1){
      if(zi_covar==1 | zi_covar==0){
        stop("This function is meant for joint testing of a covariate when it is present in both components. Please see funtion lr.twosigma_custom for other cases.")
      }
    }
    if(is.atomic(mean_covar)&length(mean_covar)==1){
      if(mean_covar==1){
        stop("This function is meant for joint testing of a covariate when it is present in both components. Please see funtion lr.twosigma_custom for other cases.")
      }
    }
    if(is.vector(mean_covar)&is.vector(zi_covar)){
      if(!identical(mean_covar,zi_covar)){
        stop("Mean covariate is not idential to zi covariate. This function is only designed for cases in which the same covariate is being designed in both components.")
      }
    }
    if(is.character(covar_to_test)){
      #if(!is.vector(mean_covar)&!is.vector(zi_covar))
      if(is.matrix(mean_covar)&is.matrix(zi_covar))
      {
        if(!(covar_to_test%in%colnames(mean_covar) & covar_to_test%in%colnames(zi_covar))){
          stop("covar_to_test not found in both matrices")
        }
      }
    }
    fit_alt<-glmmTMB(formula=formulas$mean_form
      ,ziformula=formulas$zi_form
      ,weights=weights
      ,dispformula = formulas$disp_form
      ,family=nbinom2,verbose = F
      ,control = control)
    rm(formulas)
    #If numeric we are assuming that the variable is in the same position
    # need also to point out that users have some responsibilities here
    if(is.numeric(covar_to_test)){
      if(is.matrix(mean_covar)&is.matrix(zi_covar)){
        if(covar_to_test>max(ncol(mean_covar),ncol(zi_covar))){
          stop("covar_to_test seems to be ill-defined")
        }
      }
      if(is.vector(mean_covar)& is.vector(zi_covar)){
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
        mean_covar<-mean_covar[,-covar_to_test,drop=FALSE]
        if(ncol(mean_covar)==0){mean_covar<-1}
        zi_covar<-1
        formulas<-create_model_formulas(mean_covar,zi_covar
          ,mean_form=NULL,zi_form=NULL
          ,mean_re,zi_re
          ,disp_covar)
        #formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-covar_to_test]",format(formulas$mean_form)))
        #formulas$zi_form<-as.formula(gsub("zi_covar","1",format(formulas$zi_form)))
      }
      if(is.vector(mean_covar)& is.matrix(zi_covar)){
        mean_covar<-1
        zi_covar<-zi_covar[,-covar_to_test,drop=FALSE]
        if(ncol(zi_covar)==0){zi_covar<-1}
        formulas<-create_model_formulas(mean_covar,zi_covar
          ,mean_form=NULL,zi_form=NULL
          ,mean_re,zi_re
          ,disp_covar)
        #formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
        #formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-covar_to_test]",format(formulas$zi_form)))
      }
      if(is.matrix(mean_covar) & is.matrix(zi_covar)){
        mean_covar<-mean_covar[,-covar_to_test,drop=FALSE]
        zi_covar<-zi_covar[,-covar_to_test,drop=FALSE]
        if(ncol(mean_covar)==0){mean_covar<-1}
        if(ncol(zi_covar)==0){zi_covar<-1}
        formulas<-create_model_formulas(mean_covar,zi_covar
          ,mean_form=NULL,zi_form=NULL
          ,mean_re,zi_re
          ,disp_covar)
        #formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-covar_to_test]",format(formulas$mean_form)))
        #formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-covar_to_test]",format(formulas$zi_form)))
      }
    }else{
      if(is.matrix(mean_covar)){
        index<-which(colnames(mean_covar)==covar_to_test)
        if(length(index)==0){
          stop("covar_to_test Name not found in colnames of mean model covariate data matrix")
        }
        #formulas$mean_form<-as.formula(gsub("mean_covar","mean_covar[,-index]",format(formulas$mean_form)))
        mean_covar<-mean_covar[,-index,drop=FALSE]
        if(ncol(mean_covar)==0){mean_covar<-1}
      }
      if(is.vector(mean_covar)){
        #formulas$mean_form<-as.formula(gsub("mean_covar","1",format(formulas$mean_form)))
        mean_covar<-1
      }
      if(is.matrix(zi_covar)){
        index<-which(colnames(zi_covar)==covar_to_test)
        if(is.null(index)){
          stop("covar_to_test Name not found in colnames of zi model covariate data matrix")
        }
        #formulas$zi_form<-as.formula(gsub("zi_covar","zi_covar[,-index]",format(formulas$zi_form)))
        zi_covar<-zi_covar[,-index,drop=FALSE]
        if(ncol(zi_covar)==0){zi_covar<-1}
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
    tryCatch({
      sum_null<-summary(fit_null)
      sum_alt<-summary(fit_alt)},error=function(e){})
    LR_stat<-NA
    tryCatch({
      LR_stat<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
      if(LR_stat<0 | (!fit_alt$sdr$pdHess) | (!fit_null$sdr$pdHess)){
        LR_stat<-NA}
      p.val<-1-pchisq(LR_stat,df=2)},error=function(e){})

    if(return_full_fits==TRUE){
      return(list(sum_null=sum_null,sum_alt=sum_alt
        ,fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,p.val=p.val))
    }else{
      return(list(sum_null=sum_null,sum_alt=sum_alt,LR_stat=LR_stat,p.val=p.val))
    }

  }
  stopCluster(cl)
    #if(is.character(covar_to_test)){}
    #names_null<-rownames(sum_null$coefficients$cond)
    #names_alt<-rownames(sum_alt$coefficients$cond)

    #est<-sum_alt$coefficients$cond[which(!names_alt%in%names_null),1]
    #est_zi<-sum_alt$coefficients$zi[which(!names_alt%in%names_null),1]
#browser()
for(i in 1:ngenes){
  tryCatch({
    p.val[i]<-a[[i]]$p.val
    LR_stat[i]<-a[[i]]$LR_stat
    #browser()
    sum_fit_alt[[i]]<-a[[i]]$sum_alt
    sum_fit_null[[i]]<-a[[i]]$sum_null
    if(return_full_fits==TRUE){
      fits_all_null[[i]]<-a[[i]]$fit_null
      fits_all_alt[[i]]<-a[[i]]$fit_alt
    }
  },error=function(e){})
    #if(!silent){print(paste("Finished Gene Number",i,"of",ngenes))}
  }
  names(p.val)<-rownames(count_matrix)
  names(LR_stat)<-rownames(count_matrix)
  names(sum_fit_alt)<-rownames(count_matrix)
  names(sum_fit_null)<-rownames(count_matrix)

  if(return_full_fits==TRUE){
    names(fits_all_null)<-rownames(count_matrix)
    names(fits_all_alt)<-rownames(count_matrix)
    return(list(full_fit_null=fits_all_null,full_fit_alt=fits_all_alt,summary_fit_null=sum_fit_null,summary_fit_alt=sum_fit_alt,LR_stat=LR_stat,LR_p.val=p.val))
  }else{
    return(list(summary_fit_null=sum_fit_null,summary_fit_alt=sum_fit_alt,LR_stat=LR_stat,LR_p.val=p.val))
  }
}
