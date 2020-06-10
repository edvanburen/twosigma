##' Conveinent wrapper function for performing joint likelihood ratio tests with the TWO-SIGMA model using custom user-specified formulas.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns correspoding to cells. It is recommended to make the rownames the gene names for better output.
##' @param mean_form_alt Custom two-sided model formula for the (conditional) mean model under the null. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the dependent variable matches the argument to the parameter "count."
##' @param zi_form_alt Custom one-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param mean_form_null Custom two-sided model formula for the (conditional) mean model under the null. Syntax is as in \code{mean_form_alt}.
##' @param zi_form_null Custom one-sided model formula for the zero-inflation model under the null. Syntax is as in \code{zi_form_alt}.
##' @param id Vector of individual-level (sample-level) ID's. Used for random effect prediction but required regardless of their presence in the model.
##' @param lr.df Degrees of Freedom for the constructed likelihood ratio test. Must be a non-negative integer.
##' @param return_full_fits If TRUE, full fit objects of class glmmTMB are returned.  If FALSE, only fit objects of class summary.glmmTMB are returned.  The latter requires far less memory to store.
##' @param silent If TRUE, progress is not printed.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in glmmTMB.  See \code{?glmmTMBControl}.
##' @param ncores Number of cores used for parallelization. Defaults to 1.
##' @section Details:
##' This function is a wrapper for conducting fixed effect likelihood ratio tests with twosigma.  There is no checking to make sure that the alt and null model formulas represent a valid likelihood ratio test when fit together.  Users must ensure that inputted formulas represent valid nested models. If either model fails to converge, or the LR statistic is negative, both the statistic and p-value are assigned as NA.
##'
##' @return A list containing model fit objects of class \code{glmmTMB} (only if \code{return_full_fits}=TRUE), model fit summary objects of class \code{summary.glmmTMB}, the 2 d.f. Likelihood Ratio statistic, and the p-value, and all model formulas used.
##' @export lr.twosigma_custom

lr.twosigma_custom<-function(count_matrix,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null
                      ,id,lr.df,return_full_fits=TRUE,ncores=1,
                      disp_covar=NULL
                      ,weights=rep(1,ncol(count_matrix))
                      ,control=glmmTMBControl(),silent=FALSE)
  {
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_form_alt","zi_form_alt","mean_form_null","zi_form_null","id","lr.df")
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

  cl <- makeCluster(ncores)
  vars<-unique(c(all.vars(mean_form_alt)[-1],all.vars(zi_form_alt)
    ,all.vars(mean_form_null)[-1],all.vars(zi_form_null)))
  vars<-vars[!vars=="id"]
  clusterExport(cl,list=vars)
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
  a<-foreach(i=1:ngenes,.options.snow = opts)%dopar%{
    count<-count_matrix[i,,drop=FALSE]
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
    #if(!silent){print(paste("Finished Gene Number",i,"of",ngenes))}
  },error=function(e){})
  }
  names(p.val)<-rownames(count_matrix)
  names(LR_stat)<-rownames(count_matrix)
  names(sum_fit_alt)<-rownames(count_matrix)
  names(sum_fit_null)<-rownames(count_matrix)
  if(return_full_fits==TRUE){
    names(fits_all_null)<-rownames(count_matrix)
    names(fits_all_alt)<-rownames(count_matrix)
    return(list(full_fit_null=fits_all_null,full_fit_alt=fits_all_alt,summary_fit_null=sum_fit_null,summary_fit_alt=sum_fit_alt,LR_stat=LR_stat,LR_p.val=p.val,mean_form_alt=mean_form_alt,zi_form_alt=zi_form_alt,mean_form_null=mean_form_null,zi_form_null=zi_form_null))
  }else{
  return(list(summary_fit_null=sum_fit_null,summary_fit_alt=sum_fit_alt,LR_stat=LR_stat,LR_p.val=p.val,mean_form_alt=mean_form_alt,zi_form_alt=zi_form_alt,mean_form_null=mean_form_null,zi_form_null=zi_form_null))
}
}
