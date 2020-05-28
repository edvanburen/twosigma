##' Fit the TWO-SIGMA model with custom user-specified model formulas.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns correspoding to cells. It is recommended to make the rownames the gene names for better output.
##' @param mean_form Custom two-sided model formula for the (conditional) mean model. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the LHS of the formula begins with "count."
##' @param zi_form Custom one-sided model formula for the zero-inflation model. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param id Vector of individual-level (sample-level) ID's. Used for random effect prediction but required regardless of their presence in the model.
##' @param return_summary_fits If TRUE, the package returns a summary.glmmTMB object for each gene.  If TRUE, a glmmTMB object is returned for each gene. The latter requires far more storage space.
##' @param silent If TRUE, progress is not printed.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.  See \code{?glmmTMBControl}.
##' @section Details:
##' This function is likely only needed if users wish to include random effect terms beyond random intercepts. Users should be confident in their abilities to specify random effects using the syntax of lme4.
##' @return If \code{return_summary_fits=TRUE}, returns a list of objects of class \code{summary.glmmTMB}. If \code{return_summary_fits=FALSE}, returns a list of model fit objects of class \code{glmmTMB}. In either case, the order matches the row order of \code{count_matrix}, and the names of the list elements are taken as the rownames of \code{count_matrix}.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm
##' @export twosigma_custom
twosigma_custom<-function(count_matrix,mean_form,zi_form,id,return_summary_fits=TRUE,
  silent=FALSE,disp_covar=NULL
  ,weights=rep(1,ncol(count_matrix))
  ,control = glmmTMBControl()){

  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_form","zi_form","id")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}

  ngenes<-nrow(count_matrix)
  fit<-vector('list',length=ngenes)
  #browser()
  for(i in 1:ngenes){
    count<-count_matrix[i,,drop=FALSE]
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

    f<-glmmTMB(formula=formulas$mean_form
      ,ziformula=formulas$zi_form
      ,weights=weights
      ,dispformula = formulas$disp_form
      ,family=nbinom2,verbose = F
      ,control = control)
    if(return_summary_fits){
      fit[[i]]<-summary(f)
    }else{
      fit[[i]]<-f
    }
    if(!silent){print(paste("Finished Gene Number",i,"of",ngenes))}
  }
  names(fit)<-rownames(count_matrix)
  return(fit)
}


