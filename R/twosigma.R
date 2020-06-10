##' Fit the TWO-SIGMA Model.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns correspoding to cells. It is recommended to make the rownames the gene names for better output.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or = 1 to indicate an intercept only model.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column), = 1 to indicate an intercept only model, or = 0 to indicate no zero-inflation model desired.
##' @param mean_re Should random intercepts be included in the (conditional) mean model? Ignored if adhoc=TRUE.
##' @param zi_re Should random intercepts be included in the zero-inflation model? Ignored if adhoc=TRUE.
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param adhoc Should the adhoc method be used by default to judge if random effects are needed?
##' @param adhoc_thresh Value below which the adhoc p-value is deemed significant (and thus RE are deemed necessary). Only used if adhoc==TRUE.
##' @param return_summary_fits If TRUE, the package returns a \code{summary.glmmTMB} object for each gene.  If FALSE, an object of class \code{glmmTMB} is returned for each gene. The latter requires far more memory to store.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model. Defaults to NULL for constant dispersion.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.  See \code{?glmmTMBControl}.
##' @param ncores Number of cores used for parallelization. Defaults to 1.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return If \code{return_summary_fits=TRUE}, returns a list of model fit objects of class \code{summary.glmmTMB}. If \code{return_summary_fits=FALSE}, returns a list of model fit objects of class \code{glmmTMB}. In either case, the order matches the row order of \code{count_matrix}, and the names of the list elements are taken as the rownames of \code{count_matrix}.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @import progress
##' @import doSNOW
##' @import snow
##' @import foreach
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm model.frame cor pnorm pt vcov
##' @importFrom multcomp glht mcp adjusted
##' @export twosigma

#zi_covar, mean_covar matrices must be specified? Need to figure out what to do exactly here
# should know what theyre doing for that though
twosigma<-function(count_matrix,mean_covar,zi_covar
  ,mean_re=TRUE,zi_re=TRUE
  ,id,adhoc=TRUE,adhoc_thresh=0.1
  ,return_summary_fits=TRUE,ncores=1
  ,disp_covar=NULL #need to be able to use data option?
  ,weights=rep(1,ncol(count_matrix))
  ,control = glmmTMBControl()){
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_covar","zi_covar","id")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}
  ngenes<-nrow(count_matrix)
  fit<-vector('list',length=ngenes)
  cl <- makeCluster(ncores)
  vars<-c("mean_covar","zi_covar")
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
    count<-count_matrix[i,,drop=FALSE]
    check_twosigma_input(count,mean_covar,zi_covar
      ,mean_re,zi_re
      ,disp_covar,adhoc=adhoc,id=id)
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
        message("adhoc method used to set both mean_re and zi_re to FALSE. Set adhoc=FALSE to customize mean_re and zi_re.")
      }
    }
    formulas<-create_model_formulas(mean_covar,zi_covar
      ,mean_form=NULL,zi_form=NULL
      ,mean_re,zi_re
      ,disp_covar)
    #browser()
    f<-glmmTMB(formula=formulas$mean_form
      ,ziformula=formulas$zi_form
      ,weights=weights
      ,dispformula = formulas$disp_form
      ,family=nbinom2,verbose = F
      ,control = control)

    if(return_summary_fits==TRUE){
      f<-summary(f)
    }
    return(f)
  }
  stopCluster(cl)
  for(i in 1:ngenes){
    tryCatch({
    fit[[i]]<-a[[i]]
    }
      ,error=function(e){})
  }
  names(fit)<-rownames(count_matrix)
  return(fit)
}


