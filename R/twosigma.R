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
##' @param ncores Number of cores used for parallelization. Defaults to 1, meaning no parallelization of any kind is done.
##' @param cluster_type Whether to use a "cluster of type "Fork" or "Sock". On Unix systems, "Fork" will likely improve performance. On Windows, only "Sock" will actually result in parallelized computing.
##' @param chunk_size Number of genes to be sent to each parallel environment. Parallelization is more efficient, particuarly with a large count matrix, when the count matrix is 'chunked' into some common size (e.g. 10, 50, 200). Defaults to 10.
##' @param lb Should load balancing be used for parallelization? Users will likely want to set to FALSE for improved performance.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return A list with the following elements:
##' ##' \itemize{
##' \item{\code{fit: }} If \code{return_summary_fits=TRUE}, returns a list of model fit objects of class \code{summary.glmmTMB}. If \code{return_summary_fits=FALSE}, returns a list of model fit objects of class \code{glmmTMB}. In either case, the order matches the row order of \code{count_matrix}, and the names of the list elements are taken as the rownames of \code{count_matrix}.
##' \item{\code{adhoc_include_RE: }} Logical vector indicator whether the adhoc method determined random effects needed.  If \code{adhoc=F}, then a vector of NA's.
##' ##' }
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @import parallel
##' @import pbapply
##' @import doParallel
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm model.frame cor pnorm pt vcov
##' @importFrom multcomp glht mcp adjusted
##' @export twosigma

#zi_covar, mean_covar matrices must be specified? Need to figure out what to do exactly here
# should know what theyre doing for that though
twosigma<-function(count_matrix,mean_covar,zi_covar
  ,mean_re=TRUE,zi_re=TRUE
  ,id,adhoc=TRUE,adhoc_thresh=0.1
  ,return_summary_fits=TRUE
  ,disp_covar=NULL #need to be able to use data option?
  ,weights=rep(1,ncol(count_matrix))
  ,control = glmmTMBControl(),ncores=1,cluster_type="Fork",chunk_size=1,lb=FALSE){
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_covar","zi_covar","id")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}
  ngenes<-nrow(count_matrix)
  genes<-rownames(count_matrix)
  if(is.null(genes)){genes<-1:ngenes}
  fit<-vector('list',length=ngenes)
  print("Running Gene-Level Models")
  #browser()
  fit_ts<-function(chunk,id,adhoc){
    k<-0
    num_err=0
    fit<-vector('list',length=length(chunk))
    gene_err<-rep(NA,length(chunk))
    logLik<-rep(NA,length(chunk))
    adhoc_include_RE<-rep(NA,length(chunk))
    for(l in unlist(chunk)){
      if(num_err>0){break}
      k<-k+1
      #setTxtProgressBar(pb, i)
      count<-count_matrix[l,,drop=FALSE]
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
          #message("adhoc method used to set both mean_re and zi_re to TRUE. Set adhoc=FALSE to customize mean_re and zi_re.")
          adhoc_include_RE[k]<-TRUE
        }else{
          mean_re=FALSE
          zi_re=FALSE
          adhoc_include_RE[k]<-FALSE
          #message("adhoc method used to set both mean_re and zi_re to FALSE. Set adhoc=FALSE to customize mean_re and zi_re.")
        }
      }
      formulas<-create_model_formulas(mean_covar,zi_covar
                                      ,mean_form=NULL,zi_form=NULL
                                      ,mean_re,zi_re
                                      ,disp_covar)
      fit[[k]]<-glmmTMB(formula=formulas$mean_form
                 ,ziformula=formulas$zi_form
                 ,weights=weights
                 ,dispformula = formulas$disp_form
                 ,family=nbinom2,verbose = F
                 ,control = control)
      names<-rownames(fit$coefficients$cond)
      if(sum(grepl("Intercept",names))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the mean model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the mean model are:",names),collapse=" "))}
      names_zi<-rownames(fit$coefficients$zi)
      if(!is.null(names_zi)&sum(grepl("Intercept",names_zi))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the ZI model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the ZI model are:",names_zi),collapse=" "))}
      if(return_summary_fits==TRUE){
        fit[[k]]<-summary(fit)
        logLik[k]<-as.numeric(fit[[k]]$logLik)
      }else{
        logLik[k]<-as.numeric(summary(fit[[k]])$logLik)
      }
      gene_err[k]<-(is.na(logLik[k]) | fit[[k]]$sdr$pdHess==FALSE)
    }
    return(list(fit=fit,gene_err=gene_err,adhoc_include_RE=adhoc_include_RE))
  }
  size=chunk_size
  chunks<-split(1:ngenes,ceiling(seq_along(genes)/size))
  nchunks<-length(chunks)
  cl=NULL
  if(cluster_type=="Sock" & ncores>1){
    cl<-parallel::makeCluster(ncores)
    registerDoParallel(cl)
    vars<-c("mean_covar","zi_covar")
    #vars<-vars[!vars=="id"]
    clusterExport(cl,varlist=vars,envir=environment())
  }
  if(cluster_type=="Fork"&ncores>1){
    cl<-parallel::makeForkCluster(ncores,outfile="")
    registerDoParallel(cl)
  }
  pboptions(type="timer")
  if(lb==TRUE){pboptions(use_lb=TRUE)}
  #browser()
  a<-pblapply(chunks,FUN=fit_ts,id=id,adhoc=adhoc,cl=cl)
  if(ncores>1){parallel::stopCluster(cl)}
  rm(count_matrix)
  gc()
  #browser()
  fit<-vector('list',length=ngenes)
  adhoc_include_RE<-rep(NA,length=ngenes)
  for(i in 1:nchunks){
    for(l in chunks[[i]]){
      if(!a[[i]]$gene_err[1]){
        fit[[l]]<-a[[i]]$fit[[1]]
        adhoc_include_RE[l]<-a[[i]]$adhoc_include_RE[[1]]
      }else{
        fit[[l]]<-"Model Fit Error"
      }
      # Remove fits we have used to prevent needing to store more than necessary in memory
      a[[i]]$gene_err<-a[[i]]$gene_err[-1]
      a[[i]]$adhoc_include_RE<-a[[i]]$adhoc_include_RE[-1]
      a[[i]]$fit<-a[[i]]$fit[-1]
    }
  }

  names(fit)<-genes
  return(list(fit=fit,adhoc_include_RE=adhoc_include_RE))
}


