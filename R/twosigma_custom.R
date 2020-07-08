##' Fit the TWO-SIGMA model with custom user-specified model formulas.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns correspoding to cells. It is recommended to make the rownames the gene names for better output.
##' @param mean_form Custom two-sided model formula for the (conditional) mean model. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the LHS of the formula begins with "count."
##' @param zi_form Custom one-sided model formula for the zero-inflation model. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param id Vector of individual-level (sample-level) ID's. Used for random effect prediction but required regardless of their presence in the model.
##' @param return_summary_fits If TRUE, the package returns a \code{summary.glmmTMB} object for each gene.  If FALSE, a \code{glmmTMB} object is returned for each gene. The latter requires far more storage space.
##' @param silent If TRUE, progress is not printed.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.  See \code{?glmmTMBControl}.
##' @param ncores Number of cores used for parallelization. Defaults to 1, meaning no parallelization of any kind is done.
##' @param cluster_type Whether to use a "cluster of type "Fork" or "Sock". On Unix systems, "Fork" will likely improve performance. On Windows, only "Sock" will actually result in parallelized computing.
##' @param chunk_size Number of genes to be sent to each parallel environment. Parallelization is more efficient, particuarly with a large count matrix, when the count matrix is 'chunked' into some common size (e.g. 10, 50, 200). Defaults to 10.
##' @param lb Should load balancing be used for parallelization? Users will likely want to set to FALSE for improved performance.
##' @param internal_call Not needed by users called \code{twosigma_custom} directly.
##' @section Details:
##' This function is likely only needed if users wish to include random effect terms beyond random intercepts. Users should be confident in their abilities to specify random effects using the syntax of lme4.
##' @return A list with the following elements:
##' \itemize{
##' \item{\code{fit: }} If \code{return_summary_fits=TRUE}, returns a list of model fit objects of class \code{summary.glmmTMB}. If \code{return_summary_fits=FALSE}, returns a list of model fit objects of class \code{glmmTMB}. In either case, the order matches the row order of \code{count_matrix}, and the names of the list elements are taken as the rownames of \code{count_matrix}.
##' }
##' @export twosigma_custom
twosigma_custom<-function(count_matrix,mean_form,zi_form,id,return_summary_fits=TRUE,
  silent=FALSE,disp_covar=NULL
  ,weights=rep(1,ncol(count_matrix))
  ,control = glmmTMBControl(),ncores=1,cluster_type="Fork"
  ,chunk_size=1,lb=FALSE,internal_call=FALSE){

  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_form","zi_form","id")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  ngenes<-nrow(count_matrix)
  genes<-rownames(count_matrix)
  if(is.null(genes)){genes<-1:ngenes}
  #if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}




  fit_ts<-function(chunk,id){
    k<-0
    num_err=0
    fit<-vector('list',length=length(chunk))
    gene_err<-rep(NA,length(chunk))
    logLik<-rep(NA,length(chunk))
    for(l in unlist(chunk)){
      if(num_err>0){break}
      k<-k+1
    count<-count_matrix[l,,drop=FALSE]
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
    if(return_summary_fits==TRUE){
      fit[[k]]<-summary(f)
      logLik[k]<-as.numeric(fit[[k]]$logLik)
      gene_err[k]<-(is.na(logLik[k]) | f$sdr$pdHess==FALSE)
    }else{
      fit[[k]]<-f
      logLik[k]<-as.numeric(summary(fit[[k]])$logLik)
      gene_err[k]<-(is.na(logLik[k]) | f$sdr$pdHess==FALSE)
    }
    return(list(fit=fit,gene_err=gene_err))
    }
  }
  #print("Running Gene-Level Models")
  size=chunk_size
  chunks<-split(1:ngenes,ceiling(seq_along(genes)/size))
  nchunks<-length(chunks)
  cl=NULL
  if(cluster_type=="Sock" & ncores>1){
    cl<-parallel::makeCluster(ncores)
    registerDoParallel(cl)
    vars<-unique(c(all.vars(mean_form)[-1],all.vars(zi_form)[-1]))
    #vars<-vars[!vars=="id"]
    clusterExport(cl,varlist=vars,envir=environment())
  }
  if(cluster_type=="Fork"&ncores>1){
    cl<-parallel::makeForkCluster(ncores,outfile="")
    registerDoParallel(cl)
  }
  if(internal_call==FALSE){
    print("Running Gene-Level Models")
    pboptions(type="timer")
    if(lb==TRUE){pboptions(use_lb=TRUE)}
    if(silent==TRUE){pboptions(type="none")}
    #browser()
    a<-pblapply(chunks,FUN=fit_ts,id=id,cl=cl)
    if(ncores>1){parallel::stopCluster(cl)}
  }else{
    a<-lapply(chunks,FUN=fit_ts,id=id)
  }

  rm(count_matrix)
  gc()
  #browser()
  fit<-vector('list',length=ngenes)
  adhoc_include_RE<-rep(NA,length=ngenes)
  for(i in 1:nchunks){
    for(l in chunks[[i]]){
      #if(!a[[i]]$gene_err[1]){
        fit[[l]]<-a[[i]]$fit[[1]]
      #}else{
      #  fit[[l]]<-"Model Fit Error"
     # }
      # Remove fits we have used to prevent needing to store more than necessary in memory
      a[[i]]$gene_err<-a[[i]]$gene_err[-1]
      a[[i]]$fit<-a[[i]]$fit[-1]
    }
  }
  names(fit)<-genes
  return(fit)
}

