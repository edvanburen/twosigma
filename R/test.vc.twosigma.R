##' Convenient wrapper function for performing (joint) likelihood ratio tests of variance components using the TWO-SIGMA model.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns corresponding to cells. It is recommended to make the rownames the gene names for better output.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column) or a vector if a single covariate is being tested.
##' @param mean_re Should random intercepts be tested in the (conditional) mean model?
##' @param zi_re Should random intercepts be tested in the zero-inflation model?
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param return_full_fits If TRUE, fit objects of class glmmTMB are returned. If FALSE, only objects of class summary.glmmTMB are returned. The latter require a much larger amount of memory to store.
##' @param adhoc Should the adhoc method be used by default to judge if random effects are needed?
##' @param adhoc_thresh Value below which the adhoc p-value is deemed significant (and thus RE are deemed necessary). Only used if adhoc==TRUE.
##' @param silent If TRUE, progress is not printed.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.  See \code{?glmmTMBControl}.
##' @param control Control parameters for optimization in glmmTMB.
##' @param ncores Number of cores used for parallelization. Defaults to 1, meaning no parallelization of any kind is done.
##' @param cluster_type Whether to use a "cluster of type "Fork" or "Sock". On Unix systems, "Fork" will likely improve performance. On Windows, only "Sock" will actually result in parallelized computing.
##' @param chunk_size Number of genes to be sent to each parallel environment. Parallelization is more efficient, particularly with a large count matrix, when the count matrix is 'chunked' into some common size (e.g. 10, 50, 200). Defaults to 10.
##' @param lb Should load balancing be used for parallelization? Users will likely want to set to FALSE for improved performance.
##' @section Details:
##' If either model fails to converge, or the LR statistic is negative, both the statistic and p-value are assigned as NA.
##' @return A list with the following elements:
##' \itemize{
##' \item{\code{fit_null: }} Model fits under the null hypothesis. If \code{return_summary_fits=TRUE}, returns a list of objects of class \code{summary.glmmTMB}. If \code{return_summary_fits=FALSE}, returns a list of model fit objects of class \code{glmmTMB}. In either case, the order matches the row order of \code{count_matrix}, and the names of the list elements are taken as the rownames of \code{count_matrix}.
##' \item{\code{fit_alt: }} Model fits under the alt hypothesis of the same format as \code{fit_null}.
##' \item{\code{LR_stat: }} Vector of Likelihood Ratio statistics. A value of 'NA' implies a convergence issue or other model fit problem.
##' \item{\code{LR_p.val: }} Vector of Likelihood Ratio p-values. A value of 'NA' implies a convergence issue or other model fit problem.
##' }
##' @examples
##' # Set Parameters to Simulate Some Data
##'
##'nind<-10;ncellsper<-rep(50,nind)
##'sigma.a<-.5;sigma.b<-.5;phi<-.1
##'alpha<-c(1,0,-.5,-2);beta<-c(2,0,-.1,.6)
##'beta2<-c(2,1,-.1,.6)
##'id.levels<-1:nind;nind<-length(id.levels)
##'id<-rep(id.levels,times=ncellsper)
##'sim.seed<-1234
##'
##' # Simulate individual level covariates
##'
##'t2d_sim<-rep(rbinom(nind,1,p=.4),times=ncellsper)
##'cdr_sim<-rbeta(sum(ncellsper),3,6)
##'age_sim<-rep(sample(c(20:60),size=nind,replace = TRUE),times=ncellsper)
##'
##'# Construct design matrices
##'
##'Z<-cbind(scale(t2d_sim),scale(age_sim),scale(cdr_sim))
##'colnames(Z)<-c("t2d_sim","age_sim","cdr_sim")
##'X<-cbind(scale(t2d_sim),scale(age_sim),scale(cdr_sim))
##'colnames(X)<-c("t2d_sim","age_sim","cdr_sim")
##'
##' # Simulate Data
##'
##'sim_dat<-matrix(nrow=2,ncol=sum(ncellsper))
##'for(i in 1:nrow(sim_dat)){
##'    sim_dat[i,]<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta2
##'    ,phi,sigma.a,sigma.b,id.levels=NULL)$Y
##'}
##'rownames(sim_dat)<-paste("Gene",1:2)
##'
##' # Run test.vc.twosigma
##'
##' test.vc.twosigma(sim_dat[1,,drop=FALSE],mean_covar = X,zi_covar=Z
##' ,mean_re = TRUE,zi_re=FALSE,id = id)

##' @export test.vc.twosigma
##'
test.vc.twosigma<-function(count_matrix,mean_covar,zi_covar,mean_re=TRUE,zi_re=TRUE,
id,return_full_fits=TRUE,adhoc=FALSE,adhoc_thresh=0.1
,silent=FALSE
,disp_covar=NULL
,weights=rep(1,ncol(count_matrix))
,control = glmmTMBControl(),ncores=1,cluster_type="Fork",chunk_size=1,lb=FALSE)
#,control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e5
# ,step.max=.00001,step.min=.00001
#  ,rel.tol=1e-5,x.tol=1e-5)))
{
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_covar","zi_covar","id")

  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}
  if(mean_re==FALSE & zi_re==FALSE){stop("Please set either mean_re or zi_re to TRUE to use this function.")}
  ngenes<-nrow(count_matrix)
  genes<-rownames(count_matrix)
  if(is.null(genes)){genes<-1:ngenes}
  LR_stat<-rep(NA,length=ngenes)
  p.val<-rep(NA,length=ngenes)
  # Need to duplicate originals because we construct
  # LR test manually and in doing so drop
  # tested covariates each time
  mc<-mean_covar
  zc<-zi_covar
  df_test<-as.numeric(mean_re+zi_re)
  fit_lr<-function(chunk,id,adhoc){
    k<-0
    num_err=0
    f_n<-vector('list',length=length(chunk))
    f_a<-vector('list',length=length(chunk))
    gene_err<-rep(NA,length(chunk))
    logLik<-rep(NA,length(chunk))
    for(l in unlist(chunk)){
      k<-k+1
      mean_covar<-mc
      zi_covar<-zc
      count<-count_matrix[l,,drop=FALSE]
      check_twosigma_input(count,mean_covar,zi_covar
                           ,mean_re,zi_re
                           ,disp_covar,id=id,adhoc=adhoc)
      count<-as.numeric(count)
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
          stop("Mean covariate is not identical to zi covariate. This function is only designed for cases in which the same covariate is being designed in both components.")
        }
      }
      fit_alt<-glmmTMB(formula=formulas$mean_form
                       ,ziformula=formulas$zi_form
                       ,weights=weights
                       ,dispformula = formulas$disp_form
                       ,family=nbinom2,verbose = F
                       ,control = control)
      rm(formulas)
      mean_re_null<-FALSE
      zi_re_null<-FALSE
      formulas<-create_model_formulas(mean_covar,zi_covar
                                      ,mean_form=NULL,zi_form=NULL
                                      ,mean_re_null,zi_re_null
                                      ,disp_covar)
      #If numeric we are assuming that the variable is in the same position
      # need also to point out that users have some responsibilities here


      fit_null<-glmmTMB(formula=formulas$mean_form
                        ,ziformula=formulas$zi_form
                        ,weights=weights
                        ,dispformula = formulas$disp_form
                        ,family=nbinom2,verbose = F
                        ,control = control)
      tryCatch({
        LR_stat[k]<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
        if(LR_stat[k]<0 | (!fit_alt$sdr$pdHess) | (!fit_null$sdr$pdHess)){
          LR_stat[k]<-NA}
        p.val[k]<-1-pchisq(LR_stat[k],df=df_test)},error=function(e){})

      if(return_full_fits==TRUE){
        f_n[[k]]<-fit_null
        f_a[[k]]<-fit_alt
      }else{
        tryCatch({
          f_n[[k]]<-summary(fit_null)
          f_a[[k]]<-summary(fit_alt)
        })
      }
    }
    return(list(fit_null=f_n,fit_alt=f_a,p.val=p.val,LR_stat=LR_stat))
  }
  #browser()
  size=chunk_size
  chunks<-split(1:ngenes,ceiling(seq_along(genes)/size))
  nchunks<-length(chunks)
  cl=NULL
  if(cluster_type=="Sock" & ncores>1){
    cl<-parallel::makeCluster(ncores)
    registerDoParallel(cl)
    vars<-c("mean_covar","zi_covar","mc","zc")
    #vars<-vars[!vars=="id"]
    clusterExport(cl,varlist=vars,envir=environment())
  }
  if(cluster_type=="Fork"&ncores>1){
    cl<-parallel::makeForkCluster(ncores,outfile="")
    registerDoParallel(cl)
  }
  pboptions(type="timer")
  if(lb==TRUE){pboptions(use_lb=TRUE)}
  print("Running Gene-Level Models")
  a<-pblapply(chunks,FUN=fit_lr,id=id,adhoc=adhoc,cl=cl)
  if(ncores>1){parallel::stopCluster(cl)}

  fit_null<-vector('list',length=ngenes)
  fit_alt<-vector('list',length=ngenes)
  p.val<-rep(NA,length=ngenes)
  LR_stat<-rep(NA,length=ngenes)
  for(i in 1:nchunks){
    for(l in chunks[[i]]){
      fit_null[[l]]<-a[[i]]$fit_null[[1]]
      fit_alt[[l]]<-a[[i]]$fit_alt[[1]]
      p.val[l]<-a[[i]]$p.val[1]
      LR_stat[l]<-a[[i]]$LR_stat[1]
      # Remove fits we have used to prevent needing to store more than necessary in memory
      a[[i]]$fit_null<-a[[i]]$fit_null[-1]
      a[[i]]$fit_alt<-a[[i]]$fit_alt[-1]
      a[[i]]$p.val<-a[[i]]$p.val[-1]
      a[[i]]$LR_stat<-a[[i]]$LR_stat[-1]
    }
  }

  names(p.val)<-genes
  names(LR_stat)<-genes
  names(fit_null)<-genes
  names(fit_alt)<-genes
  return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,LR_p.val=p.val))


}
