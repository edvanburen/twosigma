##' Convenient wrapper function for performing joint likelihood ratio tests using the TWO-SIGMA model.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns corresponding to cells. It is recommended to make the rownames the gene names for better output.
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
##' @param ncores Number of cores used for parallelization. Defaults to 1, meaning no parallelization of any kind is done.
##' @param cluster_type Whether to use a "cluster of type "Fork" or "Sock". On Unix systems, "Fork" will likely improve performance. On Windows, only "Sock" will actually result in parallelized computing.
##' @param chunk_size Number of genes to be sent to each parallel environment. Parallelization is more efficient, particularly with a large count matrix, when the count matrix is 'chunked' into some common size (e.g. 10, 50, 200). Defaults to 10.
##' @param lb Should load balancing be used for parallelization? Users will likely want to set to FALSE for improved performance.
##' @section Details:
##' This function assumes that the variable being tested is in both components of the model (and thus that the zero-inflation component exists and contains more than an Intercept). Users wishing to do fixed effect testing in other cases or specify custom model formulas they will need to construct the statistics themselves using either two separate calls to \code{twosigma} or the \code{lr.twosigma_custom} function. If \code{adhoc=TRUE}, any input in mean_re and zi_re will be ignored. If either model fails to converge, or the LR statistic is negative, both the statistic and p-value are assigned as NA.
##' @return A list with the following elements:
##' \itemize{
##' \item{\code{fit_null: }} Model fits under the null hypothesis. If \code{return_summary_fits=TRUE}, returns a list of objects of class \code{summary.glmmTMB}. If \code{return_summary_fits=FALSE}, returns a list of model fit objects of class \code{glmmTMB}. In either case, the order matches the row order of \code{count_matrix}, and the names of the list elements are taken as the rownames of \code{count_matrix}.
##' \item{\code{fit_alt: }} Model fits under the alt hypothesis of the same format as \code{fit_null}.
##' \item{\code{LR_stat: }} Vector of Likelihood Ratio statistics. A value of 'NA' implies a convergence issue or other model fit problem.
##' \item{\code{LR_p.val: }} Vector of Likelihood Ratio p-values. A value of 'NA' implies a convergence issue or other model fit problem.
##' \item{\code{adhoc_include_RE: }} Logical vector indicator whether the adhoc method determined random effects needed.  If \code{adhoc=F}, then a vector of NA's.
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
##' # Run lr.twosigma
##'
##' lr.twosigma(count=sim_dat[1,,drop=FALSE],mean_covar = X,zi_covar = Z,id=id,covar_to_test = 1)
##' @export lr.twosigma

lr.twosigma<-function(count_matrix,mean_covar,zi_covar,covar_to_test
                      ,mean_re=FALSE,zi_re=FALSE,
                       id,return_full_fits=TRUE,adhoc=FALSE,adhoc_thresh=0.1
                      ,silent=FALSE
                      ,disp_covar=NULL
                      ,weights=rep(1,ncol(count_matrix))
                      ,control = glmmTMBControl(),ncores=1,cluster_type="Fork",chunk_size=10,lb=FALSE)
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
  genes<-rownames(count_matrix)
  if(is.null(genes)){genes<-1:ngenes}
  # Need to duplicate originals because we construct
  # LR test manually and in doing so drop
  # tested covariates each time
  mc<-mean_covar
  zc<-zi_covar

  fit_lr<-function(chunk,id,adhoc){
    k<-0
    num_err=0
    f_n<-vector('list',length=length(chunk))
    f_a<-vector('list',length=length(chunk))
    gene_err<-rep(TRUE,length(chunk))
    logLik<-rep(NA,length(chunk))
    adhoc_include_RE<-rep(NA,length(chunk))
    LR_stat<-rep(NA,length(chunk))
    p.val<-rep(NA,length(chunk))
    for(l in unlist(chunk)){
    k<-k+1
    mean_covar<-mc
    zi_covar<-zc
    count<-count_matrix[l,,drop=FALSE]
    check_twosigma_input(count,mean_covar,zi_covar
                         ,mean_re,zi_re
                         ,disp_covar,id=id,adhoc=adhoc)
    count<-as.numeric(count)
    if(adhoc==TRUE){
      if(is.atomic(zi_covar)&length(zi_covar)==1){
        if(zi_covar==0){stop("adhoc method only implemented when ZI model contains at minimum an intercept. Please either set adhoc=FALSE or specify at minimum an intercept in the ZI model.")}}
      adhoc_p.val<-adhoc.twosigma(count=count,mean_covar=mean_covar,zi_covar = zi_covar,id=id,weights=weights)
      if(adhoc_p.val<adhoc_thresh){
        mean_re=TRUE
        zi_re=TRUE
        #message("adhoc method used to set both mean_re and zi_re to TRUE. Set adhoc=FALSE to customize mean_re and zi_re.")
        adhoc_include_RE[k]<-TRUE
      }else{
        mean_re=FALSE
        zi_re=FALSE
        adhoc_include_RE[k]<-FALSE
        #print("adhoc method used to set both mean_re and zi_re to FALSE. Set adhoc=FALSE to customize user-inputted values for mean_re and zi_re.")
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
        stop("Mean covariate is not identical to zi covariate. This function is only designed for cases in which the same covariate is being designed in both components.")
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
    tryCatch({
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
      LR_stat[k]<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
      if(!is.na(LR_stat[k])& (LR_stat[k]<0 | (!fit_alt$sdr$pdHess) | (!fit_null$sdr$pdHess)
         |is.na(fit_alt$logLik)|is.na(fit_null$logLik))){
        LR_stat[k]<-NA}
      },error=function(e){})
    p.val[k]<-1-pchisq(LR_stat[k],df=2)

    if(return_full_fits==TRUE){
      f_n[[k]]<-fit_null
      f_a[[k]]<-fit_alt
    }else{
      tryCatch({
      f_n[[k]]<-summary(fit_null)
      f_a[[k]]<-summary(fit_alt)
      })
    }},error=function(e){}) # end try-catch
    }
    return(list(fit_null=f_n,fit_alt=f_a,p.val=p.val,LR_stat=LR_stat,adhoc_include_RE=adhoc_include_RE))
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

  #fit_null<-vector('list',length=ngenes)
  #fit_alt<-vector('list',length=ngenes)
  #p.val<-rep(NA,length=ngenes)
  #LR_stat<-rep(NA,length=ngenes)
  #adhoc_include_RE<-rep(NA,length=ngenes)

  #browser()
  fit_null<-do.call(c,sapply(a,'[',1))
  fit_alt<-do.call(c,sapply(a,'[',2))
  p.val<-unlist(sapply(a,'[',3))
  LR_stat<-unlist(sapply(a,'[',4))
  adhoc_include_RE<-unlist(sapply(a,'[',5))

  names(p.val)<-genes
  names(LR_stat)<-genes
  names(fit_null)<-genes
  names(fit_alt)<-genes
  names(adhoc_include_RE)<-genes

  return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,LR_p.val=p.val,adhoc_include_RE=adhoc_include_RE))
}
