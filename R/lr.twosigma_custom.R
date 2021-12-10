##' Convenient wrapper function for performing joint likelihood ratio tests with the TWO-SIGMA model using custom user-specified formulas.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns corresponding to cells. It is recommended to make the rownames the gene names for better output.
##' @param mean_form_alt Custom two-sided model formula for the (conditional) mean model under the null. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the dependent variable matches the argument to the parameter "count."
##' @param zi_form_alt Custom one-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param mean_form_null Custom two-sided model formula for the (conditional) mean model under the null. Syntax is as in \code{mean_form_alt}.
##' @param zi_form_null Custom one-sided model formula for the zero-inflation model under the null. Syntax is as in \code{zi_form_alt}.
##' @param id Vector of individual-level (sample-level) ID's. Used for random effect prediction but required regardless of their presence in the model.
##' @param lr.df Degrees of Freedom for the constructed likelihood ratio test. Must be a non-negative integer.
##' @param return_full_fits If TRUE, full fit objects of class glmmTMB are returned.  If FALSE, only fit objects of class summary.glmmTMB are returned.  The latter requires far less memory to store.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix or = 1 to indicate an intercept only model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in glmmTMB.  See \code{?glmmTMBControl}.
##' @param ncores Number of cores used for parallelization. Defaults to 1, meaning no parallelization of any kind is done.
##' @param cluster_type Whether to use a "cluster of type "Fork" or "Sock". On Unix systems, "Fork" will likely improve performance. On Windows, only "Sock" will actually result in parallelized computing.
##' @param chunk_size Number of genes to be sent to each parallel environment. Parallelization is more efficient, particularly with a large count matrix, when the count matrix is 'chunked' into some common size (e.g. 10, 50, 200). Defaults to 10.
##' @param lb Should load balancing be used for parallelization? Users will likely want to set to FALSE for improved performance.
##' @param internal_call Not needed by users called \code{lr.twosigma_custom} directly.
##' @section Details:
##' This function is a wrapper for conducting fixed effect likelihood ratio tests with twosigma.  There is no checking to make sure that the alt and null model formulas represent a valid likelihood ratio test when fit together.  Users must ensure that inputted formulas represent valid nested models. If either model fails to converge, or the LR statistic is negative, both the statistic and p-value are assigned as NA.
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
##' # Run lr.twosigma_custom
##'
##' lr.twosigma_custom(count=sim_dat[1,,drop=FALSE]
##' ,mean_form_alt = count~X,mean_form_null = count~X[,-1]
##' ,zi_form_alt = ~0,zi_form_null = ~0,id=id,lr.df=1)
##' @export lr.twosigma_custom

lr.twosigma_custom<-function(count_matrix,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null
                      ,id,lr.df,return_full_fits=TRUE,
                      disp_covar=NULL
                      ,weights=rep(1,ncol(count_matrix))
                      ,control = glmmTMBControl(),ncores=1,cluster_type="Fork",chunk_size=10,lb=FALSE,internal_call=FALSE)
  {
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","mean_form_alt","zi_form_alt","mean_form_null","zi_form_null","id","lr.df")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}
  ngenes<-nrow(count_matrix)
  genes<-rownames(count_matrix)
  if(is.null(genes)){genes<-1:ngenes}
  sum_fit_alt<-vector('list',length=ngenes)
  sum_fit_null<-vector('list',length=ngenes)

  if(return_full_fits==TRUE){
    fits_all_null<-vector('list',length=ngenes)
    fits_all_alt<-vector('list',length=ngenes)
  }
  fit_lr<-function(chunk,id){
    k<-0
    num_err=0
    f_n<-vector('list',length=length(chunk))
    f_a<-vector('list',length=length(chunk))
    gene_err<-rep(NA,length(chunk))
    logLik<-rep(NA,length(chunk))
    LR_stat<-rep(NA,length=length(chunk))
    p.val<-rep(NA,length=length(chunk))
    for(l in unlist(chunk)){
      k<-k+1
      count<-count_matrix[l,,drop=FALSE]
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
      tryCatch({
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
        LR_stat[k]<- as.numeric(-2*(summary(fit_null)$logLik-summary(fit_alt)$logLik))
        if(!is.na(LR_stat[k])& (LR_stat[k]<0 | (!fit_alt$sdr$pdHess) | (!fit_null$sdr$pdHess)
                                |is.na(fit_alt$logLik)|is.na(fit_null$logLik))){
          LR_stat[k]<-NA}},error=function(e){}) # end try-catch
      p.val[k]<-1-pchisq(LR_stat[k],df=lr.df)
      if(return_full_fits==TRUE){
        f_n[[k]]<-fit_null
        f_a[[k]]<-fit_alt
      }else{
        tryCatch({
          f_n[[k]]<-summary(fit_null)
          f_a[[k]]<-summary(fit_alt)
        })
      }},error=function(e){})
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
  vars<-unique(c(all.vars(mean_form_alt)[-1],all.vars(zi_form_alt)
                 ,all.vars(mean_form_null)[-1],all.vars(zi_form_null)))
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
  a<-pblapply(chunks,FUN=fit_lr,id=id,cl=cl)
  if(ncores>1){parallel::stopCluster(cl)}
}else{
a<-lapply(chunks,FUN=fit_lr,id=id)
}
#browser()
fit_null<-do.call(c,sapply(a,'[',1))
fit_alt<-do.call(c,sapply(a,'[',2))
p.val<-unlist(sapply(a,'[',3))
LR_stat<-unlist(sapply(a,'[',4))

names(p.val)<-genes
names(LR_stat)<-genes
names(fit_null)<-genes
names(fit_alt)<-genes
return(list(fit_null=fit_null,fit_alt=fit_alt,LR_stat=LR_stat,LR_p.val=p.val))
}
