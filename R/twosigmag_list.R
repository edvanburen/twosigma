##' Gene set testing adjusting for inter-gene correlation using the TWO-SIGMA model of ...
##' @param count_matrix Matrix of non-negative integer read counts. If specifying custom formula(s) via the arguments mean_form and zi_form the expression in mean_form will supersede.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or = 1 to indicate an intercept only model.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column), = 1 to indicate an intercept only model, or = 0 to indicate no zero-inflation model desired.
##' @param mean_re Should random intercepts be included in the (conditional) mean model? Ignored if adhoc=TRUE.
##' @param zi_re Should random intercepts be included in the zero-inflation model? Ignored if adhoc=TRUE.
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param rho Inter-gene correlation value. If NULL (default), estimated using model residuals from TWO-SIGMAG.
##' @param adhoc Should the adhoc method be used by default to judge if random effects are needed? Defaults to FALSE. Setting to TRUE means that different genes could have different models used, which may not be desirable for interpretability.
##' @param adhoc_thresh Value below which the adhoc p-value is deemed significant (and thus RE are deemed necessary). Only used if adhoc==TRUE.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model. Defaults to NULL for constant dispersion.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return An object of class \code{glmmTMB}.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm
##' @export twosigmag_list

twosigmag_list<-function(count_matrix,index_test,index_ref=NULL,contrast,mean_covar,zi_covar
  ,mean_re=TRUE,zi_re=TRUE
  ,id,rho=NULL,adhoc=TRUE,adhoc_thresh=0.1
  ,disp_covar=NULL #need to be able to use data option?
  ,verbose_output=FALSE
  ,weights=rep(1,length(count_matrix[1,]))
  ,control = glmmTMBControl()){

  if(adhoc==TRUE){
    message("adhoc method is allowed but discouraged for gene set testing because statistics for different genes can be based on different models. Users may not wish this to occur.")
  }
  if(is.list(index_test)){
    nsets<-length(index_test)
  }else{
    nsets<-1
  }
  if(!is.null(index_ref)){
    genes<-unique(c(unlist(index_test),unlist(index_ref)))
    ngenes<-length(genes)
    ref_inputted<-TRUE
    if(is.list(index_ref)){
      if(length(index_ref)==1){
        index_ref<-rep(index_ref,nsets)
      }
    }else{
      index_ref<-rep(list(index_ref),nsets)
    }
  }else { # will need all genes
    genes<-1:nrow(count_matrix)
    ngenes<-length(genes)
    ref_inputted<-FALSE
  }
  if(!is.null(index_ref)){
    for(i in 1:nsets){
      if(sum(index_test[[i]]%in%index_ref[[i]])>0){stop(paste("A gene should not be in both the test and reference sets. Check element number",i,"in test set or reference set."))}
    }
    }

  # index_ref must be NULL, or equal in length to index_test? length 1? this would allow for a gene to be in both unless careful though
  if(max(unlist(index_test))>nrow(count_matrix) | min(unlist(index_test))<1){stop("Test Index seems to be invalid, must be numeric within the dimensions of the input count_matrix")}
  ncells<-ncol(count_matrix)

  # Fit all gene level statistics that are needed
  fits_twosigmag<-vector('list',length=nrow(count_matrix))
  residuals_all<-matrix(nrow=nrow(count_matrix),ncol=ncells)
  stats_all<-rep(NA,length=nrow(count_matrix))
  for(i in 1:ngenes){
    l<-genes[i]
    fits_twosigmag[[l]]<-lr.twosigma(count_matrix[l,],contrast = contrast
      ,mean_covar=mean_covar,zi_covar=zi_covar
      ,mean_re=mean_re,zi_re=zi_re
      ,id=id,adhoc=adhoc)

    residuals_all[l,]<-residuals(fits_twosigmag[[l]]$fit_alt)
    stats_all[l]<-fits_twosigmag[[l]]$LR_stat
    print(paste("Finished Gene Number",i,"of",ngenes))
  }
  stats_test<-vector('list',length=nsets)
  stats_ref<-vector('list',length=nsets)
  p.val<-numeric(length=nsets)
  rho_est<-numeric(length=nsets)
  if(is.null(index_ref)){index_ref<-vector('list',length=nsets)}
  for(i in 1:nsets){
    if(is.list(index_test)){
      stats_test[[i]]<-stats_all[index_test[[i]]]
      test_size<-length(index_test[[i]])
      residuals_test<-residuals_all[index_test[[i]],]
    }else{
      stats_test[[i]]<-stats_all[index_test]
      test_size<-length(index_test)
      residuals_test<-residuals_all[index_test,]
    }

    if(ref_inputted==FALSE){
    index_ref[[i]]<-setdiff(1:nrow(count_matrix),index_test[[i]])
    stats_ref[[i]]<-stats_all[index_ref[[i]]]
    ref_size<-length(index_ref[[i]])
    }else{
      if(is.list(index_ref)){
        stats_ref[[i]]<-stats_all[index_ref[[i]]]
        ref_size<-length(index_ref[[i]])
      }else{
        stats_ref[[i]]<-stats_all[index_ref]
        ref_size<-length(index_ref)
      }
    }

    if(!is.null(rho)){rho_est[i]<-rho}
    if(is.null(rho)){
      print(paste("Estimating Set-Level correlation and calculating p-value"))
      nind<-length(unique(id))
      cor_temp<-numeric(length=nind)
      unique_id<-unique(id)
      j<-0
      for(y in unique_id){
        j<-j+1
        temp<-cor(t(residuals_test[,which(id==y)])) # any checks for id behavior needed here?
        cor_temp[j]<-mean(temp[upper.tri(temp)],na.rm = T)
      }
      rho_est[i]<-mean(cor_temp)
    }
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho_est[i])+(test_size-1)*asin((rho_est[i]+1)/2))
    wilcox_stat<-sum(rank(c(stats_test[[i]],stats_ref[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    p.val[i]<-2*pnorm(-1*abs((wilcox_stat-.5*test_size*ref_size)/sqrt(var)))
  }

   if(verbose_output==TRUE){
     return(list(gene_level_fits=fits_twosigmag,LR_stats_gene_level_all=stats_all,set_p.val=p.val,corr=rho_est,index_test=index_test))
   }else{
    return(list(LR_stats_gene_level_all=stats_all,set_p.val=p.val,corr=rho_est,index_test=index_test))
  }
}
