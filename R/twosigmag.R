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
##' @export twosigmag

twosigmag<-function(count_matrix,index_test,index_ref=NULL,contrast,mean_covar,zi_covar
  ,mean_re=TRUE,zi_re=TRUE
  ,id,rho=NULL,adhoc=TRUE,adhoc_thresh=0.1
  ,disp_covar=NULL #need to be able to use data option?
  ,verbose_output=FALSE
  ,weights=rep(1,length(count_matrix[1,]))
  ,control = glmmTMBControl()){

  if(adhoc==TRUE){
  message("adhoc method is allowed but discouraged for gene set testing because statistics for different genes can be based on different models. Users may not wish this to occur.")
  }
  if(!is.null(index_ref)){
    genes<-c(index_test,index_ref)
    ngenes<-length(genes)
  }else {
    genes<-c(index_test,setdiff(1:nrow(count_matrix),index_test))
    ngenes<-length(genes)
  }
  if(!is.null(index_ref)){
    if(sum(index_test%in%index_ref)>0){stop("A gene should not be in both the test and reference sets")}
  }
  if(max(index_test)>ngenes | min(index_test)<1){stop("Test Index seems to be invalid, must be numeric within the dimensions of the input count_matrix")}
  ncells<-ncol(count_matrix)
  test_size<-length(index_test)
  ref_size<-ngenes-test_size
  fits_twosigmag<-vector('list',length=ngenes)
  residuals_test<-matrix(nrow=length(index_test),ncol=ncells)
  stats_test<-numeric(length=test_size)
  stats_ref<-numeric(length=ngenes-test_size)
  stats_all<-numeric(length=ngenes)
  j<-0 # Index for test set
  k<-0 # Index for ref set
    for(i in 1:ngenes){
      l<-genes[i]
      fits_twosigmag[[i]]<-lr.twosigma(count_matrix[l,],contrast = contrast
        ,mean_covar=mean_covar,zi_covar=zi_covar
        ,mean_re=mean_re,zi_re=zi_re
        ,id=id,adhoc=adhoc)
      if(l%in%index_test){
        j<-j+1
        residuals_test[j,]<-residuals(fits_twosigmag[[i]]$fit_alt)
        stats_test[j]<-fits_twosigmag[[i]]$LR_stat
      }else{# all genes not in test set are taken as reference set for now
        k<-k+1
        stats_ref[k]<-fits_twosigmag[[i]]$LR_stat
      }
      stats_all[i]<-fits_twosigmag[[i]]$LR_stat
      print(paste("Finished Gene Number",i,"of",ngenes))
    }
      if(is.null(rho)){
        print(paste("Estimating Set-Level correlation and calculating p-value"))
        nind<-length(unique(id))
        cor_temp<-numeric(length=nind)
        m<-1
        unique_id<-unique(id)
        j<-0
        for(y in unique_id){
          j<-j+1
          temp<-cor(t(residuals_test[,which(id==y)])) # any checks for id behavior needed here?
          cor_temp[j]<-mean(temp[upper.tri(temp)],na.rm = T)
        }
        rho<-mean(cor_temp)
      }
      var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho)+(test_size-1)*asin((rho+1)/2))

      wilcox_stat<-sum(rank(c(stats_test,stats_ref))[1:test_size]) - .5*test_size*(test_size+1)
      p.val<-2*pnorm(-1*abs((wilcox_stat-.5*test_size*ref_size)/sqrt(var)))
      if(verbose_output==TRUE){
        return(list(fits_twosigma=fits_twosigmag,stats_all=stats_all,set_p.val=p.val,corr=rho))
      }else{
        return(list(LR_stats_all=stats_all,set_p.val=p.val,corr=rho))
      }
}
