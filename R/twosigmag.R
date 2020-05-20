##' Gene set testing adjusting for inter-gene correlation using the TWO-SIGMA model of ...
##' @param count_matrix Matrix of non-negative integer read counts. If specifying custom formula(s) via the arguments mean_form and zi_form the expression in mean_form will supersede.
##' @param index_test Index corresponding to rows of the count matrix that are in the test set. Either a list with each element being a different set or a numeric vector to test only a single set.
##' @param index_ref Index corresponding to rows of the count matrix that are in the reference set.  If NULL, a reference set is randomly selected of the same size as test size using genes not in the test set or using all other genes. See \code{all_as_ref}.
##' @param all_as_ref Should all genes not in the test set be used as the reference? If FALSE, a random subset is taken of size equal to the test size.
##' @param contrast Either a string indicating the column name of the covariate to test or an integer referring to its column position in BOTH the mean_covar and zi_covar matrices (if the two matrices differ using a string name is preferred). Argument is ignored if mean_covar and zi_covar are both a single covariate (that covariate is assumed of interest).
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or = 1 to indicate an intercept only model.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column), = 1 to indicate an intercept only model, or = 0 to indicate no zero-inflation model desired.
##' @param mean_re Should random intercepts be included in the (conditional) mean model? Ignored if adhoc=TRUE.
##' @param zi_re Should random intercepts be included in the zero-inflation model? Ignored if adhoc=TRUE.
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param rho Inter-gene correlation value. If NULL (default), estimated using model residuals from TWO-SIGMAG.
##' @param allow_neg_corr Should negative correlation values be allowed? If FALSE, correlation is set to zero (leads to conservative inference).
##' @param adhoc Should the adhoc method be used by default to judge if random effects are needed? Defaults to FALSE. Setting to TRUE means that different genes could have different models used, which may not be desirable for interpretability.
##' @param adhoc_thresh Value below which the adhoc p-value is deemed significant (and thus RE are deemed necessary). Only used if adhoc==TRUE.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model. Defaults to NULL for constant dispersion.
##' @param return_fits Should complete model fits be returned? Use cautiously because fit objects can be very large.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return An object of class \code{glmmTMB}.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm
##' @export twosigmag

twosigmag<-function(count_matrix,index_test,index_ref=NULL,all_as_ref=FALSE,contrast,mean_covar,zi_covar
  ,mean_re=TRUE,zi_re=TRUE
  ,id,rho=NULL
  ,allow_neg_corr=FALSE
  ,adhoc=FALSE,adhoc_thresh=0.1
  ,disp_covar=NULL #need to be able to use data option?
  ,return_fits=FALSE
  ,weights=rep(1,length(count_matrix[1,]))
  ,control = glmmTMBControl()){

  if(!(adhoc==FALSE)){print("The adhoc method is not recommended for gene set testing due to interpretability.")}


  count_matrix<-as.matrix(count_matrix)
  if(is.list(index_test)){
    nsets<-length(index_test)
    list_lengths<-lapply(index_test,FUN=length)
    if(sum(list_lengths<2)>0){stop("All test sets must have at least two genes. Please remove singleton or empty sets.")}
  }else{
    nsets<-1
    if(length(index_test)<2){stop("All test sets must have at least two genes. Please remove singleton or empty sets.")}
  }
  if(all_as_ref==TRUE & !is.null(index_ref)){stop("Please specify either all_as_ref=TRUE or index_ref as a non-NULL input. If all_as_ref is TRUE, then index_ref must be NULL.")}
  if(!is.null(index_ref)){
    if(is.list(index_ref)){
      for(i in 1:nsets){
        if(sum(index_test[[i]]%in%index_ref[[i]])>0){stop(paste("A gene should not be in both the test and reference sets. Check element number",i,"in test set or reference set."))}
      }

      if(is.list(index_test) & length(index_ref)!=length(index_test)){
        stop("If index_test and index_ref are both lists they should be the same length.")
      }
    }else{
      for(i in 1:nsets){
        if(sum(index_test[[i]]%in%index_ref)>0){stop(paste("A gene should not be in both the test and reference sets. Check element number",i,"in test set or reference set."))}
      }

      index_ref<-rep(list(index_ref),nsets)
    }
    genes<-unique(c(unlist(index_test),unlist(index_ref)))
    ngenes<-length(genes)
    ref_inputted<-TRUE
  }else {# will need to construct reference set
    if(is.list(index_test)){
      index_ref<-vector('list',length=nsets)
      for(i in 1:nsets){
        if(all_as_ref==FALSE){
          index_ref[[i]]<-sample(setdiff(1:nrow(count_matrix),index_test[[i]]),size=length(index_test[[i]]))
        }else{#all_as_ref==TRUE
          index_ref[[i]]<-setdiff(1:nrow(count_matrix),index_test[[i]])
        }
      }
      genes<-unique(c(unlist(index_test),unlist(index_ref)))
    }else{#index_test is just a single gene set
      if(all_as_ref==FALSE){
        index_ref<-sample(setdiff(1:nrow(count_matrix),index_test),size=length(index_test))
      }else{
        index_ref<-setdiff(1:nrow(count_matrix),index_test)
      }

      genes<-unique(c(unlist(index_test),unlist(index_ref)))
    }
    ngenes<-length(genes)
    ref_inputted<-FALSE
  }

  if(max(unlist(index_test))>nrow(count_matrix) | min(unlist(index_test))<1){stop("Test Index seems to be invalid, must be numeric within the dimensions of the input count_matrix")}
  ncells<-ncol(count_matrix)

  # Fit all gene level statistics that are needed
  if(return_fits==FALSE){
    fit_twosigmag<-list()
  }else{
    fit_twosigmag<-vector('list',length=nrow(count_matrix))
  }
  residuals_all<-matrix(nrow=nrow(count_matrix),ncol=ncells)
  stats_all<-rep(NA,length=nrow(count_matrix))
  p.vals_gene_level<-rep(NA,length=nrow(count_matrix))
  avg_logFC_gene_level<-rep(NA,length=nrow(count_matrix))
  #browser()
  for(i in 1:ngenes){
    l<-genes[i]
    if(return_fits==TRUE){
      fit_twosigmag[[l]]<-lr.twosigma(count_matrix[l,],contrast = contrast
        ,mean_covar=mean_covar,zi_covar=zi_covar
        ,mean_re=mean_re,zi_re=zi_re
        ,id=id,adhoc=adhoc)
      residuals_all[l,]<-residuals(fit_twosigmag[[l]]$fit_alt)
      stats_all[l]<-fit_twosigmag[[l]]$LR_stat
      p.vals_gene_level[l]<-fit_twosigmag[[l]]$LR_p.val
      avg_logFC_gene_level[l]<-fit_twosigmag[[l]]$mean_comp_logFC
    }else{
      fit_twosigmag<-lr.twosigma(count_matrix[l,],contrast = contrast
        ,mean_covar=mean_covar,zi_covar=zi_covar
        ,mean_re=mean_re,zi_re=zi_re
        ,id=id,adhoc=adhoc)
      residuals_all[l,]<-residuals(fit_twosigmag$fit_alt)
      stats_all[l]<-fit_twosigmag$LR_stat
      p.vals_gene_level[l]<-fit_twosigmag$LR_p.val
      avg_logFC_gene_level[l]<-fit_twosigmag$mean_comp_logFC
    }

    print(paste("Finished Gene Number",i,"of",ngenes))
  }
  #browser()
  stats_test<-vector('list',length=nsets)
  stats_ref<-vector('list',length=nsets)
  direction<-vector(length=nsets)
  p.val<-numeric(length=nsets)
  p.val_ttest<-numeric(length=nsets)
  rho_est<-numeric(length=nsets)
  if(is.null(index_ref)){index_ref<-vector('list',length=nsets)}
  for(i in 1:nsets){
    if(is.list(index_test)){
      stats_test[[i]]<-stats_all[index_test[[i]]]
      stats_test[[i]]<-stats_test[[i]][!is.na(stats_test[[i]])]
      test_size<-length(stats_test[[i]])
      residuals_test<-residuals_all[index_test[[i]],]
      direction[i]<-ifelse(sign(mean(avg_logFC_gene_level[index_test[[i]]]))==1,"Up","Down")
    }else{
      stats_test[[i]]<-stats_all[index_test]
      stats_test[[i]]<-stats_test[[i]][!is.na(stats_test[[i]])]
      test_size<-length(stats_test[[i]])
      residuals_test<-residuals_all[index_test,]
      direction[i]<-ifelse(sign(mean(avg_logFC_gene_level[index_test]))==1,"Up","Down")
    }

    # if(ref_inputted==FALSE){
    # #index_ref[[i]]<-setdiff(1:nrow(count_matrix),index_test[[i]])
    # stats_ref[[i]]<-stats_all[-index_test[[i]]]
    # ref_size<-length(index_ref[[i]])
    # }else{
      if(is.list(index_ref)){
        stats_ref[[i]]<-stats_all[index_ref[[i]]]
        stats_ref[[i]]<-stats_ref[[i]][!is.na(stats_ref[[i]])]
        ref_size<-length(stats_ref[[i]])
      }else{
        stats_ref[[i]]<-stats_all[index_ref]
        stats_ref[[i]]<-stats_ref[[i]][!is.na(stats_ref[[i]])]
        ref_size<-length(stats_ref[[i]])
      }
   # }

    if(!is.null(rho)){rho_est[i]<-rho}
    if(is.null(rho)){
      #print(paste("Estimating Set-Level correlation and calculating p-value"))
      nind<-length(unique(id))
      cor_temp<-numeric(length=nind)
      unique_id<-unique(id)
      j<-0
      for(y in unique_id){
        j<-j+1
        temp<-cor(t(residuals_test[,which(id==y)])) # any checks for id behavior needed here?
        cor_temp[j]<-mean(temp[upper.tri(temp)],na.rm = T)
      }
      # There will be an NA here if residuals for an individual are zero for a given gene
      rho_est[i]<-mean(cor_temp,na.rm=T)
    }
    #browser()
    if(!allow_neg_corr & rho_est[i]<0){rho_est[i]<-0}
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho_est[i])+(test_size-1)*asin((rho_est[i]+1)/2))
    wilcox_stat<-sum(rank(c(stats_test[[i]],stats_ref[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    p.val[i]<-2*pnorm(-1*abs((wilcox_stat-.5*test_size*ref_size)/sqrt(var)))
    n_genes_tested<-test_size+ref_size
    delta<-n_genes_tested/ref_size*(mean(stats_test[[i]],na.rm=T)-mean(c(stats_test[[i]],stats_ref[[i]]),na.rm=T))
    vif<-1+(test_size-1)*rho_est[i]
    varStatPooled<-((n_genes_tested-1)*var(c(stats_test[[i]],stats_ref[[i]]),na.rm=T)-delta^2*test_size*ref_size/n_genes_tested)/(n_genes_tested-2)
    two.sample.t <- delta / sqrt( varStatPooled * (vif/test_size + 1/ref_size) )
    p.val_ttest[i]<-2*pt(-1*abs(two.sample.t),df=n_genes_tested-2)
    # This only happens if all genes in the test set are NA
    # In this case p-value will report as zero when it's really NA
    if(test_size==0|ref_size==0){p.val[i]<-NA;p.val_ttest[i]<-NA}
    if(i%%100==0){print(paste0("Set ",i," of ",nsets," Finished"))}
  }
  names(p.val)<-names(index_test)
  names(p.val_ttest)<-names(index_test)
  names(stats_all)<-rownames(count_matrix)
  names(p.vals_gene_level)<-rownames(count_matrix)
  names(avg_logFC_gene_level)<-rownames(count_matrix)
  names(rho_est)<-names(index_test)
  names(direction)<-names(index_test)
   if(return_fits==TRUE){
     return(list(gene_level_fits=fit_twosigmag,LR_stats_gene_level_all=stats_all,set_p.val=p.val,set_p.val_ttest=p.val_ttest,direction=direction,p.vals_gene_level=p.vals_gene_level,corr=rho_est,avg_logFC_gene_level=avg_logFC_gene_level,test_sets=index_test,ref_sets=index_ref))
   }else{
    return(list(LR_stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,set_p.val_ttest=p.val_ttest,direction=direction,corr=rho_est,avg_logFC_gene_level=avg_logFC_gene_level,test_sets=index_test,ref_sets=index_ref,set_p.val_ttest=p.val_ttest))
  }
}
