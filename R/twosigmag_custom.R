##' Gene set testing adjusting for inter-gene correlation using the TWO-SIGMA model of ... with custom user-specified model formulas.
##' @param count_matrix Matrix of non-negative integer read counts. If specifying custom formula(s) via the arguments mean_form and zi_form the expression in mean_form will supersede.
##' @param index_test Index corresponding to rows of the count matrix that are in the test set. Either a list with each element being a different set or a numeric vector to test only a single set.
##' @param index_ref Index corresponding to rows of the count matrix that are in the reference set.  If NULL, a reference set is randomly selected of the same size as test size using genes not in the test set or using all other genes. See all_as_ref.
##' @param all_as_ref Should all genes not in the test set be used as the reference? If FALSE, a random subset is taken of size equal to the test size.
##' @param mean_form_alt Custom two-sided model formula for the (conditional) mean model under the null. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the dependent variable matches the argument to the parameter "count."
##' @param zi_form_alt Custom one-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param mean_form_null Custom two-sided model formula for the (conditional) mean model under the null. Syntax is as in \code{mean_form_alt}.
##' @param zi_form_null Custom one-sided model formula for the zero-inflation model under the null. Syntax is as in \code{zi_form_alt}.
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param lr.df degrees of freedom for the asymptotic chi-square approximation to the liklihood ratio statistic.
##' @param rho Inter-gene correlation value. If NULL (default), estimated using TWO-SIGMA model residuals.
##' @param allow_neg_corr Should negative correlation values be allowed? If FALSE, correlation is set to zero (leads to conservative inference).
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model. Defaults to NULL for constant dispersion.
##' @param return_fits Should complete model fits be returned? Use cautiously because fit objects can be very large.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return An object of class \code{glmmTMB}.
##' @param LR_stats_gene_level_all: Gives all gene-level likelihood ratio statistics.  Order matches the order of the inputted count matrix
##' @param p.vals_gene_level: Gives p-values associated with `LR_stats_gene_level_all`.
##' @param set_p.val: Vector of unadjusted set-level p-values. Order matches the order of inputted test sets.
##' @param corr: Vector of estimated inter-gene correlations for each test set. Order matches the order of inputted test sets.
##' @param test_sets: Vector of numeric indices corresponding to genes in each test set.
##' @param ref_sets: Vector of numeric indices corresponding to the genes in each reference set.
##' @param gene_level_fits: List of gene-level fits, returned only if `return_fits` is TRUE.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm cor pnorm vcov
##' @export twosigmag_custom

twosigmag_custom<-function(count_matrix,index_test,index_ref=NULL,all_as_ref=FALSE,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null
  ,id,lr.df,
  covar_logFC
  ,rho=NULL
  ,allow_neg_corr=FALSE
  ,disp_covar=NULL #need to be able to use data option?
  ,return_fits=FALSE
  ,weights=rep(1,length(count_matrix[1,]))
  ,control = glmmTMBControl()){

  #if(!(adhoc==FALSE)){print("The adhoc method is not recommended for gene set testing due to interpretability.")}

  count_matrix<-as.matrix(count_matrix)
  if(is.list(index_test)){
    nsets<-length(index_test)
    list_lengths<-lapply(index_test,FUN=length)
    if(sum(list_lengths<2)>0){stop("All test sets must have at least two genes. Please remove singleton or empty sets.")}
  }else{ #index.test is a single numeric vector
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
      fit_twosigmag[[l]]<-lr.twosigma_custom(count=count_matrix[l,]
        ,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null,id=id
        ,lr.df = lr.df,covar_logFC=covar_logFC)
      residuals_all[l,]<-residuals(fit_twosigmag[[l]]$fit_alt)
      stats_all[l]<-fit_twosigmag[[l]]$LR_stat
      p.vals_gene_level[l]<-fit_twosigmag[[l]]$LR_p.val
      avg_logFC_gene_level[l]<-fit_twosigmag[[l]]$mean_comp_logFC
    }else{
      fit_twosigmag<-lr.twosigma_custom(count=count_matrix[l,]
        ,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null,id=id
        ,lr.df = lr.df,covar_logFC=covar_logFC)
      residuals_all[l,]<-residuals(fit_twosigmag$fit_alt)
      stats_all[l]<-fit_twosigmag$LR_stat
      p.vals_gene_level[l]<-fit_twosigmag$LR_p.val
      avg_logFC_gene_level[l]<-fit_twosigmag$mean_comp_logFC
    }
    print(paste("Finished Gene Number",i,"of",ngenes))
  }
  stats_test<-vector('list',length=nsets)
  stats_ref<-vector('list',length=nsets)
  direction<-vector(length=nsets)
  p.val<-numeric(length=nsets)
  rho_est<-numeric(length=nsets)
  #browser()
  if(is.null(index_ref)){index_ref<-vector('list',length=nsets)}
  print(paste("Estimating Set-Level correlations and calculating p-values"))
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
    #   #index_ref[[i]]<-setdiff(1:nrow(count_matrix),index_test[[i]])
    #   stats_ref[[i]]<-stats_all[index_ref[[i]]]
    #   ref_size<-length(index_ref[[i]])
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
      #browser()
      for(y in unique_id){
        j<-j+1
        temp<-cor(t(residuals_test[,which(id==y)])) # any checks for id behavior needed here?
        cor_temp[j]<-mean(temp[upper.tri(temp)],na.rm = T)
      }
      rho_est[i]<-mean(cor_temp)
    }
    #browser()
    if(!allow_neg_corr & rho_est[i]<0){rho_est[i]<-0}
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho_est[i])+(test_size-1)*asin((rho_est[i]+1)/2))
    wilcox_stat<-sum(rank(c(stats_test[[i]],stats_ref[[i]]))[1:test_size]) - .5*test_size*(test_size+1)
    p.val[i]<-2*pnorm(-1*abs((wilcox_stat-.5*test_size*ref_size)/sqrt(var)))
  }
  names(p.val)<-names(index_test)
  names(stats_all)<-rownames(count_matrix)
  names(p.vals_gene_level)<-rownames(count_matrix)
  names(avg_logFC_gene_level)<-rownames(count_matrix)
  names(rho_est)<-names(index_test)
  names(direction)<-names(index_test)
  #browser()
  if(return_fits==TRUE){
    return(list(LR_stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,direction=direction,avg_logFC_gene_level=avg_logFC_gene_level,corr=rho_est,test_sets=index_test,ref_sets=index_ref,gene_level_fits=fit_twosigmag))
  }else{
    return(list(LR_stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,direction=direction,avg_logFC_gene_level=avg_logFC_gene_level,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
  }
}
