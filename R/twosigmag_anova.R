##' Gene set testing adjusting for inter-gene correlation using the TWO-SIGMA model of ... with custom user-specified model formulas.
##' @param count_matrix Matrix of non-negative integer read counts. If specifying custom formula(s) via the arguments mean_form and zi_form the expression in mean_form will supersede.
##' @param index_test Index corresponding to rows of the count matrix that are in the test set. Either a list with each element being a different set or a numeric vector to test only a single set.
##' @param index_ref Index corresponding to rows of the count matrix that are in the reference set.  If NULL, a reference set is randomly selected of the same size as test size using genes not in the test set or using all other genes. See all_as_ref.
##' @param all_as_ref Should all genes not in the test set be used as the reference? If FALSE, a random subset is taken of size equal to the test size.
##' @param mean_form Custom two-sided model formula for the (conditional) mean model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the dependent variable matches the argument to the parameter "count."
##' @param zi_form Custom one-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param contrast Either a string indicating the column name of the covariate to test or an integer referring to its column position in BOTH the mean_covar and zi_covar matrices (if the two matrices differ using a string name is preferred). Argument is ignored if mean_covar and zi_covar are both a single covariate (that covariate is assumed of interest).
##' @param id Vector of individual-level ID's. Used for random effect prediction and the adhoc method but required regardless.
##' @param rho Inter-gene correlation value. If NULL (default), estimated using TWO-SIGMA model residuals.
##' @param allow_neg_corr Should negative correlation values be allowed? If FALSE, correlation is set to zero (leads to conservative inference).
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model. Defaults to NULL for constant dispersion.
##' @param return_fits Should complete model fits be returned?
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return An object of class \code{glmmTMB}.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm
##' @importFrom  multcomp glht mcp adjusted
##' @export twosigmag_anova

twosigmag_anova<-function(count_matrix,index_test,index_ref=NULL,all_as_ref=FALSE,mean_form,zi_form
  ,id,contrast=NULL
  ,rho=NULL
  ,allow_neg_corr=FALSE
  ,disp_covar=NULL #need to be able to use data option?
  ,return_fits=FALSE
  ,weights=rep(1,length(count_matrix[1,]))
  ,control = glmmTMBControl()){
  count_matrix<-as.matrix(count_matrix)
  ncomps<-nrow(contrast)
  #browser()
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
  stats_all<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast))
  p.vals_gene_level_raw<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast))
  p.vals_gene_level_adjust<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast))
  estimates_gene_level<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast))
  #browser()
  re_present<-any(grepl("id",mean_form[[3]])>0)
  gene_error<-rep(FALSE,nrow(count_matrix))
  gene_err<-FALSE
  if(re_present){
    re_sigma_est<-rep(NA,length=nrow(count_matrix))
  }
  #browser()
  for(i in 1:ngenes){
    l<-genes[i]
    if(return_fits==TRUE){
      fit_twosigmag[[l]]<-twosigma_custom(count=count_matrix[l,]
        ,mean_form=mean_form,zi_form=zi_form,id=id)
      if(re_present){
        re_sigma_est[l]<-exp(fit_twosigmag$sdr$par.fixed['theta'])
      }
      residuals_all[l,]<-residuals(fit_twosigmag[[l]])
      temp2<-glht_glmmTMB(fit_twosigmag[[l]],
        linfct = mcp(combined_num2_factor = contrast))
      temp<-summary(temp2,test=adjusted("none"))
      stats_all[l,]<-temp$test$tstat
      p.vals_gene_level_raw[l,]<-2*pnorm(-1*abs(temp$test$tstat))
      temp<-summary(temp2)
      estimates_gene_level[l,]<-temp$test$coefficients
      #p.vals_gene_level_adjust[l,]<-as.numeric(temp$test$pvalues)
    }else{
      tryCatch({
      fit_twosigmag<-twosigma_custom(count=count_matrix[l,]
        ,mean_form=mean_form,zi_form=zi_form,id=id)
      if(re_present){
        re_sigma_est[l]<-exp(fit_twosigmag$sdr$par.fixed['theta'])
      }
      residuals_all[l,]<-residuals(fit_twosigmag)
      temp2<-glht_glmmTMB(fit_twosigmag,
          linfct = mcp(combined_num2_factor = contrast))
      temp<-summary(temp2,test=adjusted("none"))
      stats_all[l,]<-temp$test$tstat
      p.vals_gene_level_raw[l,]<-2*pnorm(-1*abs(temp$test$tstat))
      temp<-summary(temp2)
      estimates_gene_level[l,]<-temp$test$coefficients
      #p.vals_gene_level_adjust[l,]<-as.numeric(temp$test$pvalues)
      }
      ,error=function(e){
        #assign("gene_err",TRUE,env=parent.frame())
        print(paste0("Error at Gene ",l, " in the Dataset"))
        #print(gene_err)
        })
    }
    #gene_error[l]<-gene_err
    print(paste("Finished Gene Number",i,"of",ngenes))
    #browser()
  }
  #browser()
  #residuals_all[which(gene_error),]<-rep(NA,ncells)
  #stats_all[which(gene_error),]<-rep(NA,nrow(contrast))
  #p.vals_gene_level_raw[which(gene_error)]<-NA
  #estimates_gene_level[which(gene_error),]<-rep(NA,nrow(contrast))
  #re_sigma_est[which(gene_error)]<-NA


  #browser()
  stats_test<-vector('list',length=nsets)
  stats_ref<-vector('list',length=nsets)
  p.val<-matrix(NA,nrow=nsets,ncol=ncomps)
  rho_est<-numeric(length=nsets)
  #browser()
  #if(is.null(index_ref)){index_ref<-vector('list',length=nsets)}
  for(i in 1:nsets){
    if(is.list(index_test)){
      stats_test_temp<-matrix(NA,nrow=length(index_test[[i]]),ncol=nrow(contrast))
      if(ref_inputted==FALSE){
        stats_ref_temp<-matrix(NA,nrow=length(index_ref[[i]]),ncol=nrow(contrast))
      }else{
        if(is.list(index_ref)){
          stats_ref_temp<-matrix(NA,nrow=length(index_ref[[i]]),ncol=nrow(contrast))
        }else{
          stats_ref_temp<-matrix(NA,nrow=length(index_ref[[i]]),ncol=nrow(contrast))
        }
      }
    }else{ # testing a single gene set
      stats_test_temp<-matrix(NA,nrow=length(index_test),ncol=nrow(contrast))
      if(ref_inputted==FALSE){
        stats_ref_temp<-matrix(NA,nrow=length(index_ref),ncol=nrow(contrast))
      }else{
        if(is.list(index_ref)){
          stats_ref_temp<-matrix(NA,nrow=length(index_ref),ncol=nrow(contrast))
        }else{
          stats_ref_temp<-matrix(NA,nrow=length(index_ref),ncol=nrow(contrast))
        }
      }
    }
    if(is.list(index_test)){
      stats_test_temp<-stats_all[index_test[[i]],]
      stats_test[[i]]<-stats_test_temp
      stats_test[[i]]<-stats_test[[i]][!is.na(stats_test[[i]])]
      test_size<-length(stats_test[[i]])

      residuals_test<-residuals_all[index_test[[i]],]
      stats_ref_temp<-stats_all[index_ref[[i]],]
      stats_ref[[i]]<-stats_ref_temp
      stats_ref[[i]]<-stats_ref[[i]][!is.na(stats_ref[[i]])]
      ref_size<-length(stats_ref[[i]])
    }else{# then both index_test and index_ref should be numeric vectors
      stats_test_temp<-stats_all[index_test,]
      stats_test[[i]]<-stats_test_temp
      stats_test[[i]]<-stats_test[[i]][!is.na(stats_test[[i]])]
      test_size<-length(stats_test[[i]])

      residuals_test<-residuals_all[index_test,]
      stats_ref_temp<-stats_all[index_ref,]
      stats_ref[[i]]<-stats_ref_temp
      stats_ref[[i]]<-stats_ref[[i]][!is.na(stats_ref[[i]])]
      ref_size<-length(stats_ref[[i]])


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
    if(!allow_neg_corr & rho_est[i]<0){rho_est[i]<-0}
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho_est[i])+(test_size-1)*asin((rho_est[i]+1)/2))
    for(b in 1:ncomps){
      wilcox_stat<-sum(rank(c(stats_test[[i]][,b],stats_ref[[i]][,b]))[1:test_size]) - .5*test_size*(test_size+1)
      p.val[i,b]<-2*pnorm(-1*abs((wilcox_stat-.5*test_size*ref_size)/sqrt(var)))
    }
    #need to add code to take test and ref set test statistics by column for test and output set level p-values as a matrix
  }
  #browser()
  if(re_present){
    if(return_fits==TRUE){
      return(list(gene_level_fits=fit_twosigmag,z_stats_gene_level_all=stats_all,set_p.val=p.val,p.vals_gene_level_raw=p.vals_gene_level_raw
        #,p.vals_gene_level_adjust_anova=p.vals_gene_level_adjust
        ,estimates_gene_level=estimates_gene_level
        ,re_sigma_est=re_sigma_est
        ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
    }else{
      return(list(z_stats_gene_level_all=stats_all,set_p.val=p.val,p.vals_gene_level_raw=p.vals_gene_level_raw
        #,p.vals_gene_level_adjust_anova=p.vals_gene_level_adjust
        ,estimates_gene_level=estimates_gene_level
        ,re_sigma_est=re_sigma_est
        ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
    }
  }else{
      if(return_fits==TRUE){
        return(list(gene_level_fits=fit_twosigmag,z_stats_gene_level_all=stats_all,set_p.val=p.val,p.vals_gene_level_raw=p.vals_gene_level_raw
          #,p.vals_gene_level_adjust_anova=p.vals_gene_level_adjust
          ,estimates_gene_level=estimates_gene_level
          ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
      }else{
        return(list(z_stats_gene_level_all=stats_all,set_p.val=p.val,p.vals_gene_level_raw=p.vals_gene_level_raw
          #,p.vals_gene_level_adjust_anova=p.vals_gene_level_adjust
          ,estimates_gene_level=estimates_gene_level
          ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
      }}
}
