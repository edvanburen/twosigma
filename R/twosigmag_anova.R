##' Gene set testing of ANOVA-style pairwise contrasts for single-cell RNA-sequencing data adjusting for inter-gene correlation using custom user-specified model formulas.
##' @param count_matrix Matrix of non-negative integer read counts, with rows corresponding to genes and columns correspoding to cells. It is recommended to make the rownames the gene names for better output.
##' @param index_test List of indices corresponding to rows of the count matrix that are in the test set. Names of each list element (i.e. Gene Set Names) are carried forward to output if present.
##' @param index_ref List of indices corresponding to rows of the count matrix that are in the reference set.  If NULL, a reference set is randomly selected of the same size as the test size using genes not in the test set (if all_as_ref=FALSE) or using all other genes (if all_as_ref=TRUE). See \code{all_as_ref}. Must be either NULL or a list with the same length as index_test.
##' @param all_as_ref Should all genes not in the test set be used as the reference? If FALSE, a random subset is taken of size equal to the test size.
##' @param mean_form Custom two-sided model formula for the (conditional) mean model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the LHS of the formula begins with "count."
##' @param zi_form Custom one-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param id Vector of individual-level (sample-level) ID's. Used to estimate inter-gene correlation and for random effect prediction (if present) and required.
##' @param contrast_matrix Contrast matrix to test using the glht and mcp functions. Rows correspond to different contrasts, and columns refer to the ordered levels of \code{fact_name}. Rownames are carried forward to output if present.
##' @param fact_name Factor variable to test using \code{contrast_matrix}.
##' @param rho Inter-gene correlation value. If NULL (default), estimated using TWO-SIGMA model residuals.
##' @param allow_neg_corr Should negative correlation values be allowed? If FALSE, correlation is set to zero (leads to conservative inference).
##' @param return_summary_fits Should complete model fits be returned?
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model. Defaults to NULL for constant dispersion.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.  See \code{?glmmTMBControl}.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return A list with the following elements:
##' \itemize{
##' \item{\code{gene_summary_fits: }}{Summary.glmmTMB objects for each gene from the alternative model (if return_summary_fits=TRUE)}
##' \item{\code{z_stats_gene_level_all: }}{Matrix of z statistics testing the contrast_matrix specified in \code{contrast_matrix}. Row order matches the order of inputted test sets, and column order matches the order of \code{contrast_matrix}.}
##' \item{\code{p.vals_gene_level: }}{Matrix of p-values testing the contrast_matrix specified in \code{contrast_matrix}. Row order matches the order of inputted test sets, and column order matches the order of \code{contrast_matrix}.}
##' \item{\code{set_p.val: }}{Matrix of unadjusted set-level p-values. Row order matches the order of inputted test sets, and column order matches the order of \code{contrast_matrix}.}
##' \item{\code{set_p.val_ttest: }}{Matrix of unadjusted set-level p-values using the t-test. Row order matches the order of inputted test sets, and column order matches the order of \code{contrast_matrix}.}
##' \item{\code{estimates_gene_level: }}{Matrix of gene-level contrast_matrix estimates. Row order matches the order of inputted test sets, and column order matches the order of \code{contrast_matrix}}
##' \item{\code{direction: }}{Reports whether the test set tends to be Up or Down Regulated based on the covariate specified in the mean_covar_logFC argument.}
##' \item{\code{corr: }}{Vector of estimated inter-gene correlations for each test set. Order matches the order of inputted test sets.}
##' \item{\code{test_sets: }}{Vector of numeric indices corresponding to genes in each test set.}
##' \item{\code{ref_sets: }}{Vector of numeric indices corresponding to the genes in each reference set.}
##' }
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm pt
##' @importFrom  multcomp glht mcp adjusted
##' @export twosigmag_anova

twosigmag_anova<-function(count_matrix,index_test,index_ref=NULL,all_as_ref=FALSE,mean_form,zi_form
  ,id,contrast_matrix,fact_name
  ,rho=NULL
  ,allow_neg_corr=FALSE
  ,return_summary_fits=TRUE
  ,disp_covar=NULL #need to be able to use data option?
  ,weights=rep(1,ncol(count_matrix))
  ,control = glmmTMBControl()){
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","index_test","mean_form","zi_form","id","contrast_matrix","fact_name")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.matrix(count_matrix)){stop("Please ensure the input count_matrix is of class matrix.")}
  if(!is.matrix(contrast_matrix)){stop("Please ensure the input contrast_matrix is of class matrix.")}
  if(!is.list(index_test)){stop("Please ensure the input index_test is a list, even if only testing one set.")}
  if(!is.null(index_ref) & !is.list(index_ref)){stop("Please ensure the input index_ref is a list, even if only testing one set.")}
  if(length(id)!=ncol(count_matrix)){stop("Argument id should be a numeric vector with length equal to the number of columns of count_matrix (i.e. the number of cells).")}
  ncomps<-nrow(contrast_matrix)
  ncells<-ncol(count_matrix)
  nsets<-length(index_test)
  list_lengths<-lapply(index_test,FUN=length)
  if(sum(list_lengths<2)>0){stop("All test sets must have at least two genes. Please remove singleton or empty sets.")}
  if(all_as_ref==TRUE & !is.null(index_ref)){stop("Please specify either all_as_ref=TRUE or index_ref as a non-NULL input. If all_as_ref is TRUE, then index_ref must be NULL.")}

  if(!is.null(index_ref)){
    for(i in 1:nsets){
      if(sum(index_test[[i]]%in%index_ref[[i]])>0){stop(paste("A gene should not be in both the test and reference sets. Check element number",i,"in test set or reference set."))}}
    if(is.list(index_test) & length(index_ref)!=nsets){
      stop("index_test and index_ref should be lists of the same length.")}
    genes<-unique(c(unlist(index_test),unlist(index_ref)))
    ngenes<-length(genes)
    ref_inputted<-TRUE
  }else {# will need to construct reference set
    index_ref<-vector('list',length=nsets)
    for(i in 1:nsets){
      if(all_as_ref==FALSE){
        index_ref[[i]]<-sample(setdiff(1:nrow(count_matrix),index_test[[i]]),size=length(index_test[[i]]))
      }else{#all_as_ref==TRUE
        index_ref[[i]]<-setdiff(1:nrow(count_matrix),index_test[[i]])
      }
    }
    genes<-unique(c(unlist(index_test),unlist(index_ref)))
    ngenes<-length(genes)
    ref_inputted<-FALSE
  }

  if(max(unlist(index_test))>nrow(count_matrix) | min(unlist(index_test))<1){stop("index_test seems to be invalid, indices must be numeric within the row dimensions of the input count_matrix")}

  # Fit all gene level statistics that are needed
  if(return_summary_fits==TRUE){
    fit<-vector('list',length=nrow(count_matrix))
  }
  zi_stat<-rep(NA,length=nrow(count_matrix))
  residuals_all<-matrix(nrow=nrow(count_matrix),ncol=ncells)
  stats_all<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast_matrix))
  p.vals_gene_level_raw<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast_matrix))
  p.vals_gene_level_adjust<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast_matrix))
  estimates_gene_level<-matrix(NA,nrow=nrow(count_matrix),ncol=nrow(contrast_matrix))
  re_present<-any(grepl("id",mean_form[[3]])>0)
  gene_error<-rep(FALSE,nrow(count_matrix))
  gene_err<-FALSE
  if(re_present){
    re_sigma_est<-rep(NA,length=nrow(count_matrix))
  }
  for(i in 1:ngenes){
    #browser()
    l<-genes[i]
      tryCatch({
      fit_twosigmag<-twosigma_custom(count_matrix[l,,drop=FALSE]
        ,mean_form=mean_form,zi_form=zi_form,id=id,return_summary_fits = FALSE,silent=TRUE)
      if(return_summary_fits==TRUE){
        fit[[l]]<-summary(fit_twosigmag[[1]])
      }
      if(re_present){
        re_sigma_est[l]<-exp(fit_twosigmag[[1]]$sdr$par.fixed['theta'])
      }
      residuals_all[l,]<-residuals(fit_twosigmag[[1]])
      mcp<- mcp(fact_name = contrast_matrix)
      names(mcp)<-fact_name
      temp2<-glht_glmmTMB(fit_twosigmag[[1]],
          linfct = mcp)
      temp<-summary(temp2,test=adjusted("none"))
      stats_all[l,]<-temp$test$tstat
      p.vals_gene_level_raw[l,]<-2*pnorm(-1*abs(temp$test$tstat))
      temp<-summary(temp2)
      estimates_gene_level[l,]<-temp$test$coefficients
      }
      ,error=function(e){
        #assign("gene_err",TRUE,env=parent.frame())
        print(paste0("Error at Gene ",l, " in the Dataset:",e))
        #print(gene_err)
        })

    #gene_error[l]<-gene_err
    print(paste("Finished Gene Number",i,"of",ngenes))
  }
  #residuals_all[which(gene_error),]<-rep(NA,ncells)
  #stats_all[which(gene_error),]<-rep(NA,nrow(contrast_matrix))
  #p.vals_gene_level_raw[which(gene_error)]<-NA
  #estimates_gene_level[which(gene_error),]<-rep(NA,nrow(contrast_matrix))
  #re_sigma_est[which(gene_error)]<-NA


  stats_test<-vector('list',length=nsets)
  stats_ref<-vector('list',length=nsets)
  p.val<-matrix(NA,nrow=nsets,ncol=ncomps)
  p.val_ttest<-matrix(NA,nrow=nsets,ncol=ncomps)
  rho_est<-numeric(length=nsets)
  direction<-matrix(NA,nrow=nsets,ncol=ncomps)
  for(i in 1:nsets){
      stats_test[[i]]<-stats_all[index_test[[i]],,drop=FALSE]
      residuals_test<-residuals_all[index_test[[i]],]
      stats_ref[[i]]<-stats_all[index_ref[[i]],,drop=FALSE]

    if(!is.null(rho)){rho_est[i]<-rho}
    if(is.null(rho)){
      nind<-length(unique(id))
      cor_temp<-numeric(length=nind)
      unique_id<-unique(id)
      j<-0
      for(y in unique_id){
        j<-j+1
        temp<-cor(t(residuals_test[,which(id==y)])) # any checks for id behavior needed here?
        cor_temp[j]<-mean(temp[upper.tri(temp)],na.rm = T)
      }
      rho_est[i]<-mean(cor_temp,na.rm=T)
    }
    if(!allow_neg_corr & rho_est[i]<0){rho_est[i]<-0}
    for(b in 1:ncomps){
      direction[i,b]<-ifelse(sign(mean(estimates_gene_level[index_test[[i]],b],na.rm=T))==1,"Up","Down")
      # Missing values will get dropped in rank statement
      # so make sure sizes of test and reference sets don't include them
      # Don't want to drop missing values because want them documented in output
      test_size<-length(stats_test[[i]][,b][!is.na(stats_test[[i]][,b])])
      ref_size<-length(stats_ref[[i]][,b][!is.na(stats_ref[[i]][,b])])
      var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho_est[i])+(test_size-1)*asin((rho_est[i]+1)/2))
      wilcox_stat<-sum(rank(c(stats_test[[i]][,b],stats_ref[[i]][,b]),na.last = NA)[1:test_size]) - .5*test_size*(test_size+1)
      p.val[i,b]<-2*pnorm(-1*abs((wilcox_stat-.5*test_size*ref_size)/sqrt(var)))

      # T-test output below

      n_genes_tested<-test_size+ref_size
      delta<-(n_genes_tested)/ref_size*(mean(stats_test[[i]][,b],na.rm=T)-mean(c(stats_test[[i]][,b],stats_ref[[i]][,b]),na.rm=T))
      vif<-1+(test_size-1)*rho_est[i]
      varStatPooled<-((n_genes_tested-1)*var(c(stats_test[[i]][,b],stats_ref[[i]][,b]),na.rm=T)-delta^2*test_size*ref_size/n_genes_tested)/(n_genes_tested-2)
      two.sample.t <- delta / sqrt( varStatPooled * (vif/test_size + 1/ref_size) )
      p.val_ttest[i,b]<-2*pt(-1*abs(two.sample.t),df=n_genes_tested-2)
      # This only happens if all genes in the test set are NA
      # In this case p-value will report as zero when it's really NA
      if(test_size==0 | ref_size==0){p.val[i,b]<-NA;p.val_ttest[i,b]<-NA}
    }
    #browser()
    if(i%%100==0){print(paste0("Set ",i," of ",nsets," Finished"))}
  }
  colnames(p.val)<-rownames(contrast_matrix)
  rownames(p.val)<-names(index_test)
  colnames(p.val_ttest)<-rownames(contrast_matrix)
  rownames(p.val_ttest)<-names(index_test)
  rownames(stats_all)<-rownames(count_matrix)
  colnames(stats_all)<-rownames(contrast_matrix)
  rownames(p.vals_gene_level_raw)<-rownames(count_matrix)
  colnames(p.vals_gene_level_raw)<-rownames(contrast_matrix)
  rownames(estimates_gene_level)<-rownames(count_matrix)
  colnames(estimates_gene_level)<-rownames(contrast_matrix)
  rownames(direction)<-names(index_test)
  colnames(direction)<-rownames(contrast_matrix)
  names(rho_est)<-names(index_test)
  if(re_present){
    if(return_summary_fits==TRUE){
      names(fit)<-rownames(count_matrix)
      return(list(summary_gene_level_fits=fit,z_stats_gene_level_all=stats_all,set_p.val=p.val,p.vals_gene_level_raw=p.vals_gene_level_raw
        ,estimates_gene_level=estimates_gene_level,zi_stat=zi_stat
        ,re_sigma_est=re_sigma_est,set_p.val_ttest=p.val_ttest
        ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
    }else{
      return(list(z_stats_gene_level_all=stats_all,set_p.val=p.val,p.vals_gene_level_raw=p.vals_gene_level_raw
        ,estimates_gene_level=estimates_gene_level,zi_stat=zi_stat
        ,re_sigma_est=re_sigma_est,set_p.val_ttest=p.val_ttest
        ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
    }
  }else{
      if(return_summary_fits==TRUE){
        names(fit)<-rownames(count_matrix)
        return(list(summary_gene_level_fits=fit,z_stats_gene_level_all=stats_all,set_p.val=p.val,set_p.val_ttest=p.val_ttest,p.vals_gene_level_raw=p.vals_gene_level_raw
          ,estimates_gene_level=estimates_gene_level,direction=direction,zi_stat=zi_stat,set_p.val_ttest=p.val_ttest
          ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
      }else{
        return(list(z_stats_gene_level_all=stats_all,set_p.val=p.val,set_p.val_ttest=p.val_ttest,p.vals_gene_level_raw=p.vals_gene_level_raw
          ,estimates_gene_level=estimates_gene_level,direction=direction,zi_stat=zi_stat,set_p.val_ttest=p.val_ttest
          ,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
      }}
}
