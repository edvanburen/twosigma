##' Gene set testing for single-cell RNA-sequencing data adjusting for inter-gene correlation using custom user-specified model formulas.
##' @param count_matrix Matrix of non-negative integer read counts. If specifying custom formula(s) via the arguments mean_form and zi_form the expression in mean_form will supersede. It is recommended to make the rownames the gene names for better output.
##' @param index_test List of indices corresponding to rows of the count matrix that are in the test set. Names of each list element (i.e. Gene Set Names) are carried forward to output if present.
##' @param index_ref List of indices corresponding to rows of the count matrix that are in the reference set.  If NULL, a reference set is randomly selected of the same size as the test size using genes not in the test set (if all_as_ref=FALSE) or using all other genes (if all_as_ref=TRUE). See \code{all_as_ref}. Must be either NULL or a list with the same length as index_test.
##' @param all_as_ref Should all genes not in the test set be used as the reference? If FALSE, a random subset is taken of size equal to the test size.
##' @param mean_form_alt Custom two-sided model formula for the (conditional) mean model under the null. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the LHS of the formula begins with "count."
##' @param zi_form_alt Custom one-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in lme4.
##' @param mean_form_null Custom two-sided model formula for the (conditional) mean model under the null. Syntax is as in \code{mean_form_alt}. Users should ensure that the LHS of the formula begins with "count."
##' @param zi_form_null Custom one-sided model formula for the zero-inflation model under the null. Syntax is as in \code{zi_form_alt}.
##' @param id Vector of individual-level (sample-level) ID's. Used to estimate inter-gene correlation and random effect prediction (if present) and is currently required.
##' @param lr.df degrees of freedom for the asymptotic chi-square approximation to the liklihood ratio statistic.
##' @param mean_covar_logFC Covariate used for reporting direction (as Up or Down) of the test set. Either a string indicating the name of the covariate to use or an integer giving its associated position in the RHS of the mean_form_alt argument.
##' @param rho Inter-gene correlation value. If NULL (default), estimated using TWO-SIGMA model residuals.
##' @param allow_neg_corr Should negative correlation values be allowed? If FALSE, correlation is set to zero (leads to conservative inference).
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model. Defaults to NULL for constant dispersion.
##' @param return_summary_fits If TRUE, returns a summary.glmmTMB object for the Alternative model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.  See \code{?glmmTMBControl}.
##' @section Details:  If adhoc=TRUE, any input in mean_re and zi_re will be ignored.
##' @return A list with the following elements:
##' ##' \itemize{
##' \item{\code{gene_summary_fits: }}{Summary.glmmTMB objects for each gene from the alternative model (if return_summary_fits=TRUE)}
##' \item{\code{LR_stats_gene_level_all: }}{Gives all gene-level likelihood ratio statistics.  Order matches the order of the inputted count matrix}
##' \item{\code{p.vals_gene_level: }}{Gives p-values associated with \code{LR_stats_gene_level_all}.}
##' \item{\code{set_p.val: }}{Vector of unadjusted set-level p-values. Order matches the order of inputted test sets.}
##' \item{\code{set_p.val_ttest: }}{Vector of unadjusted set-level p-values using the t-test. Order matches the order of inputted test sets.}
##' \item{\code{avg_logFC_gene_level: }}{Gives the average logFC for the covariate specified in the mean_covar_logFC argument.}
##' \item{\code{direction: }}{Reports whether the test set tends to be Up or Down Regulated based on the covariate specified in the mean_covar_logFC argument.}
##' \item{\code{corr: }}{Vector of estimated inter-gene correlations for each test set. Order matches the order of inputted test sets.}
##' \item{\code{test_sets: }}{Vector of numeric indices corresponding to genes in each test set.}
##' \item{\code{ref_sets: }}{Vector of numeric indices corresponding to the genes in each reference set.}
##' }
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm cor pnorm vcov pt
##' @export twosigmag_custom

twosigmag_custom<-function(count_matrix,index_test,index_ref=NULL,all_as_ref=FALSE,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null
  ,id,lr.df,mean_covar_logFC
  ,rho=NULL
  ,allow_neg_corr=FALSE
  ,disp_covar=NULL #need to be able to use data option?
  ,return_summary_fits=TRUE
  ,weights=rep(1,ncol(count_matrix))
  ,control = glmmTMBControl()){

  #if(!(adhoc==FALSE)){print("The adhoc method is not recommended for gene set testing due to interpretability.")}

  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","index_test","mean_form_alt","zi_form_alt","mean_form_null","zi_form_null","id"
    ,"lr.df","mean_covar_logFC")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.list(index_test)){stop("Please ensure the input index_test is a list, even if only testing one set.")}
  if(!is.null(index_ref) & !is.list(index_ref)){stop("Please ensure the input index_ref is a list, even if only testing one set.")}
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
  residuals_all<-matrix(nrow=nrow(count_matrix),ncol=ncells)
  stats_all<-rep(NA,length=nrow(count_matrix))
  p.vals_gene_level<-rep(NA,length=nrow(count_matrix))
  avg_logFC_gene_level<-rep(NA,length=nrow(count_matrix))
  for(i in 1:ngenes){
    l<-genes[i]
      fit_twosigmag<-lr.twosigma_custom(count_matrix[l,,drop=FALSE],silent=TRUE
        ,mean_form_alt,zi_form_alt,mean_form_null,zi_form_null,id=id,return_full_fits = TRUE
        ,lr.df = lr.df)
      residuals_all[l,]<-residuals(fit_twosigmag$full_fit_alt[[1]])
      stats_all[l]<-fit_twosigmag$LR_stat[1]
      p.vals_gene_level[l]<-fit_twosigmag$LR_p.val[1]
      if(return_summary_fits==TRUE){
        fit[[l]]<-fit_twosigmag$summary_fit_alt[[1]]
      }
      sum_fit_alt<-fit_twosigmag$summary_fit_alt[[1]]$coefficients$cond
      names<-rownames(sum_fit_alt)
      if(is.numeric(mean_covar_logFC)){
        #+1 for intercept
        avg_logFC_gene_level[l]<-sum_fit_alt[mean_covar_logFC+1,'Estimate']
      }
      if(is.character(mean_covar_logFC)){
        avg_logFC_gene_level[l]<-sum_fit_alt[grepl(mean_covar_logFC,names),'Estimate']
      }
    print(paste("Finished Gene Number",i,"of",ngenes))
  }
  stats_test<-vector('list',length=nsets)
  stats_ref<-vector('list',length=nsets)
  direction<-vector(length=nsets)
  p.val<-numeric(length=nsets)
  p.val_ttest<-numeric(length=nsets)
  rho_est<-numeric(length=nsets)
  if(is.null(index_ref)){index_ref<-vector('list',length=nsets)}
  for(i in 1:nsets){
      stats_test[[i]]<-stats_all[index_test[[i]]]
      residuals_test<-residuals_all[index_test[[i]],]
      direction[i]<-ifelse(sign(mean(avg_logFC_gene_level[index_test[[i]]]))==1,"Up","Down")
      stats_ref[[i]]<-stats_all[index_ref[[i]]]
    # }

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
    # Missing values will get dropped in rank statement
    # so make sure sizes of test and reference sets don't include them
    # Don't want to drop missing values because want them documented in output
    test_size<-length(stats_test[[i]][!is.na(stats_test[[i]])])
    ref_size<-length(stats_ref[[i]][!is.na(stats_ref[[i]])])
    var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*rho_est[i])+(test_size-1)*asin((rho_est[i]+1)/2))
    wilcox_stat<-sum(rank(c(stats_test[[i]],stats_ref[[i]]),na.last = NA)[1:test_size]) - .5*test_size*(test_size+1)
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

  #browser()
  if(return_summary_fits==TRUE){
    names(fit)<-rownames(count_matrix)
    return(list(gene_summary_fits=fit,LR_stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,set_p.val_ttest=p.val_ttest,avg_logFC_gene_level=avg_logFC_gene_level,direction=direction,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
  }else{
    return(list(LR_stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,set_p.val_ttest=p.val_ttest,avg_logFC_gene_level=avg_logFC_gene_level,direction=direction,corr=rho_est,test_sets=index_test,ref_sets=index_ref))
  }
}
