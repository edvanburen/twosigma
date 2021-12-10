##' Gene set testing for single-cell RNA-sequencing data adjusting for inter-gene correlation.
##' @param count_matrix Matrix of non-negative integer read counts. It is recommended to make the rownames the gene names for better output. No missing values can be present in the data.
##' @param index_test List of indices corresponding to rows of the count matrix that are in the test set. Names of each list element (i.e. Gene Set Names) are carried forward to output if present.
##' @param index_ref List of indices corresponding to rows of the count matrix that are in the reference set.  If \code{NULL}, a reference set is randomly selected of the same size as the test size using genes not in the test set (if \code{all_as_ref=FALSE}) or using all other genes (if \code{all_as_ref=TRUE}). See \code{all_as_ref}. Must be either \code{NULL} or a list with the same length as \code{index_test}.
##' @param all_as_ref Should all genes not in the test set be used as the reference? If \code{FALSE}, a random subset is taken of size equal to the test size.
##' @param mean_form Two-sided model formula for the (conditional) mean model. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package. Users should ensure that the LHS of the formula contains  '\code{count} '.
##' @param zi_form One-sided model formula for the zero-inflation model under the alternative. Formula is passed directly into glmmTMB with random effects specified as in the lme4 package.
##' @param mean_form_null Two-sided model formula for the (conditional) mean model under the null. Needed if and only if \code{statistic='LR'}. Syntax is as in \code{mean_form}. Users should ensure that the LHS of the formula contains '\code{count} '.
##' @param zi_form_null One-sided model formula for the zero-inflation model under the null. Needed if and only if \code{statistic='LR'}. Syntax is as in \code{zi_form}.
##' @param id Vector of individual-level (sample-level) ID's. Used to estimate inter-gene correlation and random effect prediction (if present) and is currently required.
##' @param statistic Which gene-level statistic should be used. Options are Likelihood Ratio ("LR", default), Z-statistic from the mean model ("Z"),the  Stouffer's method combined Z-statistic ("Stouffer"), or a contrast of regression parameters ("contrast"). If "Stouffer", covar_to_test must be in both components. If "contrast", covar_to_test is not used and must be \code{NULL}.
##' @param lr.df degrees of freedom for the asymptotic chi-square approximation to the likelihood ratio statistic. Needed if and only if \code{statistic='LR'}.
##' @param covar_to_test Covariate used for reporting direction (as Up or Down) of the test set and for collecting gene-level statistics. Either a string indicating the name of the covariate to use or an integer giving its associated position in the RHS of the mean_form argument. If a string, the name is matched to the predictors of the mean model, so users should ensure such a match would be unique. Not required and should be \code{NULL} if \code{statistic='contrast'.}
##' @param contrast_matrix Matrix of contrasts of regression parameters from the mean model to be tested. Each row will have separate gene-level and set-level statistics.  Rownames of \code{contrast_matrix} should correspond to a meaningful name of the hypothesis for nicely formatted output. If testing a factor, must have a number of columns exactly equal to the number of levels of the factor.  Otherwise, must have one column per parameter in the mean model (including a column for the intercept.)
##' @param factor_name Name of the factor being tested by \code{contrast_matrix}. Needed if and only if \code{statistic='contrast'} and \code{contrast_matrix} is testing a factor variable in the mean model.
##' @param rho Inter-gene correlation value. If \code{NULL} (default), estimated using TWO-SIGMA model residuals.
##' @param allow_neg_corr Should negative correlation values be allowed? If FALSE, negative correlations are set to zero (leads to conservative inference)..
##' @param return_summary_fits If \code{TRUE}, returns a list containing objects of class \code{summary.glmmTMB} for each gene.
##' @param weights weights, as in \code{glm}. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.  See \code{?glmmTMBControl}.
##' @param ncores Number of cores used for parallelization. Defaults to 1, meaning no parallelization of any kind is done.
##' @param cluster_type Whether to use a "cluster of type "Fork" or "Sock". On Unix systems, "Fork" will likely improve performance. On Windows, only "Sock" will actually result in parallelized computing.
##' @param chunk_size Number of genes to be sent to each parallel environment. Parallelization is more efficient, particularly with a large count matrix, when the count matrix is 'chunked' into some common size (e.g. 10, 50, 200). Defaults to 10.
##' @param lb Should load balancing be used for parallelization? Users will likely want to set to FALSE for improved performance.
##' @return A list with the following elements:
##' ##' \itemize{
##' \item{\code{stats_gene_level_all: }}{Gives all gene-level statistics.  Order matches the order of the inputted count matrix.}
##' \item{\code{p.vals_gene_level: }}{Gives raw (unadjusted) p-values associated with \code{stats_gene_level_all}.}
##' \item{\code{set_p.val: }}{Unadjusted set-level p-values. Order matches the order of inputted test sets.}
##' \item{\code{set_p.val_FDR: }}{FDR-corrected (using the Benjamini-Hochberg procedure) set-level p-values. Order matches the order of inputted test sets.}
##' \item{\code{estimates_gene_level: }}{Gives the average logFC or contrast estimate for each gene.}
##' \item{\code{se_gene_level: }}{Standard error of the gene-level logFC values. Useful to construct gene-level summary statistics.}
##'\item{\code{estimates_set_level: }}{Gives the set-level average of the gene-level logFC or contrast estimates.}
##' \item{\code{direction: }}{Reports whether the test set tends to be Up or Down Regulated based on the sign of \code{estimates_set_level}.}
##' \item{\code{corr: }}{Vector of estimated inter-gene correlations for each test set. Order matches the order of inputted test sets.}
##' \item{\code{gene_level_loglik: }}{Vector of log-likelihood values for each gene. Values of NA indicates a model fitting or convergence problem for that gene.}
##' \item{\code{gene_error:}} Vector indicating whether the particular gene produced an error during model fitting (TRUE) or not (FALSE).
##' \item{\code{test_sets: }}{Vector of numeric indices corresponding to genes in each test set.}
##' \item{\code{ref_sets: }}{Vector of numeric indices corresponding to the genes in each reference set.}
##' \item{\code{gene_summary_fits: }}{Summary.glmmTMB objects for each gene from the alternative model (if return_summary_fits=TRUE)}
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
##' # Simulate Data, half under null half under alternative
##'
##'sim_dat<-matrix(nrow=4,ncol=sum(ncellsper))
##'for(i in 1:nrow(sim_dat)){
##'  if(i<2){# Gene Sets Under the Null
##'    sim_dat[i,]<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta2
##'    ,phi,sigma.a,sigma.b,id.levels=NULL)$Y
                                                               ##'  }else{# Gene Sets Under the Alternative
                                                               ##'    sim_dat[i,]<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta
                                                               ##'    ,phi,sigma.a,sigma.b,id.levels=NULL)$Y
                                                               ##'  }
                                                               ##'}
##'rownames(sim_dat)<-paste("Gene",1:4)
##'
##' # Run twosigmag
##'
##' twosigmag(sim_dat,index_test = list(c(1,3)),all_as_ref = TRUE,mean_form = count~X
##' ,zi_form = ~0,id=id,covar_to_test  = "t2d_sim",statistic = "Z")
##' @export twosigmag

twosigmag<-function(count_matrix,index_test,index_ref=NULL,all_as_ref=FALSE,mean_form,zi_form,mean_form_null=NULL,zi_form_null=NULL
  ,id,statistic,lr.df=NULL,covar_to_test=NULL
  ,contrast_matrix=NULL,factor_name=NULL,rho=NULL
  ,allow_neg_corr=FALSE
  ,return_summary_fits=FALSE
  ,weights=NULL
  ,control = glmmTMBControl(),ncores=1,cluster_type="Fork",chunk_size=10,lb=FALSE){

  #if(!(adhoc==FALSE)){print("The adhoc method is not recommended for gene set testing due to interpretability.")}
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("count_matrix","index_test","mean_form","zi_form","statistic","id")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(!is.list(index_test)){stop("Please ensure the input index_test is a list, even if only testing one set.")}
  if(!is.null(index_ref) & !is.list(index_ref)){stop("Please ensure the input index_ref is a list, even if only testing one set.")}
  #browser()
  if(statistic=="LR" & (sum(is.null(mean_form_null),is.null(zi_form_null),is.null(lr.df))>0)){stop("If statistic=LR, then arguments mean_form_null, zi_form_null, and lr.df must be specified.")}
  if(statistic!="contrast"& is.null(covar_to_test)){stop("Argument covar_to_test is required and cannot be NULL if statistic = 'Z', 'LR', or 'Stouffer'.")}
  if(statistic=="contrast"& !is.null(covar_to_test)){stop("Argument covar_to_test is not required and must be NULL if statistic = 'contrast'. Gene-level estimates and set-level directions will come from the contast matrix specified in argument 'contrast_matrix'.")}
  if(statistic!="contrast"& (!is.null(contrast_matrix) |!is.null(factor_name))){stop("Arguments 'contrast_matrix' and 'factor_name' are not required and should be NULL unless statistic = 'contrast'.")}
  if(statistic=="contrast"& is.null(contrast_matrix)){stop("Argument contrast_matrix is required if statistic = 'contrast'.")}
  if(statistic!="LR"& (!is.null(mean_form_null) | !is.null(zi_form_null) | !is.null(lr.df))){stop("Arguments mean_form_null, zi_form_null, and lr.df are not needed and should be NULL unless statistic = 'LR'.")}
  if(mean_form[[2]]!="count"){stop("Please begin the two-sided formula mean_form with the name 'count'. Failure to do so will cause an error downstream.")}
  if(length(zi_form)>2){stop("ZI Formula should be one-sided (e.g.  '~1 ') and thus only have a length of 2.")}
  #browser()
  ngenes_total<-nrow(count_matrix)
  gene_names<-rownames(count_matrix)
  genes<-1:nrow(count_matrix)
  ncomps<-ifelse(statistic=="contrast",nrow(contrast_matrix),1)
  ncells<-ncol(count_matrix)
  if(!"weights"%in%passed_args){weights<-rep(1,ncells)}
  nsets<-length(index_test)
  list_lengths<-lapply(index_test,FUN=length)
  if(sum(list_lengths<2)>0){stop("All test sets must have at least two genes. Please remove singleton or empty sets.")}
  if(all_as_ref==TRUE & !is.null(index_ref)){stop("Please specify either all_as_ref=TRUE or index_ref as a non-NULL input. If all_as_ref is TRUE, then index_ref must be NULL.")}
  if(!is.null(index_ref)){
    if(length(index_ref)!=nsets){
      stop("index_test and index_ref should be lists of the same length.")}
    for(i in 1:nsets){
      if(sum(index_test[[i]]%in%index_ref[[i]])>0){stop(paste("A gene should not be in both the test and reference sets. Check element number",i,"in test set or reference set."))}}
    ref_inputted<-TRUE
  }else {# will need to construct reference set
    index_ref<-vector('list',length=nsets)
    for(i in 1:nsets){
      if(all_as_ref==FALSE){
        index_ref[[i]]<-sample(setdiff(1:ngenes_total,index_test[[i]]),size=length(index_test[[i]]))
      }else{#all_as_ref==TRUE
        index_ref[[i]]<-setdiff(1:ngenes_total,index_test[[i]])
      }
    }
    ref_inputted<-FALSE
  }

  if(max(unlist(index_test))>ngenes_total | min(unlist(index_test))<1){stop("index_test seems to be invalid, indices must be numeric within the row dimensions of the input count_matrix")}
  #options(warning.length = 2000L)
  #browser()
  fit_tsg<-function(chunk,statistic,covar_to_test=NULL
    ,factor_name=NULL,contrast_matrix=NULL,id,ncomps){
    gc()
    k<-0
    num_err=0
    residuals_all<-matrix(NA,nrow=length(chunk),ncol=ncells)
    stats_all<-matrix(NA,nrow=length(chunk),ncol=ncomps)
    p.vals_gene_level<-matrix(NA,nrow=length(chunk),ncol=ncomps)
    estimates_gene_level<-matrix(NA,nrow=length(chunk),ncol=ncomps)
    se_gene_level<-matrix(NA,nrow=length(chunk),ncol=ncomps)
    logLik<-numeric(length=length(chunk))
    gene_err<-rep(NA,length(chunk))
    if(return_summary_fits==TRUE){
      fits<-vector('list',length(chunk))
    }
    for(l in unlist(chunk)){
      #print(l)
      k<-k+1
      #setTxtProgressBar(pb, i)
      counts<-count_matrix[l,,drop=FALSE]
      #counts<-count_matrix[l,,drop=FALSE]
      if(num_err>0){break}
      if(statistic=="LR"){
        fit_twosigmag<-lr.twosigma_custom(counts
          ,mean_form_alt=mean_form,zi_form_alt=zi_form
          ,mean_form_null=mean_form_null,zi_form_null=zi_form_null
          ,id=id,return_full_fits = TRUE
          ,lr.df = lr.df,weights=weights,internal_call=TRUE,chunk_size = 1)
        residuals_all[k,]<-residuals(fit_twosigmag$fit_alt[[1]])
        stats_all[k,]<-fit_twosigmag$LR_stat[1]
        p.vals_gene_level[k,]<-fit_twosigmag$LR_p.val[1]
        fit<-summary(fit_twosigmag$fit_alt[[1]])
        if(return_summary_fits==TRUE){
          fits[[k]]<-fit
        }
        logLik[k]<-as.numeric(fit$logLik)
        gene_err[k]<-(is.na(logLik[k]) | fit_twosigmag$fit_alt[[1]]$sdr$pdHess==FALSE)
        sum_fit_alt<-summary(fit_twosigmag$fit_alt[[1]])$coefficients$cond
        names<-rownames(fit$coefficients$cond)
        if(sum(grepl("Intercept",names))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the mean model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the mean model are:",names),collapse=" "))}
        names_zi<-rownames(fit$coefficients$zi)
        if(!is.null(names_zi)&sum(grepl("Intercept",names_zi))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the ZI model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the ZI model are:",names_zi),collapse=" "))}
        if(is.numeric(covar_to_test)){
          #+1 for intercept
          estimates_gene_level[k,]<-sum_fit_alt[covar_to_test+1,'Estimate']
          se_gene_level[k,]<-sum_fit_alt[covar_to_test+1,'Std. Error']
        }
        if(is.character(covar_to_test)){
          if(sum(grepl(covar_to_test,names))>1){num_err=1;stop(paste(c("covar_to_test matches to multiple variables in mean model. Please rename these other variables to not contain the name of covar_to_test. Variable names in the mean model are:",names)),collapse=" ")}
          if(sum(grepl(covar_to_test,names))==0){num_err=1;stop(paste(c("covar_to_test not found in mean model. Variable names in the mean model are:",names),collapse=" "))}
          estimates_gene_level[k,]<-sum_fit_alt[grepl(covar_to_test,names),'Estimate']
          se_gene_level[k,]<-sum_fit_alt[grepl(covar_to_test,names),'Std. Error']
        }
      }
      if(statistic=="Z"){
        fit_twosigmag<-twosigma_custom(counts,silent=TRUE
          ,mean_form=mean_form,zi_form=zi_form
          ,id=id,return_summary_fits = FALSE,weights=weights,internal_call=TRUE,chunk_size=1)
        fit<-summary(fit_twosigmag$fit[[1]])
        if(return_summary_fits==TRUE){
          fits[[k]]<-fit
        }
        logLik[k]<-as.numeric(fit$logLik)
        gene_err[k]<-(is.na(logLik[k]) | fit_twosigmag$fit[[1]]$sdr$pdHess==FALSE)
        residuals_all[k,]<-residuals(fit_twosigmag$fit[[1]])
        names<-rownames(fit$coefficients$cond)
        if(sum(grepl("Intercept",names))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the mean model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the mean model are:",names),collapse=" "))}
        names_zi<-rownames(fit$coefficients$zi)
        if(!is.null(names_zi)&sum(grepl("Intercept",names_zi))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the ZI model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the ZI model are:",names_zi),collapse=" "))}
        if(is.numeric(covar_to_test)){
          #+1 for intercept
          estimates_gene_level[k,]<-fit$coefficients$cond[covar_to_test+1,'Estimate']
          se_gene_level[k,]<-fit$coefficients$cond[covar_to_test+1,'Std. Error']
          stats_all[k,]<-fit$coefficients$cond[covar_to_test,'z value']
          p.vals_gene_level[k,]<-fit$coefficients$cond[covar_to_test,'Pr(>|z|)']
        }
        if(is.character(covar_to_test)){
          if(sum(grepl(covar_to_test,names))>1){num_err=1;stop(paste(c("covar_to_test matches to multiple variables in mean model. Please rename these other variables to not contain the name of covar_to_test. Variable names in the mean model are:",names),collapse=" "))}
          if(sum(grepl(covar_to_test,names))==0){num_err=1;stop(paste(c("covar_to_test not found in mean model. Variable names in the mean model are:",names),collapse=" "))}
          estimates_gene_level[k,]<-fit$coefficients$cond[grepl(covar_to_test,names),'Estimate']
          se_gene_level[k,]<-fit$coefficients$cond[grepl(covar_to_test,names),'Std. Error']
          stats_all[k,]<-fit$coefficients$cond[grepl(covar_to_test,names),'z value']
          p.vals_gene_level[k,]<-fit$coefficients$cond[grepl(covar_to_test,names),'Pr(>|z|)']
        }
      }
      if(statistic=="Stouffer"){
        fit_twosigmag<-twosigma_custom(counts,silent=TRUE
          ,mean_form=mean_form,zi_form=zi_form
          ,id=id,return_summary_fits = FALSE,weights=weights,internal_call=TRUE,chunk_size=1)
        fit<-summary(fit_twosigmag$fit[[1]])
        if(return_summary_fits==TRUE){
          fits[[k]]<-fit
        }
        logLik[k]<-as.numeric(fit$logLik)
        gene_err[k]<-(is.na(logLik[k]) | fit_twosigmag$fit[[1]]$sdr$pdHess==FALSE)
        residuals_all[k,]<-residuals(fit_twosigmag$fit[[1]])
        names_cond<-rownames(fit$coefficients$cond)
        if(sum(grepl("Intercept",names_cond))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the mean model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the mean model are:",names_cond),collapse=" "))}
        names_zi<-rownames(fit$coefficients$zi)
        if(is.null(names_zi)){num_err=1;stop("If statistic='Stouffer', the ZI component must exist.")}
        if(!is.null(names_zi)&sum(grepl("Intercept",names_zi))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the ZI model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the ZI model are:",names_zi),collapse=" "))}
        if(is.numeric(covar_to_test)){
          #+1 for intercept
          estimates_gene_level[k,]<-fit$coefficients$cond[covar_to_test+1,'Estimate']
          se_gene_level[k,]<-fit$coefficients$cond[covar_to_test+1,'Std. Error']
          stats_all[k,]<-(fit$coefficients$cond[covar_to_test,'z value']+fit$coefficients$zi[covar_to_test,'z value'])/sqrt(2)
          p.vals_gene_level[k,]<-fit$coefficients$cond[covar_to_test,'Pr(>|z|)']
        }
        if(is.character(covar_to_test)){
          if(sum(grepl(covar_to_test,names_cond))>1){num_err=1;stop(paste(c("covar_to_test matches to multiple variables in mean model. Please rename these other variables to not contain the name of covar_to_test. Variable names in the mean model are: ",names_cond),collapse=" "))}
          if(sum(grepl(covar_to_test,names_zi))>1){num_err=1;stop(paste(c("covar_to_test matches to multiple variables in ZI model. Please rename these other variables to not contain the name of covar_to_test. Variable names in the ZI model are: ",names_zi),collapse=" "))}
          if(sum(grepl(covar_to_test,names_cond))==0){num_err=1;stop(paste(c("covar_to_test not found in mean model. Variable names in the mean model are: ",names_cond),collapse=" "))}
          if(sum(grepl(covar_to_test,names_zi))==0){num_err=1;stop(paste(c("covar_to_test not found in ZI model. Variable names in the ZI model are: ",names_zi),collapse=" "))}
          estimates_gene_level[k,]<-fit$coefficients$cond[grepl(covar_to_test,names_cond),'Estimate']
          se_gene_level[k,]<-fit$coefficients$cond[grepl(covar_to_test,names_cond),'Std. Error']
          stats_all[k,]<-(fit$coefficients$cond[grepl(covar_to_test,names_cond),'z value']+fit$coefficients$zi[grepl(covar_to_test,names_zi),'z value'])/sqrt(2)
          p.vals_gene_level[k,]<-fit$coefficients$cond[grepl(covar_to_test,names_cond),'Pr(>|z|)']
        }
      }
      if(statistic=="contrast"){
        fit_twosigmag<-twosigma_custom(counts,silent=TRUE
          ,mean_form=mean_form,zi_form=zi_form
          ,id=id,return_summary_fits = FALSE,weights=weights,internal_call=TRUE,chunk_size=1)
        #if(!fit_twosigmag[[1]]$sdr$pdHess){break}
        fit<-summary(fit_twosigmag$fit[[1]])
        if(return_summary_fits==TRUE){
          fits[[k]]<-fit
        }
        logLik[k]<-as.numeric(fit$logLik)
        gene_err[k]<-(is.na(logLik[k]) | fit_twosigmag$fit[[1]]$sdr$pdHess==FALSE)
        names<-rownames(fit$coefficients$cond)
        names_zi<-rownames(fit$coefficients$zi)
        if(sum(grepl("Intercept",names))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the mean model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the mean model are:",names),collapse=" "))}
        names_zi<-rownames(fit$coefficients$zi)
        if(!is.null(names_zi)&sum(grepl("Intercept",names_zi))>1){num_err=1;stop(paste(c("There seems to be two intercept terms present in the ZI model. Please remove the intercept from either the argument mean_form or from the model.matrix inputted. Variable names in the ZI model are:",names_zi),collapse=" "))}
        if(!is.null(factor_name)){
          frame<-model.frame(fit_twosigmag$fit[[1]])
          index<-grepl(factor_name,colnames(frame))
          if(!is.factor(frame[,index])){num_err=1;stop("Variable factor_name does not appear to be a factor. Please ensure it is of class factor or change other model options if you wish to test a numeric variable (e.g. change to statistic='Z' and set covar_to_test as the numeric variable to test).")}
          if(sum(grepl(factor_name,names))==0){num_err=1;stop("factor_name not found in mean model.")}
          if(sum(grepl('Intercept',names))==0){num_err=1;stop("An intercept is needed in the model for the factor contrast to be successfully calculated using the glht and mcp functions")}
          if((sum(grepl('Intercept',names))==1 & (1+sum(grepl(factor_name,names)))!=ncol(contrast_matrix))){num_err=1;stop(paste(c("Argument contrast_matrix should have number of columns exactly equal to the number of levels of the factor to be tested.  This is true even though the model contains an Intercept, as in this case one level of the factor is automatically dropped by R for model fitting but the contrast of the factor is still tested properly. There should not be columns corresponding to any other covariates or the Intercept in argument 'contrast_matrix'. Leave argument 'factor_name' as NULL if you wish to specify a contrast consisting of more than the factor given by argument 'factor_name'. Variable names in the mean model are: ",names),collapse=" "))}
          mcp<- mcp(factor_name = contrast_matrix)
          names(mcp)<-factor_name
          temp2<-glht_glmmTMB(fit_twosigmag$fit[[1]],
            linfct = mcp)
          temp<-summary(temp2,test=adjusted("none"))
        }else{
          if(length(names)!=ncol(contrast_matrix)){num_err=1;stop(paste0(c("Argument 'contrast_matrix' has ",ncol(contrast_matrix)," columns but should have ",length(names)," columns for the model specified (it must include a column for all covariates including the Intercept). If you are only testing levels of a factor, consider setting the argument 'factor_name'.Variable names in the mean model are:",names),collapse=" "))}
          temp2<-glht_glmmTMB(fit_twosigmag$fit[[1]],
            linfct = contrast_matrix)
          temp<-summary(temp2,test=adjusted("none"))
        }
        stats_all[k,]<-temp$test$tstat
        p.vals_gene_level[k,]<-2*pnorm(-1*abs(temp$test$tstat))
        #temp<-summary(temp2)
        estimates_gene_level[k,]<-temp$test$coefficients
        se_gene_level[k,]<-temp$test$sigma
        residuals_all[k,]<-residuals(fit_twosigmag$fit[[1]])
      }
      gc()

    }
    if(return_summary_fits==TRUE){
      return(list(stats_all=stats_all,p.vals_gene_level=p.vals_gene_level,se_gene_level=se_gene_level,
            estimates_gene_level=estimates_gene_level,residuals_all=residuals_all,logLik=logLik,gene_err=gene_err,fits=fits))
    }else{
      return(list(stats_all=stats_all,p.vals_gene_level=p.vals_gene_level,se_gene_level=se_gene_level,
                  estimates_gene_level=estimates_gene_level,residuals_all=residuals_all,logLik=logLik,gene_err=gene_err))
    }

  }
  size=chunk_size
  chunks<-split(genes,ceiling(seq_along(genes)/size))
  nchunks<-length(chunks)
  cl=NULL

  vars<-unique(c(all.vars(mean_form)[-1],all.vars(mean_form_null)[-1]
                 ,all.vars(zi_form),all.vars(zi_form_null)),count_matrix)
  check_forNA<-function(x){ifelse(any(is.na(get(x))),TRUE,FALSE)}
  #browser()
  if(any(unlist(lapply(vars,FUN=check_forNA)))){stop("NA values are not allowed in covariate vectors or count_matrix for gene set testing. Please remove observations (in both count_matrix and covariates) associated with NA values.")}
  if(cluster_type=="Sock" & ncores>1){
    cl<-parallel::makeCluster(ncores)
    registerDoParallel(cl)
    vars<-unique(c(all.vars(mean_form)[-1],all.vars(mean_form_null)[-1]
                   ,all.vars(zi_form),all.vars(zi_form_null)))
    #vars<-vars[!vars=="id"]
    clusterExport(cl,varlist=vars,envir=environment())
  }
 if(cluster_type=="Fork"&ncores>1){
   cl<-parallel::makeForkCluster(ncores,outfile="")
   registerDoParallel(cl)
 }
  #browser()
  pboptions(type="timer")
  if(lb==TRUE){pboptions(use_lb=TRUE)}

   a<-pblapply(chunks,FUN=fit_tsg,statistic=statistic,
          covar_to_test=covar_to_test,factor_name=factor_name,contrast_matrix=contrast_matrix,
          id=id,ncomps=ncomps,cl=cl)
  if(ncores>1){parallel::stopCluster(cl)}
   #browser()
  rm(count_matrix)
  gc()
  #browser()
  stats_all<-do.call(rbind,sapply(a,'[',1))
  p.vals_gene_level<-do.call(rbind,sapply(a,'[',2))
  se_gene_level<-do.call(rbind,sapply(a,'[',3))
  estimates_gene_level<-do.call(rbind,sapply(a,'[',4))
  residuals_all<-do.call(rbind,sapply(a,'[',5))
  logLik<-unlist(sapply(a,'[',6))
  gene_err<-unlist(sapply(a,'[',7))

  # Fit all gene level statistics that are needed
  if(return_summary_fits==TRUE){
    fit<-do.call(c,sapply(a,'[',8))
  }
rm(a)
gc()
#browser()
  p.val<-matrix(NA,nrow=nsets,ncol=ncomps)
  #p.val_FDR<-matrix(NA,nrow=nsets,ncol=ncomps)
  rho_est<-rep(NA,length=nsets)
  direction<-matrix(NA,nrow=nsets,ncol=ncomps)
  estimates_set_level<-matrix(NA,nrow=nsets,ncol=ncomps)
  for(i in 1:nsets){
    stats_test<-stats_all[index_test[[i]],,drop=FALSE]
    residuals_test<-residuals_all[index_test[[i]],]
    stats_ref<-stats_all[index_ref[[i]],,drop=FALSE]
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
    #tryCatch({if(!allow_neg_corr & rho_est[i]<0){rho_est[i]<-0}},error=function(e){})
    est_corr<-rho_est[i]
    tryCatch({if(!allow_neg_corr & est_corr<0){est_corr<-0}},error=function(e){})
    # Missing values will get dropped in rank statement
    # so make sure sizes of test and reference sets don't include them
    # Don't want to drop missing values because want them documented in output
    for(b in 1:ncomps){
      estimates_set_level[i,b]<-mean(estimates_gene_level[index_test[[i]],b],na.rm=T)
      direction[i,b]<-ifelse(sign(estimates_set_level[i,b])==1,"Up","Down")

      # Missing values will get dropped in rank statement
      # so make sure sizes of test and reference sets don't include them
      # Don't want to drop missing values because want them documented in output
      test_size<-length(stats_test[,b][!is.na(stats_test[,b])])
      ref_size<-length(stats_ref[,b][!is.na(stats_ref[,b])])

      var<-(1/(2*pi))*test_size*ref_size*(asin(1)+(ref_size-1)*asin(.5)+(test_size-1)*(ref_size-1)*asin(.5*est_corr)+(test_size-1)*asin((est_corr+1)/2))
      wilcox_stat<-sum(rank(c(stats_test[,b],stats_ref[,b]),na.last = NA)[1:test_size]) - .5*test_size*(test_size+1)
      p.val[i,b]<-2*pnorm(-1*abs((wilcox_stat-.5*test_size*ref_size)/sqrt(var)))

      # T-test output below

      #n_genes_tested<-test_size+ref_size
      #delta<-(n_genes_tested)/ref_size*(mean(stats_test[,b],na.rm=T)-mean(c(stats_test[,b],stats_ref[,b]),na.rm=T))
      #vif<-1+(test_size-1)*rho_est[i]
      #varStatPooled<-((n_genes_tested-1)*var(c(stats_test[,b],stats_ref[,b]),na.rm=T)-delta^2*test_size*ref_size/n_genes_tested)/(n_genes_tested-2)
      #two.sample.t <- delta / sqrt( varStatPooled * (vif/test_size + 1/ref_size) )
      #p.val_ttest[i,b]<-2*pt(-1*abs(two.sample.t),df=n_genes_tested-2)
      # This only happens if all genes in the test set are NA
      # In this case p-value will report as zero when it's really NA
      if(test_size==0 | ref_size==0){p.val[i,b]<-NA
      #p.val_ttest[i,b]<-NA
      }
    }
    if(i%%100==0){print(paste0("Set ",i," of ",nsets," Finished"))}
  }
  #browser()
  p.val_FDR<-matrix(apply(p.val,MARGIN=2,FUN=p.adjust,method="fdr"),nrow=nsets,ncol=ncomps)
  rownames(stats_all)<-genes
  colnames(p.val)<-rownames(contrast_matrix)
  rownames(p.val)<-names(index_test)
  colnames(p.val_FDR)<-rownames(contrast_matrix)
  rownames(p.val_FDR)<-names(index_test)
  #rownames(p.val_ttest)<-names(index_test)
  rownames(direction)<-names(index_test)
  names(rho_est)<-names(index_test)
  rownames(p.vals_gene_level)<-gene_names
  rownames(estimates_gene_level)<-gene_names
  rownames(se_gene_level)<-gene_names
  rownames(estimates_set_level)<-names(index_test)
  names(logLik)<-gene_names
  names(gene_err)<-gene_names
  #browser()
  if(statistic=="contrast"){
    colnames(se_gene_level)<-rownames(contrast_matrix)
    colnames(estimates_gene_level)<-rownames(contrast_matrix)
    colnames(estimates_set_level)<-rownames(contrast_matrix)
    colnames(direction)<-rownames(contrast_matrix)
    colnames(p.vals_gene_level)<-rownames(contrast_matrix)
    #colnames(p.val_ttest)<-rownames(contrast_matrix)
    colnames(stats_all)<-rownames(contrast_matrix)
  }
  # if(return_summary_fits==TRUE){
  #   names(fit)<-gene_names
  #   return(list(gene_summary_fits=fit,stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,set_p.val_ttest=p.val_ttest,estimates_gene_level=estimates_gene_level,se_gene_level=se_gene_level,estimates_set_level=estimates_set_level,direction=direction,corr=rho_est,gene_level_logLik=logLik,gene_error=gene_err,test_sets=index_test,ref_sets=index_ref))
  # }else{
  #   return(list(stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,set_p.val_ttest=p.val_ttest,estimates_gene_level=estimates_gene_level,se_gene_level=se_gene_level,estimates_set_level=estimates_set_level,direction=direction,corr=rho_est,gene_level_logLik=logLik,gene_error=gene_err,test_sets=index_test,ref_sets=index_ref))
  # }
  if(return_summary_fits==TRUE){
    names(fit)<-gene_names
    return(list(gene_summary_fits=fit,stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,set_p.val_FDR=p.val_FDR,estimates_gene_level=estimates_gene_level,se_gene_level=se_gene_level,estimates_set_level=estimates_set_level,direction=direction,corr=rho_est,gene_level_logLik=logLik,gene_error=gene_err,test_sets=index_test,ref_sets=index_ref))
  }else{
    return(list(stats_gene_level_all=stats_all,p.vals_gene_level=p.vals_gene_level,set_p.val=p.val,set_p.val_FDR=p.val_FDR,estimates_gene_level=estimates_gene_level,se_gene_level=se_gene_level,estimates_set_level=estimates_set_level,direction=direction,corr=rho_est,gene_level_logLik=logLik,gene_error=gene_err,test_sets=index_test,ref_sets=index_ref))
  }
}
