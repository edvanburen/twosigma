##' TWO-SIGMA: Fit the TWO-component SInGle cell Model-based Association method of...
##' @param count Vector of non-negative integer read counts. If specifying custom formula(s) via the arguments mean_form and zi_form the expression in mean_form will supersede.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or = 1 to indicate an intercept only model.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column), = 1 to indicate an intercept only model, or = 0 to indicate no zero-inflation model desired.
##' @param mean_re Should random intercepts be included in the (conditional) mean model?
##' @param zi_re Should random intercepts be included in the zero-inflation model?
##' @param id Vector of individual-level ID's. Used for random effect prediction.
##' @param disp_covar Covariates for a log-linear model for the dispersion. Either a matrix of covariates or = 1 to indicate an intercept only model. Random effect terms are not permitted in the dispersion model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed.
##' @param control Control parameters for optimization in \code{glmmTMB}.
##' @section Details: Note that the overdispersion estimate returned is the inverse of phi in the parameterization used in TWO-SIGMA.  See ?nbinom2
##' @return An object of class \code{glmmTMB}. See Details for a note about the estimated overdispersion parameter.
##' @import glmmTMB
##' @import methods
##' @import pscl
##' @importFrom stats anova as.formula lm pchisq rbinom residuals rnbinom rnorm
##' @export twosigma

#zi_covar, mean_covar matrices must be specified? Need to figure out what to do exactly here
# users can input own model formulas to overwrite random effects specification
# should know what theyre doing for that though
twosigma_noadhoc<-function(count,mean_covar,zi_covar
                  ,mean_re=TRUE,zi_re=TRUE
                  ,id
                   ,disp_covar=NULL #need to be able to use data option?
                   ,weights=rep(1,length(count))
                   ,control = glmmTMBControl()){
  #if(!grepl("nbinom2",family$family)){
  #  stop("Only the Negative Binomial Distribution is implemented in TWO-SIGMA")
  #}
  # Check that response is only valid counts if using Negative Binomial
    #if(grepl("nbinom2",family$family) & (sum(!as.matrix(count,ncol=1)%%1==0)>0 | min(count)<0)){

  check_twosigma_input(count,mean_covar,zi_covar
                       ,mean_re,zi_re
                       ,disp_covar)
  formulas<-create_model_formulas(mean_covar,zi_covar
    ,mean_form=NULL,zi_form=NULL
    ,mean_re,zi_re
    ,disp_covar)

  fit<-glmmTMB(formula=formulas$mean_form
            ,ziformula=formulas$zi_form
            ,weights=weights
            ,dispformula = formulas$disp_form
            ,family=nbinom2,verbose = F
            ,control = control)

return(fit)
}


