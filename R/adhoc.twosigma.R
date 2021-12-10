##' adhoc.twosigma: Perform the ad hoc method described in TWO-SIGMA paper
##' @export adhoc.twosigma
##' @param count Vector of non-negative integer read counts.
##' @param mean_covar Covariates for the (conditional) mean model. Must be a matrix (without an intercept column) or = 1 to indicate an intercept only model.
##' @param zi_covar Covariates for the zero-inflation model. Must be a matrix (without an intercept column), = 1 to indicate an intercept only model, or = 0 to indicate no zero-inflation model desired.
##' @param id Vector of individual-level ID's. Used as predictor in ANOVA model.
##' @param weights weights, as in glm. Defaults to 1 for all observations and no scaling or centering of weights is performed. Passed into zeroinfl function.
##' @return P-value from the ANOVA F test.
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
##' # Run adhoc.twosigma
##'
##' adhoc.twosigma(sim_dat[1,],mean_covar = X,zi_covar=Z,id = id)

adhoc.twosigma<-function(count,mean_covar,zi_covar,id
                        ,weights=rep(1,length(count))){
#browser()
  check_twosigma_input(count,mean_covar,zi_covar
    ,mean_re=TRUE,zi_re=TRUE,id=id
    ,disp_covar=NULL)

  form<-create_adhoc_formulas(count,mean_covar,zi_covar)
resid<-residuals(zeroinfl(form,dist="negbin",link="logit"),type="pearson",weights=weights)
     lm1<-lm(resid~id)
     p.val_adhoc<-anova(lm1)$`Pr(>F)`[1]
     return(p.val_adhoc)
}
