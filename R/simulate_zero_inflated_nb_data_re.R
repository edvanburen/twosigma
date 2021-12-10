##' Simulated zero-inflated negative binomial data with random effects
##' @param ncellsper Vector giving the number of cells per individual.  Length of the vector is taken as the number of individuals.
##' @param X Covariate matrix (without intercept) for the (conditional) mean model.
##' @param Z Covariate matrix (without intercept) for the zero-inflation model.
##' @param alpha Column vector of true parameters from the zero-inflation model. Number of rows must match number of columns in Z.
##' @param beta Column vector of true parameters from the (conditional) mean model. Number of rows must match number of columns in X.
##' @param phi Overdispersion parameter for the negative binomial distribution (see details for more about parameterization).
##' @param sigma.a Standard deviation for the zero-inflation model random intercept.
##' @param sigma.b Standard deviation for the (conditional) mean random intercept.
##' @param id.levels Individual-level IDs. If NULL set as 1,2,... up to the number of individuals.
##' @param sim.seed Random seed to be used. If NULL one will be randomly chosen.
##' @return Y Simulated counts
##' @return X Covariate matrix (without intercept) for the (conditional) mean model.
##' @return Z Covariate matrix (without intercept) for the zero-inflation model.
##' @return a Random effects for the zero-inflation model.
##' @return b Random effects for the (conditional) mean model.
##' @return alpha Column vector of true parameters from the zero-inflation model. Number of rows must match number of columns in Z.
##' @return beta Column vector of true parameters from the (conditional) mean model. Number of rows must match number of columns in X.
##' @return phi Overdispersion parameter for the negative binomial distribution (see details for more about parameterization).
##' @return sigma.a Standard deviation for the zero-inflation model random intercept.
##' @return sigma.b Standard deviation for the (conditional) mean random intercept.
##' @return nind Number of individuals.
##' @return ncellsper Vector giving the number of cells per individual.
##' @return id.levels Individual-level IDs.
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
##' @export simulate_zero_inflated_nb_random_effect_data



simulate_zero_inflated_nb_random_effect_data<-function(ncellsper,X,Z,alpha,beta,phi,sigma.a,sigma.b,
                                                      id.levels=NULL,sim.seed=NULL)
{
  if(is.null(sim.seed)){sim.seed<-sample.int(1e8,1)}
  if(phi<=0){
    stop("phi must be >0")
  }
  if(sigma.a<0 | sigma.b<0){
    stop("sigma.a and sigma.b cannot be less than zero")
  }
  if(!is.null(sim.seed)){
    set.seed(sim.seed)
  }
  # Define for use with rnbinom statement below
  phiinv<-1/phi
  id.levels<-1:length(ncellsper)
  nind<-length(id.levels)
  id<-rep(id.levels,times=ncellsper)

  names<-colnames(Z)
  Z<-cbind(1,Z)
  colnames(Z)<-c("Intercept",names)

  names<-colnames(X)
  X<-cbind(1,X)
  colnames(X)<-c("Intercept",names)

  #simulate random effects at sample level
  # ".rep" vectors ensure the dimensions match
  # the number of cells for each individual
  a <- as.matrix(rnorm(nind,mean=0,sd=sigma.a))
  a.rep <- rep(a,times=ncellsper)
  b <- as.matrix(rnorm(nind,mean=0,sd=sigma.b))
  b.rep <- rep(b,times=ncellsper)

  # drop-out probability and mean
  logit.p<-Z%*%as.matrix(alpha,ncol=1) +a.rep
  log.mu<-X%*%as.matrix(beta,ncol=1)+b.rep
  p  <- exp(logit.p)/(1+exp(logit.p)) # inverse logit function
  mu <- exp(log.mu)

  #Y gives the simulated counts
  Y<-rep(NA,sum(ncellsper))
  ind.dropout <- rbinom(length(Y), 1, p)
  for (i in 1:length(Y)){
    if(ind.dropout[i] == 1){Y[i]=0}
    if(ind.dropout[i] == 0)
    {
        Y[i]<-rnbinom(1,size=(1/phiinv),prob=(1/(1+phiinv*mu[i])))
      #Watch out for negative binomial parameterizations
    }
  }
  # return X and Z without the intercept for convenience
  return(list(Y=Y,X=X[,-1],Z=Z[,-1],a=a,b=b,alpha=alpha,beta=beta,sigma.a=sigma.a,sigma.b=sigma.b,
    nind=nind,ncellsper=ncellsper,id=id))
}

