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
##' @param sim.seed Random seed to be used.
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

  #simulate random effects
  a <- as.matrix(rnorm(nind,mean=0,sd=sigma.a))
  a.rep <- rep(a,times=ncellsper)
  b <- as.matrix(rnorm(nind,mean=0,sd=sigma.b))
  b.rep <- rep(b,times=ncellsper)
  logit.p<-Z%*%as.matrix(alpha,ncol=1) +a.rep
  log.mu<-X%*%as.matrix(beta,ncol=1)+b.rep

  p  <- exp(logit.p)/(1+exp(logit.p)) # inverse logit function
  mu <- exp(log.mu)

  Y<-rep(NA,sum(ncellsper))
  ind.dropout <- rbinom(length(Y), 1, p)
  for (i in 1:length(Y)){
    if(ind.dropout[i] == 1){Y[i]=0}
    if(ind.dropout[i] == 0)
    {
      if(phi>0)
      {
        Y[i]<-rnbinom(1,size=(1/phiinv),prob=(1/(1+phiinv*mu[i])))
      }
      #Watch out for negative binomial parameterizations
    }
  }
  return(list(Y=Y,X=X[,-1],Z=Z[,-1],a=a,b=b,alpha=alpha,beta=beta,sigma.a=sigma.a,sigma.b=sigma.b,
    nind=nind,ncellsper=ncellsper,id=id))
}

