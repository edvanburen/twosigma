library(testthat)
library(twosigma)


# Simulate data to test exported function simulate_zero_inflated_nb_data_re.R

nind<-10;ncellsper<-rep(50,nind)
sigma.a<-.5;sigma.b<-.5;phi<-.1
alpha<-c(1,0,-.5,-2);beta<-c(2,0,-.1,.6);beta2<-c(2,1,-.1,.6)
id.levels<-1:nind;nind<-length(id.levels)
id<-rep(id.levels,times=ncellsper)
sim.seed<-1234

# individual level covariates
t2d_sim<-rep(rbinom(nind,1,p=.4),times=ncellsper)
cdr_sim<-rbeta(sum(ncellsper),3,6)
age_sim<-rep(sample(c(20:60),size=nind,replace = TRUE),times=ncellsper)

#Construct design matrices
Z<-cbind(scale(t2d_sim),scale(age_sim),scale(cdr_sim))
colnames(Z)<-c("t2d_sim","age_sim","cdr_sim")
X<-cbind(scale(t2d_sim),scale(age_sim),scale(cdr_sim))
colnames(X)<-c("t2d_sim","age_sim","cdr_sim")

# Test exported function simulate_zero_inflated_nb_random_effect_data
# Use output to test other functions

sim_dat<-matrix(nrow=4,ncol=sum(ncellsper))
for(i in 1:nrow(sim_dat)){
  if(i<2){# Null Gene Sets
    sim_dat[i,]<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta2,phi,sigma.a,sigma.b,id.levels=NULL)$Y
  }else{# Alternative Gene Sets
    sim_dat[i,]<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta,phi,sigma.a,sigma.b,id.levels=NULL)$Y
  }
}
rownames(sim_dat)<-paste("Gene",1:4)

# Test exported function adhoc.twosigma
adhoc.twosigma(sim_dat[1,],mean_covar = X,zi_covar=Z,id = id)

# Test exported function lr.twosigma
lr.twosigma(count=sim_dat[1,,drop=FALSE],mean_covar = X,zi_covar = Z,id=id,covar_to_test = 1)


# Test exported function lr.twosigma_custom
lr.twosigma_custom(count=sim_dat[1,,drop=FALSE],mean_form_alt = count~X,mean_form_null = count~X[,-1],zi_form_alt = ~0,zi_form_null = ~0,id=id,lr.df=1)

# Test exported function test.vc.twosigma
test.vc.twosigma(sim_dat[1,,drop=FALSE],mean_covar = X,zi_covar=Z,mean_re = TRUE,zi_re=FALSE,id = id)

# Test exported function twosigma
twosigma(sim_dat[1:2,],mean_covar = X,zi_covar=1,id = id)

# Test exported function twosigma_custom
twosigma_custom(sim_dat[1:2,],mean_form = count~X,zi_form = ~0,id=id)

# Test exported function twosigmag
twosigmag(sim_dat,index_test = list(c(1,3)),all_as_ref = TRUE,mean_form = count~X,zi_form = ~0,id=id,covar_to_test  = "t2d_sim",statistic = "Z")

