# TWOSIGMA (TWO-component SInGle cell Model-based Association method)

## Introduction

<code>twosigma</code> is an R package for differential expression (DE) analysis and gene set testing (GST) in single-cell RNA-seq (scRNA-seq) data.  At the gene-level, DE can be assessed by fitting our proposed TWO-component SInGle cell Model-based Association method (TWO-SIGMA). The first component models the drop-out probability with a mixed effects logistic regression, and the second component models the (conditional) mean read count with a mixed-effects log-linear negative binomial regression. Our approach thus allows for both excess zero counts and overdispersed counts while also accommodating dependency in both drop-out probability and mean mRNA abundance.  TWO-SIGMA is especially useful in its flexibility to analyze DE beyond a two-group comparison while simultaneously controlling for additional subject-level or cell-level covariates including batch effects.  At the set-level, the package can perform competitive gene set testing using our proposed TWO-SIGMA-G method. Users can specify the number of cores to be used for parallelization in all functions using the ncores argument.
## Installation
We recommend installing from Github for the latest version of the code:
```r
install.packages("devtools")
devtools::install_github("edvanburen/twosigma")
library(twosigma)
```

## Gene-Level Model
TWO-SIGMA is based on the following parameterization of the negative binomial distribution: 
<img src="inst/image/nb3.jpg" width="600" align="center">

A point mass at zero is added to the distribution to account for dropout.  The result is the probability mass function for the zero-inflated negative binomial distribution:

<img src="inst/image/zinb.jpg" width="600" align="center">

The full TWO-SIGMA specification is therefore as follows:

<img src="inst/image/twosigma.jpg" width="600" align="center">

## Usage  
The workhorse function is twosigma, which can be easiest called as
```r
twosigma(count_matrix, mean_covar, zi_covar,id,ncores=1)
```

- **count_matrix**: A vector of non-negative integer counts. No normalization is done.
- **mean_covar**: A matrix (such as from model.matrix) of covariates for the (conditional) mean model without an intercept term. Columns give covariates and the number of rows should correspond to the number of cells.
- **zi_covar**: A matrix (such as from model.matrix) of covariates for the zero-inflation model without an intercept term. Columns give covariates and the number of rows should correspond to the number of cells.
- **id**: Vector of individual-level ID's (length equal to the total number of cells). Used for random effect prediction and the ad hoc method and is currently required even if neither is being used.
- **ncores**: Number of cores to use for parallelization. **Multiple cores are recommended if possible.**

By default, we employ our ad hoc procedure to determine if random effects are needed. If users wish to specify their own random effect specifications, they can set adhoc=FALSE, and use the following inputs:

- **mean_re**: Should random intercept terms be included in the (conditional) mean model?
- **zi_re**: Should random intercept terms be included in the zero-inflation model?

If <code> adhoc=TRUE</code>, mean_re and zi_re are ignored and a warning is printed.

If users wish to customize the random effect or fixed effect specification, they may do so via the function <code> twosigma_custom </code>, which has the following basic syntax:
```r
twosigma_custom(count_matrix, mean_form, zi_form, id)
````
- **count_matrix**: A matrix of non-negative integer counts containing no missing values. No normalization is done. Batch can be controlled for by inclusion in the design matrices.
- **mean_form** a two-sided formula for the (conditional) mean model. Left side specifies the response and right side includes fixed and random effect terms. Users should ensure that the response has the name "count", e.g. <code> mean_form = count ~ 1 </code>
- **zi_form** a one-sided formula for the zero-inflation model including fixed and random effect terms, e.g. <code>  ~ 1 </code>
- **id**: Vector of individual-level ID's. Used for random effect prediction.

Some care must be taken, however, because these formulas are used directly. **It is therefore the user's responsibility to ensure that formulas being inputted will operate as expected**. Syntax is identical to the <code> lme4 </code> package.

For example, each of the following function calls reproduces the default TWO-SIGMA specification with random intercepts in both components:

```r
fits<-twosigma(count_matrix, mean_covar=mean_covar_matrix, zi_covar=zi_covar_matrix, mean_re = TRUE, zi_re = TRUE, id=id,adhoc=F)
twosigma_custom(count, mean_form=count~mean_covar_matrix+(1|id),zi_form=~zi_covar_matrix+(1|id),id=id)
```
## Fixed Effect Testing  

# Likelihood Ratio Test
If users wish to jointly test a fixed effect using the twosigma model via the likelihood ratio test, they may do so using the <code> lr.twosigma </code> or <code> lr.twosigma_custom </code> functions:
```r
lr.twosigma(count_matrix, mean_covar, zi_covar, covar_to_test, mean_re = TRUE,zi_re = TRUE, disp_covar = NULL,adhoc=TRUE)
lr.twosigma_custom(count_matrix, mean_form_alt, zi_form_alt, mean_form_null, zi_form_null, id, lr.df)
```
- **covar_to_test**: Either a string indicating the column name of the covariate to test or an integer referring to its column position in BOTH the mean_covar and zi_covar matrices. If an integer is specified there is no check that it corresponds to the same covariate in both the mean_covar and zi_covar matrices. 
- **lr.df** If custom formulas are input users must provide the degrees of freedom from which the likelihood ratio p-value can be calculated. Must be a non-negative integer. 

The <code> lr.twosigma </code> function assumes that the variable being tested is in both components of the model (and thus that the zero-inflation component exists and contains more than an Intercept). Users wishing to do fixed effect testing in other cases can use the <code> lr.twosigma_custom </code> function with custom formulas or construct the test themselves using two calls to <code>twosigma</code> or <code> twosigma_custom</code>. The formula inputs <code> mean_form_alt </code>, <code> mean_form_null</code>, <code> zi_form_alt</code>, and `zi_form_null` should be specified as in the <code> lr.twosigma_custom</code> function and once again **users must ensure custom formulas represent a valid likelihood ratio test**.  One part of this responsibility is specifying the argument \code{lr.df} giving the degrees of freedom of the likelihood ratio test.

# Z-statistic or Stouffer statistic
Assume \code{fits} is an object returned from \code{twosigma} or \code{twosigma_custom}.  Then, we can get some gene-level staistics using:
```r
calc_logFC<-function(x){if(class(x)=="glmmTMB"){x<-summary(x)};x$coefficients$cond['t2d_sim','Estimate']}
calc_Z<-function(x){if(class(x)=="glmmTMB"){x<-summary(x)};x$coefficients$cond['t2d_sim','Estimate']}
calc_p.vals<-function(x){if(class(x)=="glmmTMB"){x<-summary(x)};x$coefficients$cond['t2d_sim',4]}

logFC<-unlist(lapply(fits,calc_logFC))
Zstats<-unlist(lapply(fits,calc_Z))
p.val_Zstat<-unlist(lapply(fits,calc_p.vals)

calc_Stouffer<-function(x){(x$coefficients$cond['t2d_sim','Estimate'] + x$coefficients$zi['t2d_sim','Estimate'])/sqrt(2)}
```

# Testing of a contrast matrix
Functionality to do this is ongoing.  For now, individuals can use the \code{twosigmag} to test custom contrast matrices (even if not interested in gene set testing, simply set \code{index_test=list(c(1,2))} and \code{all_as_ref=TRUE} and look at gene-level output).

## Ad hoc method
As mentioned in the paper, we mention a method that can be useful in selecting genes that may benefit from the inclusion of random effect terms. This method fits a zero-inflated negative binomial model without random effects and uses a one-way ANOVA regressing the Pearson residuals on the individual ID to look for differences between individuals.

```r
adhoc.twosigma(count, mean_covar, zi_covar, id)
```
The p-value from the ANOVA F test is returned, and can be used as a screening for genes that are most in need of random effects. This functionality is built into the <code> twosigma </code> function so users likely do not need to call directly themselves.

## Gene-set Testing
Competitive gene set testing can be performed using the function `twosigmag`. Gene-level statistics currently implemented include likelihood ratio, Z-statistic from the mean model, Stouffer's combination of the Z-statistics from the mean and ZI model, or a test of a custom contrast matrix.  If a contrast matrix is input, set-level results are returned for each row of the contrast. **Multiple cores are once again recommended if possible, particularly if using the likelihood ratio test.**
 
 **The adhoc procedure is not recommended for use in gene set testing**.  This is because geneset testing relies on a common gene-level null hypothesis being tested.  When some genes have random effects and others do not, it is not clear that this requirement is met. Arguments which require more explanation over above are given as follows:

- **index_test**: A list of numeric vectors corresponding to indices (row numbers from ) belonging to the test set(s) of interest.  
- **index_ref**: A list of numeric vectors, corresponding to indices (row numbers) belonging to desired reference set(s). Most users should not need to modify this option.
- **all_as_ref**: Should all genes not in the test set be taken as the reference set? Defaults to true.  If \code{FALSE}, a random sample of size identical to the size of the test set is taken as the reference. 
- **allow_neg_corr**: Should negative correlations be allowed for a test set? By default negative correlations are set to zero to be conservative. Most users should not need to modify this option.
- **statistic**: Gene-level statistic that should be used for determining set-level enrichment.  Options include "LR" for likelihood ratio, "Z" for the Z-statistic from the mean model, "Stouffer" for Stouffer's combination of the Z-statistics from the mean and ZI model, or 'contrast' for a test of a custom contrast matrix.
- **covar_to_test**: Which covariate should be used to determine gene-level significance if statistic = "LR", "Z", or "Stouffer".
- **contrast_matrix**: Contrast matrix to be used if statistic = "contrast". Each row of the matrix will have separate gene-level and set-level statistics.  Rownames of \code{contrast_matrix} should correspond to a meaningful name of the hypothesis for nicely formatted output. If testing a factor, must have a number of columns exactly equal to the number of levels of the factor.  Otherwise, must have one column per parameter in the mean model (including a column for the intercept.)
- **factor_name**: Name of the factor being tested by \code{contrast_matrix}. Needed if and only if \code{statistic='contrast'} and \code{contrast_matrix} is testing a factor variable in the mean model.

**Multiple cores are once again recommended if possible, particularly if using the likelihood ratio test.**

## Examples
```r

#--------------------------------------------------
#--- Simulate Data
#--------------------------------------------------


#Set parameters for the simulation
set.seed(1)
nind<-10
ncellsper<-rep(1000,nind)
sigma.a<-.1
sigma.b<-.1
alpha<-c(1,0,-.5,-2)
beta<-c(2,0,-.1,.6)
phi<-.1
id.levels<-1:nind


# Simulate some individual level covariates as well as cell-level Cellular Detection Rate (CDR)
t2d_ind<-rbinom(nind,1,p=.4)
t2d_sim<-rep(t2d_ind,times=ncellsper)
nind<-length(id.levels)
id<-rep(id.levels,times=ncellsper)
cdr_sim<-rbeta(sum(ncellsper),3,6)
age_sim_ind<-sample(c(20:60),size=nind,replace = TRUE)
age_sim<-rep(age_sim_ind,times=ncellsper)

#Construct design matrices
Z<-cbind(t2d_sim,age_sim,cdr_sim)
colnames(Z)<-c("t2d_sim","age_sim","cdr_sim")
X<-cbind(t2d_sim,age_sim,cdr_sim)
colnames(X)<-c("t2d_sim","age_sim","cdr_sim")

sim_dat<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta,phi,sigma.a,sigma.b,
                                                      id.levels=NULL)
sim_dat2<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta,phi,sigma.a,sigma.b,
                                                      id.levels=NULL)
#--------------------------------------------------
#--- Fit TWO-SIGMA to simulated data
#--------------------------------------------------
id<-sim_dat$id
counts<-matrix(rbind(sim_dat$Y,sim_dat2$Y),nrow=2,byrow=FALSE)
rownames(counts)<-paste0("Gene ",1:2)

fit<-twosigma(counts,zi_covar=Z,mean_covar = X,id=id,mean_re=FALSE,zi_re=FALSE,adhoc=F,ncores=1)
fit2<-twosigma_custom(counts, mean_form=count~t2d_sim+age_sim+cdr_sim
                      ,zi_form=~t2d_sim+age_sim+cdr_sim,id=id,ncores=1)

#fit and fit2 are the same for the both genes

fit[[1]];fit2[[1]]
fit$`Gene 2`;fit2$`Gene 2`

#--- Fit TWO-SIGMA without a zero-inflation component
fit_noZI<-twosigma(counts,zi_covar=0,mean_covar = X,id=id,mean_re=F,zi_re=F,adhoc=F),ncores=1
fit_noZ2I<-twosigma_custom(counts,zi_form=~0,mean_form=count~X,id=id,ncores=1)

#--- Fit TWO-SIGMA with an intercept only zero-inflation component and no random effects
fit_meanZI<-twosigma(counts,zi_covar=1,mean_covar = X,id=id,mean_re=F,zi_re=F,adhoc=F,ncores=1)
fit_meanZI2<-twosigma_custom(counts, mean_form=count~t2d_sim+age_sim+cdr_sim,zi_form=~1,id=id,ncores=1)

fit_noZI[[1]]
fit_meanZI[[1]]

# Perform Likelihood Ratio Test on variable "t2d_sim"
          
lr.fit<-lr.twosigma(counts,covar_to_test="t2d_sim",mean_covar = X,zi_covar=Z,id=id)
lr.fit$LR_stat
lr.fit$LR_p.val

lr.fit_custom<-lr.twosigma_custom(counts,mean_form_alt=count~t2d_sim+age_sim+cdr_sim
, zi_form_alt=~t2d_sim+age_sim+cdr_sim, mean_form_null=count~age_sim+cdr_sim
,zi_form_null=~age_sim+cdr_sim,id=id,lr.df=2)

lr.fit_custom$LR_stat
lr.fit_custom$LR_p.val

# Collect gene-level Z-statistics for variable 't2d_sim':


calc_logFC<-function(x){if(class(x)=="glmmTMB"){x<-summary(x)};x$coefficients$cond[2,'Estimate']}
calc_Z<-function(x){if(class(x)=="glmmTMB"){x<-summary(x)};x$coefficients$cond[2,'z value']}
calc_p.vals<-function(x){if(class(x)=="glmmTMB"){x<-summary(x)};x$coefficients$cond[2,4]}

logFC<-unlist(lapply(fit,calc_logFC))
Zstats<-unlist(lapply(fit,calc_Z))
p.val_Zstat<-unlist(lapply(fit,calc_p.vals))

#--------------------------------------------------
# Perform Gene-Set Testing
#--------------------------------------------------

# First, simulate some DE genes and some non-DE genes

sim_dat2<-matrix(nrow=10,ncol=sum(ncellsper))
beta<-c(2,0,-.1,.6)
beta2<-c(2,.5,-.1,.6)
for(i in 1:nrow(sim_dat2)){
  if(i<5){
    sim_dat2[i,]<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta2,phi,sigma.a,sigma.b=.5,
      id.levels=NULL)$Y
  }else{
    sim_dat2[i,]<-simulate_zero_inflated_nb_random_effect_data(ncellsper,X,Z,alpha,beta,phi,sigma.a,sigma.b,
      id.levels=NULL)$Y
  }
}

# Use LR statistic
gst<-twosigmag(sim_dat2,index_test = list(c(6:10)),index_ref = list(c(1:5)),
mean_form=count~t2d_sim+age_sim+cdr_sim,zi_form=~t2d_sim+age_sim+cdr_sim
,mean_form_null=count~age_sim+cdr_sim,
zi_form_null=~age_sim+cdr_sim,statistic="LR",
covar_to_test='t2d_sim',lr.df=2,id=id)

gst$set_p.val

# Use Z-statistic

gst2<-twosigmag(sim_dat2,index_test = list("Set 1" = c(6:10)),mean_form = count~t2d_sim+age_sim+cdr_sim
,zi_form = ~t2d_sim+age_sim+cdr_sim,id=id,covar_to_test  = "t2d_sim",ncores = 1,statistic = "Z",ncores=1)

gst2$set_p.val

# Testing a simple contrast equivalent to the Z statistic for 't2d_sim'

gst3<-twosigmag(sim_dat2,index_test = list("Set 1" = c(6:10)),mean_form = count~t2d_sim+age_sim+cdr_sim
,zi_form = ~t2d_sim+age_sim+cdr_sim,id=id,statistic = "contrast"
,contrast_matrix = matrix(c(0,1,0,0),nrow=1),ncores = 1)

# Same result as using Z test

gst3$set_p.val


# Testing a contrast of a factor variable
# set seed to make sure factor has all three levels 
# so contrast matrix is properly defined
set.seed(1234) 
fact<-factor(rep(sample(c(0,1,2),nind,replace=T),times=ncellsper))
cont_matrix<-matrix(c(-1,1,0,-1,0,1),nrow=2,byrow = T)
rownames(cont_matrix)<-c("Test 1","Test 2")

gst4<-twosigmag(sim_dat2,index_test = list("Set 1" = c(6:10))
,mean_form = count~t2d_sim+age_sim+cdr_sim+fact
,zi_form = ~t2d_sim+age_sim+cdr_sim,id=rep(id.levels,times=ncellsper)
,statistic = "contrast",contrast_matrix = cont_matrix
  ,factor_name="fact",ncores = 1,return_summary_fits = T,ncores=1)

# Finally, test the factor "manually to show results are the same
cont_matrix2<-matrix(c(0,0,0,0,1,0,0,0,0,0,0,1),nrow=2,byrow = T)
rownames(cont_matrix2)<-c("Test 1","Test 2")

gst5<-twosigmag(sim_dat2,index_test = list("Set 1" = c(6:10)),mean_form = count~t2d_sim+age_sim+cdr_sim+fact
,zi_form = ~t2d_sim+age_sim+cdr_sim,id=id,statistic = "contrast",contrast_matrix = cont_matrix2
,ncores = 1,return_summary_fits = T,ncores=1)
  
#Two give the same results

gst4$set_p.val
gst5$set_p.val

```

