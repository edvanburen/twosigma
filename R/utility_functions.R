# Should we have NULL for covars if requiring their input?
# Should we extend ad hoc to non ZI case?
check_twosigma_input<-function(count,mean_covar=NULL,zi_covar=NULL
                               ,mean_form=NULL,zi_form=NULL
                               ,mean_re=TRUE,zi_re=TRUE
                               ,disp_covar=NULL){

  # Override count with mean_form if specified then check the inputs
  if(!is.null(mean_form)){
      count<-mean_form[[1]]
  }

  if(sum(!as.matrix(count,ncol=1)%%1==0)>0 | min(count)<0){
    stop("When using the Negative Binomial Distribution data must contain only non-negative integers")
  }

  if(is.null(mean_covar)&is.null(mean_form))
  {
    stop("At minimum must specify Mean Model Covariates or provide model formula")
  }

  #Check for intercept being input
  if(class(mean_covar)=="matrix"){
    check_no_intercept<-function(x){mean(x==1)}
    if(any(apply(mean_covar,2,check_no_intercept)==1)){
      stop("Please remove intercept from mean model covariate matrix")
    }
  }
  if(class(zi_covar)=="matrix"){
      check_no_intercept<-function(x){mean(x==1)}
      if(any(apply(zi_covar,2,check_no_intercept)==1)){
      stop("Please remove intercept from ZI model covariate matrix")
      }
  }
      if(class(disp_covar)=="matrix"){
        check_no_intercept<-function(x){mean(x==1)}
        if(any(apply(disp_covar,2,check_no_intercept)==1)){
          stop("Please remove intercept from dispersion covariate matrix")
        }
      }

  if(mean(count==0)>.9){
    warning("More than 90% of data are zeros. Mean model results may be misleading for such sparse data")
  }
  if(mean(count==0)<.1){
    if(!(is.atomic(zi_covar)&length(zi_covar)==1)){
      warning("Less than 10% of data are zeros. Zero-Inflation model results may be misleading or unnecessary")
    }else{
      if(is.atomic(zi_covar) & length(zi_covar)==1& zi_covar==1){
        warning("Less than 10% of data are zeros. Zero-Inflation model results may be misleading or unnecessary")
      }
    }
  }

  #if(!is.null(zi_form)){

# }
  # if(mean(count==0)<.1 & !is.atomic(zi_covar)){
  #     warning("Less than 10% of data are zeros. Zero-Inflation model results may be misleading or unnecessary")
  #   }else if(mean(count==0)<.1 & is.atomic(zi_covar) & length(zi_covar)==1& zi_covar==1){
  #     warning("Less than 10% of data are zeros. Zero-Inflation model results may be misleading or unnecessary")
  # }

}
create_model_formulas<-function(mean_covar=NULL,zi_covar=NULL
                                ,mean_form=NULL,zi_form=NULL
                                ,mean_re=TRUE,zi_re=TRUE
                                ,disp_covar=NULL,disp_form=NULL){
if(is.null(zi_form)){
    if(is.null(zi_covar)){
      if(zi_re==TRUE){
        zi_form<- ~1+(1|id)} #Default is intercept only
      if(zi_re==FALSE){
        zi_form<- ~1} #Default is intercept only
    }else{
      if(is.atomic(zi_covar)&length(zi_covar)==1){
        if(zi_covar==0){zi_form<- ~0}else{
          if(zi_covar==1){
            if(zi_re==TRUE){
              zi_form<-~1+(1|id)
            }else{
              zi_form<-~1
            }
          }else{
            stop("Invalid Zero-Infation covariate option") #No zero-inflation component)
          }
        }
      }
    }
  }

if(is.null(mean_form)){
  if(is.null(mean_covar)){
    if(mean_re==TRUE){
      mean_form<-count ~1+(1|id)} #Default is intercept only
    if(mean_re==FALSE){
      mean_form<-count ~1} #Default is intercept only
  }else{
    if(is.atomic(mean_covar)&length(mean_covar)==1){
      if(mean_covar==0){stop("At minimum an intercept is needed in the mean model")}else{
        if(mean_covar==1){
          if(mean_re==TRUE){
            mean_form<-count~1+(1|id)
          }else{
            mean_form<-count~1
          }
        }else{
          stop("Invalid Mean Covariate Option") #No zero-inflation component)
        }
      }
    }
  }
}

  if(is.null(disp_covar)){
      disp_form<- ~1 #Default is intercept only
  }else{
    if(is.atomic(disp_covar)&length(disp_covar)==1){
      if(disp_covar==0){stop("Invalid dispersion covariate option")}else{
        if(disp_covar==1){
            disp_form<-~1
        }else{
          stop("Invalid dispersion covariate option") #No zero-inflation component)
        }
      }
    }
  }

  if(is.null(disp_form)){
    disp_form<- ~disp_covar}


  if(is.null(mean_form)){
    if(mean_re==TRUE){
      mean_form<-count~mean_covar+ (1|id)
    }
    if(mean_re==FALSE){
      mean_form<-count~mean_covar
    }
  }
  if(is.null(zi_form)){
    if(zi_re==TRUE){
      zi_form<-~zi_covar+ (1|id)
    }
    if(zi_re==FALSE){
      zi_form<-~zi_covar
    }
  }
  mean_form<-as.formula(mean_form)
  zi_form<-as.formula(zi_form)
  return(list(mean_form=mean_form,zi_form=zi_form,disp_form=disp_form))

}
create_adhoc_formulas<-function(mean_covar,zi_covar){
  if(is.matrix(mean_covar)&is.matrix(zi_covar)){
    form<-count~mean_covar|zi_covar
  }
  if(class(zi_covar)=="numeric" & is.matrix(mean_covar)){
    if(zi_covar==0){
      stop("Ad hoc method only implemented when ZI model contains at minimum an intercept")
      }else if(zi_covar==1){
      form<-count~mean_covar|1
      }else{
      stop("Invalid zero-inflation covariates specified")
      }
  }
  if(class(mean_covar)=="numeric" & is.matrix(zi_covar)){
    if(mean_covar==0){
      stop("Ad hoc method only implemented when Mean model contains at minimum an intercept")
    }else if(mean_covar==1){
      form<-count~1|zi_covar
    }else{
      stop("Invalid mean model covariates specified")
    }
  }
  if(class(mean_covar)=="numeric" & class(zi_covar)=="numeric"){
    if(mean_covar==0 | zi_covar==0){
      stop("Ad hoc method only implemented when both mean and zi model contain at minimum an intercept")
    }
    if(mean_covar==1 & zi_covar==1){
      form<-count~1|1
    }
    else{
      stop("Ad hoc method only implemented when both mean and zi model contain at minimum an intercept")
    }
  }

  return(form)
}
