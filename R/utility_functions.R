
check_twosigma_input<-function(count,mean_covar=NULL,zi_covar=NULL
                               ,mean_form=NULL,zi_form=NULL
                               ,mean_re=TRUE,zi_re=TRUE
                               ,disp_covar=NULL){

  if(sum(!as.matrix(count,ncol=1)%%1==0)>0 | min(count)<0){
    stop("When using the Negative Binomial Distribution data must contain only non-negative integers")
  }

  if(is.null(mean_covar)&is.null(mean_form))
  {
    stop("At minimum must specify Mean Model Covariates or provide model formula")
  }

}
create_model_formulas<-function(mean_covar=NULL,zi_covar=NULL
                                ,mean_form=NULL,zi_form=NULL
                                ,mean_re=TRUE,zi_re=TRUE
                                ,disp_covar=NULL){

  if(is.null(zi_covar)){
    if(zi_re==TRUE){
      zi_form<- ~1+(1|id)} #Default is intercept only
    if(zi_re==FALSE){
      zi_form<- ~1} #Default is intercept only
  }else{
    if(class(zi_covar)=="numeric"){
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


  if(is.null(mean_covar)){
    if(mean_re==TRUE){
      mean_form<-count ~1+(1|id)} #Default is intercept only
    if(mean_re==FALSE){
      mean_form<-count ~1} #Default is intercept only
  }else{
    if(class(mean_covar)=="numeric"){
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


  if(is.null(disp_covar)){
    disp_form<- ~1}else{
      disp_form<-~disp_covar
    }
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
