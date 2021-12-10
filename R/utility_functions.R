# Should we have NULL for covars if requiring their input?
# Should we extend ad hoc to non ZI case?
check_twosigma_input<-function(count,mean_covar,zi_covar
                               ,mean_re=TRUE,zi_re=TRUE
                               ,disp_covar=NULL,adhoc,id){

  # Override count with mean_form if specified then check the inputs
  if(length(unique(id))==1){
    if(adhoc==TRUE | (adhoc==FALSE & (mean_re==TRUE|zi_re==TRUE))){
      stop("The id variable only has one level. The adhoc method or any random effect inclusion is not appropriate.")
    }
  }
  if(any(is.na(count))){stop("Missing values not allowed in count matrix, please remove NA values.")}
  if(any(is.na(mean_covar))){stop("Missing values not allowed in mean covariate matrices, please remove NA values.")}
  if(any(is.na(zi_covar))){stop("Missing values not allowed in ZI covariate matrices, please remove NA values.")}
  if(sum(!as.matrix(count,ncol=1)%%1==0)>0 | min(count)<0){
    stop("When using the Negative Binomial Distribution data must contain only non-negative integers")
  }

  if(is.null(mean_covar))
  {
    stop("At minimum must specify Mean Model Covariates or provide model formula")
  }

  #Check for intercept being input
  # if(class(mean_covar)=="matrix"){
  #   check_no_intercept<-function(x){mean(x==1)}
  #   if(any(apply(mean_covar,2,check_no_intercept)==1)){
  #     stop("Please remove intercept from mean model covariate matrix")
  #   }
  # }
  # if(class(zi_covar)=="matrix"){
  #     check_no_intercept<-function(x){mean(x==1)}
  #     if(any(apply(zi_covar,2,check_no_intercept)==1)){
  #     stop("Please remove intercept from ZI model covariate matrix")
  #     }
  # }
  #     if(class(disp_covar)=="matrix"){
  #       check_no_intercept<-function(x){mean(x==1)}
  #       if(any(apply(disp_covar,2,check_no_intercept)==1)){
  #         stop("Please remove intercept from dispersion covariate matrix")
  #       }
  #     }

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
}

check_twosigma_custom_input<-function(count
  ,mean_form,zi_form
  ,mean_re=TRUE,zi_re=TRUE
  ,disp_covar=NULL){

  # Override count with mean_form if specified then check the inputs
  if(mean_form[[2]]!="count"){stop("Please begin the two-sided formula mean_form with the name 'count'. Failure to do so will cause an error downstream.")}
  if(length(zi_form)>2){stop("ZI Formula should be one-sided and thus only have a length of 2.")}
  if(any(is.na(count))){stop("Missing values not allowed in count matrix, please remove NA values.")}
  if(sum(!as.matrix(count,ncol=1)%%1==0)>0 | min(count)<0){
    stop("When using the Negative Binomial Distribution data must contain only non-negative integers")
  }

  # if(mean(count==0)>.9){
  #   warning("More than 90% of data are zeros. Mean model results may be misleading for such sparse data")
  # }
  # if(mean(count==0)<.1){
  #   if(!(zi_form==~0)){
  #     warning("Less than 10% of data are zeros. Zero-Inflation model results may be misleading or unnecessary")
  #   }
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
      #mean_form<-count~mean_covar+ (1|id)
      if(!is.null(colnames(mean_covar))){
        #mean_form<-as.formula(paste0("count","~",paste(unlist(strsplit(colnames(mean_covar),split=" ")),collapse="+"),"+","(1|id)"))
        mean_form<-count~mean_covar+ (1|id)
      }else{
        mean_form<-count~mean_covar+ (1|id)
      }
    }
    if(mean_re==FALSE){
      if(!is.null(colnames(mean_covar))){
      #mean_form<-as.formula(paste0("count","~",paste(unlist(strsplit(colnames(mean_covar),split=" ")),collapse="+")))
        mean_form<-count~mean_covar
      }else{
        mean_form<-count~mean_covar
      }
    }
  }
  if(is.null(zi_form)){
    if(zi_re==TRUE){
      if(!is.null(colnames(zi_covar))){
        #zi_form<-as.formula(paste0("~",paste(unlist(strsplit(colnames(zi_covar),split=" ")),collapse="+"),"+","(1|id)"))
        zi_form<-~zi_covar+ (1|id)
      }else{
        zi_form<-~zi_covar+ (1|id)
      }

    }
    if(zi_re==FALSE){
      if(!is.null(colnames(zi_covar))){
      #zi_form<-as.formula(paste0("~",paste(unlist(strsplit(colnames(zi_covar),split=" ")),collapse="+")))
        zi_form<-~zi_covar
      }else{
        zi_form<-~zi_covar
      }
    }
  }
  mean_form<-as.formula(mean_form)
  zi_form<-as.formula(zi_form)
  return(list(mean_form=mean_form,zi_form=zi_form,disp_form=disp_form))

}

create_adhoc_formulas2<-function(count,mean_covar,zi_covar){
  if(is.matrix(mean_covar)&is.matrix(zi_covar)){
    form<-count~mean_covar|zi_covar
  }
  if((class(zi_covar)=="numeric"| (class(zi_covar)=="integer")) & length(zi_covar)==1 & is.matrix(mean_covar)){
    if(length(zi_covar)==1){
      if(zi_covar==0){
        stop("Ad hoc method only implemented when ZI model contains at minimum an intercept")
      }else if(zi_covar==1){
        form<-count~mean_covar|1
      }else{
        stop("Invalid zero-inflation covariates specified")
      }
    }else{
      stop("Invalid zero-inflation covariates specified")
    }
    form<-count~mean_covar|zi_covar
  }
  if((class(mean_covar)=="numeric"| (class(mean_covar)=="integer")) & length(mean_covar)==1 & is.matrix(zi_covar)){
    if(length(mean_covar==1)){
      if(mean_covar==0){
        stop("Ad hoc method only implemented when Mean model contains at minimum an intercept")
      }else if(mean_covar==1){
        form<-count~1|zi_covar
      }else{
        stop("Invalid mean model covariates specified")
      }
    }else{
      stop("Invalid mean model covariates specified")
    }
    form<-count~mean_covar|zi_covar
  }
  if((class(mean_covar)=="numeric"| (class(mean_covar)=="integer")) & length(mean_covar)>1 & is.matrix(zi_covar)){
    form<-count~mean_covar|zi_covar
  }

  if((class(mean_covar)=="numeric"| (class(mean_covar)=="integer")) & length(mean_covar)>1
    & (class(zi_covar)=="numeric"| (class(zi_covar)=="integer")) & length(zi_covar)>1){
      if(length(mean_covar)==1 & length(zi_covar)==1){
        if(mean_covar==0 | zi_covar==0){
          stop("Ad hoc method only implemented when both mean and zi model contain at minimum an intercept")
        }
        if(mean_covar==1 & zi_covar==1){
          form<-count~1|1
        }
        else{
          stop("Ad hoc method only implemented when both mean and zi model contain at minimum an intercept")
        }
      }else{
        form<-count~mean_covar|zi_covar
      }

    }
  if((class(zi_covar)=="numeric"| (class(zi_covar)=="integer")) & length(zi_covar)>1 & is.matrix(mean_covar)){
    form<-count~mean_covar|zi_covar
  }
  return(form)
}

create_adhoc_formulas<-function(count,mean_covar,zi_covar){
  if(length(zi_covar)==1){
    if(zi_covar==0){
      stop("Ad hoc method only implemented when both mean and zi model contain at minimum an intercept")}
    if(zi_covar!=1){
      stop("Invalid zi model covariates specified")
    }
  }
  if(length(mean_covar)==1){
    if(mean_covar==0){
      stop("Ad hoc method only implemented when both mean and zi model contain at minimum an intercept")}
    if(mean_covar!=1){
      stop("Invalid mean model covariates specified")
    }
  }
  if(length(mean_covar)==1 & length(zi_covar)>1){
    if(mean_covar==1){form<-count~1|zi_covar}
  }
  if(length(mean_covar)>1 & length(zi_covar)==1){
    if(zi_covar==1){form<-count~mean_covar|1}
  }
  if(length(mean_covar)==1 & length(zi_covar)==1){
    if(mean_covar==1 & zi_covar==1){form<-count~1|1}
  }
  if(!exists("form")){
    form<-count~mean_covar|zi_covar
  }
  return(form)
}
glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
    coef. = function(x) fixef(x)[[component]],
    vcov. = function(x) vcov(x)[[component]],
    df = NULL)
}

# This is a useful, CRAN-friendly way to access an unexported function
# which is required in the twosigmag function
modelparm.default<-utils::getFromNamespace("modelparm.default","multcomp")

modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
  modelparm.default(model, coef. = coef., vcov. = vcov.,
                    df = df, ...)
}

# modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
#   vcov. = function(x) vcov(x)[[component]],
#   df = NULL, component="cond", ...) {
#   multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
#     df = df, ...)
# }

