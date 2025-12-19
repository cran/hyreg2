

#' M-step driver to be used in flexmix
#'
#' @description Function used in flexmix M-Step to estimate hybrid model
#'
#' @param formula model `formula`, automatically provided by `hyreg2` and [flexmix::flexmix()]
#' @param family default `"hyreg"`, needed for [flexmix::flexmix()]
#' @param type `character` vector containing the indicator whether that datapoint (row)
#'  contains continuous or dichotomous data, see Details of [hyreg2]
#' @param type_cont value of `type` referring to continuous `data`, see Details of [hyreg2]
#' @param type_dich value of `type` referring to dichotomous `data`, see Details of [hyreg2]
#' @param variables_both `character` vector; variables to be fitted on both continuous and dichotomous data.
#'   see Details of [hyreg2]
#' @param variables_cont `character` vector; variables to be fitted only on continous data. see Details of [hyreg2]
#' @param variables_dich character vactor; variables to be fitted only on dichotomous data. see Details of [hyreg2]
#' @param stv `named vector` or `list` of named vectors containing start values for all coefficients from
#'  `formula`, including `theta`, see Details of [hyreg2::hyreg2]
#' @param offset offset as in [flexmix::flexmix()], default `NULL`
#' @param optimizer `character`, optimizer to be used in [bbmle::mle2()], default `"optim"`
#' @param opt_method `character`, optimization method to be used in `optimizer`, default `"BFGS"`
#' @param lower  lower bound for censored data. If this is used, `opt_method` must be
#'  set to `"L-BFGS-B"`, default `-INF`,
#' @param upper  upper bound for censored data. If this is used, `opt_method` must be
#'  set to `"L-BFGS-B"`,default `INF`
#' @param formula_type_classic `logical`; is the provided `formula` a typical R formula containing only variables
#'   or does it include variables and parameters? default `TRUE`
#' @param ... additional arguments for [flexmix::flexmix()] or [bbmle::mle2()]
#'
#' @return a `model` object, that can be used in [hyreg2] as input for parameter `model` in [flexmix::flexmix()]
#'
#'
#' @return a model object, that can be used in hyreg2 as input for parameter model in flexmix::flexmix
#'
#' @author Svenja Elkenkamp and Kim Rand
#'
#' @importFrom methods new
#' @importFrom stats as.formula dnorm model.matrix pnorm setNames
#'
#' @examples
#'
#'formula <- y ~  -1 + x1 + x2 + x3
#'the$k <- 2
#'stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'
#'
#'
#' x <- model.matrix(formula,simulated_data_norm)
#' y <- simulated_data_norm$y
#' w <- 1

#'model <- FLXMRhyreg(formula = formula,
#'                     family=c("hyreg"),
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = list(iter.max = 1000, verbose = 4),
#'                     offset = NULL,
#'                     optimizer = "optim",
#'                     variables_both =  names(stv)[!is.element(names(stv),c("sigma","theta"))],
#'                     variables_cont = NULL,
#'                     variables_dich = NULL,
#'                     lower = -Inf,
#'                     upper = Inf,
#')


#' @importFrom flexmix flexmix
#' @importFrom  bbmle mle2
#' @importFrom  bbmle summary
#' @export




### FLXMRhyreg ###


FLXMRhyreg <- function(formula= . ~ . ,
                       family=c("hyreg"),
                       type = NULL,
                       type_cont = NULL,
                       type_dich = NULL,
                       variables_both = NULL,
                       variables_cont = NULL,
                       variables_dich = NULL,
                       stv = NULL, # has to include starting values for sigma and theta as well
                       offset = NULL,
                       opt_method = "BFGS",
                       optimizer = "optim",
                       lower = -Inf,
                       upper = Inf,
                       formula_type_classic = TRUE,
                       ...
)
{
  family <- match.arg(family)

  # refit function has to depend on x,y,w.
  hyregrefit <- function(x, y, w) {
    warning(paste0("Not defined"))
    return(NA)
  }

  z <- new("FLXMRglm", weighted=TRUE, formula=formula, #  change formula for formula_type_classic?
           name=paste("FLXMRhyreg"), offset = offset,
           family="hyreg", refit=hyregrefit)

  z@preproc.y <- function(x){
    if (ncol(x) > 1)
      stop(paste("for the", family, "family y must be univariate"))
    x
  }

  if(family=="hyreg"){
    z@defineComponent <- function(para) {  # we get para from the first estimation with start values

      ### NEW ###
      predict <- function(x, type, ...) {
        dotarg = list(...)
        if("offset" %in% names(dotarg)) offset <- dotarg$offset


        if(formula_type_classic == TRUE){
          # CLASSIC
          if(type == type_cont){
            p <- x %*% para$coef[is.element(names(para$coef),c(variables_cont,variables_both))]  # Xb in xreg
          }
          if(type == type_dich){
            p <- (x %*% para$coef[is.element(names(para$coef),c(variables_dich,variables_both))]) * para$theta
          }
        }else{

          # NON-CLASSIC
          if(type == type_cont){
            p <- eval_formula_non(the$formula_non,the$data,para$coef[is.element(names(para$coef),c(variables_cont,variables_both))])
          }
          if(type == type_dich){
            p <- eval_formula_non(the$formula_non,the$data,para$coef[is.element(names(para$coef),c(variables_dich,variables_both))]) * para$theta
          }
        }


        if (!is.null(offset)) p <-  p + offset
        p
      }



      logLik <- function(x, y, sigma = para$sigma, theta = para$theta, return_vector = TRUE, ...){

        # prepare data
        # choose subset of x and y depending on type

        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]


        sigma <- exp(sigma)

        if(isTRUE(formula_type_classic)){

          # classic
          x1 <- x[type == type_cont,c(variables_cont,variables_both)]
          x2 <- x[type == type_dich,c(variables_dich,variables_both)]

          if(length(c(variables_cont,variables_both)) == 1){
            Xb1 <- as.matrix(x1 * para$coef[c(variables_cont,variables_both)])
            colnames(Xb1) <- c(variables_cont,variables_both)
          }else{
            Xb1 <- x1 %*% para$coef[colnames(x1)]
          }

          if(length(c(variables_dich,variables_both)) == 1){
            Xb2 <- as.matrix( (x2 * para$coef[c(variables_dich,variables_both)]) * exp(theta) )
            colnames(Xb2) <- c(variables_dich,variables_both)
          }else{
            Xb2 <-(x2 %*% para$coef[colnames(x2)]) * exp(theta)
          }


          #  Xb2 <- (x2[variables_both] %*% para$coef[variables_both]) * exp(theta) +
          #         (x2[variables_dich] %*% para$coef[variables_dich])  # theta only for variables_both


        }else{

          # non-classic
          Xb1 <- as.matrix(eval_formula_non(the$formula_cont, # formula_non only for cont
                                            the$data[type == type_cont,],
                                            para$coef[c(variables_cont,variables_both)]))
          Xb2 <- as.matrix(eval_formula_non(the$formula_dich, # formula_non only for dich
                                            the$data[type == type_dich,],
                                            para$coef[c(variables_dich,variables_both)])) * exp(theta)  # only dich and both variables

        }

        # pvals and likelihood
        logistic_tmp <- .5+.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic

        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal


        # for censored data
        if(upper != Inf ){
          censV <- y1 == upper
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-upper)/sigma,0,1) - pnorm((Xb1[censV]-Inf)/sigma,0,1))
        }

        if(lower != -Inf ){
          censV <- y1 == lower
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-(-Inf))/sigma,0,1) - pnorm((Xb1[censV]-lower)/sigma,0,1))
        }


        pvals <- c(pvals1,pvals2)

        pvals[pvals == -Inf] <- log(.Machine$double.xmin)
        pvals[pvals == Inf] <- log(.Machine$double.xmax)

        if(return_vector == TRUE){
          return(pvals)
        }else{
          return(-sum(pvals))         # get neg log L for optimizer
        }

      }

      new("FLXcomponent",
          parameters=list(coef=para$coef,
                          sigma=para$sigma,
                          theta = para$theta,
                          fit_mle = para$fit_mle),
          logLik=logLik, predict=predict,
          df=para$df)
    }


    z@fit <- function(x, y, w, component, ...){


      ### LIKELIHOOD FUNCTION to be used in ML Estimation ###

      # function to be used in mle, same as logLik but depending on stv and giving out the neg logL directly
      logLik2 <- function(stv){

        # prepare data

        y1 <- y[type == type_cont]
        y2 <-  y[type == type_dich]

        # prepare stv and sigma and theta
        sigma <- exp(stv[is.element(names(stv),c("sigma"))][[1]])
        theta <- exp(stv[is.element(names(stv),c("theta"))][[1]])
        stv_cont <- stv[!is.element(names(stv),c("sigma","theta", variables_dich))]
        stv_dich <- stv[!is.element(names(stv),c("sigma","theta", variables_cont))]


        if(isTRUE(formula_type_classic)){

          # classic
          x1 <- x[type == type_cont,c(variables_cont,variables_both)]
          x2 <-  x[type == type_dich,c(variables_dich,variables_both)]

          if(length(stv_cont) == 1){
            Xb1 <- as.matrix(x1 * stv_cont)
            colnames(Xb1) <- stv_cont
          }else{
             Xb1 <- x1 %*% stv_cont[colnames(x1)]
          }

          if(length(stv_dich) == 1){
            Xb2 <- as.matrix(x2 * stv_dich)
            colnames(Xb2) <- stv_dich
          }else{
            Xb2 <- x2 %*% stv_dich[colnames(x2)]
          }


          Xb2 <- Xb2*theta


          # for variables_dich use theta only with variables_both?
          # Xb2 <- (x2[variables_both] %*% stv_dich[variables_both]) * theta +
          #  (x2[variables_dich] %*% stv_dich[variables_dich])

        }else{

          # non-classic
          Xb1 <- as.matrix(eval_formula_non(the$formula_cont, # formula_non only for cont
                                            the$data[type == type_cont,],
                                            stv_cont))
          Xb2 <- as.matrix(eval_formula_non(the$formula_dich, # formula_non only for dich
                                            the$data[type == type_dich,],
                                            stv_dich))

          Xb2 <- Xb2*theta

        }

        # pvals and likelihood
        logistic_tmp <- 0.5 + 0.5*tanh(Xb2/2)
        pvals2 <- log(y2 *logistic_tmp + (1-y2)* (1-logistic_tmp)) # dich_logistic

        pvals1 <- dnorm(y1, mean=Xb1, sd=sigma, log=TRUE) # cont_normal


        # for box constraints, censored data
        if(upper != Inf ){
          censV <- y1 == upper
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-upper)/sigma,0,1) - pnorm((Xb1[censV]-Inf)/sigma,0,1))
        }
        if(lower != -Inf ){
          censV <- y1 == lower
          pvals1[which(censV)] <- log(pnorm((Xb1[censV]-(-Inf))/sigma,0,1) - pnorm((Xb1[censV]-lower)/sigma,0,1))
        }


        pvals <- c(pvals1,pvals2)


        pvals[pvals == -Inf] <- log(.Machine$double.xmin)
        pvals[pvals == Inf] <- log(.Machine$double.xmax)



        # multiply weights for use in EM algo
        pvals_w <- pvals * w


        return(-sum(pvals_w))
      }

      # set names of inputs for logLik2
      if(isTRUE(formula_type_classic)){
        bbmle::parnames(logLik2) <- c(colnames(x),"sigma","theta")
      }else{
        if(is.list(stv)){
          bbmle::parnames(logLik2) <- c(names(the$stv[[1]]))
        }else{
          bbmle::parnames(logLik2) <- c(names(the$stv))
        }
      }

      if(!exists("counter", envir = the)){

        the$counter <- 1


        if(is.list(stv)){
          stv_in <- stv[[the$counter]]
        }else{
          stv_in <- stv
        }

        fit_mle <- bbmle::mle2(minuslogl = logLik2,
                               start = stv_in,
                               optimizer = optimizer,
                               method = opt_method,
                               # control?,
                               lower = -Inf, # or upper and lower from input??
                               upper = Inf)



      }else{
        if(the$counter < the$k){
          the$counter <- the$counter + 1

          if(is.list(stv)){
            stv_in <- stv[[the$counter]]
          }else{
            stv_in <- stv
          }

          fit_mle <- bbmle::mle2(minuslogl = logLik2,
                                 start = stv_in,
                                 optimizer = optimizer,
                                 method = opt_method,
                                 # control?,
                                 lower = -Inf, # or upper and lower from input?? # depending on method?
                                 upper = Inf)



        }else{
          if(isTRUE(formula_type_classic)){
            stv_new <- setNames(c(component$coef,component$sigma,component$theta),c(colnames(x),"sigma","theta"))
          }else{
            if(is.list(stv)){
              stv_new <- setNames(c(component$coef,component$sigma,component$theta),c(names(the$stv[[1]]))) # order important, maybe flexiblize?
            }else{
              stv_new <- setNames(c(component$coef,component$sigma,component$theta),c(names(the$stv))) # order important, maybe flexiblize?
            }
          }

          fit_mle <- bbmle::mle2(minuslogl = logLik2,
                                 start = stv_new,
                                 optimizer = optimizer,
                                 method = opt_method,
                                 # control?,
                                 lower = -Inf, # or upper and lower from input??
                                 upper = Inf)
        }
      }


      z@defineComponent(para = list(coef = fit_mle@coef[!is.element(names(fit_mle@coef),c("sigma","theta"))],
                                    df = ncol(x)+1, # not changed yet
                                    sigma = fit_mle@coef[is.element(names(fit_mle@coef),c("sigma"))],
                                    theta = fit_mle@coef[is.element(names(fit_mle@coef),c("theta"))],
                                    fit_mle = fit_mle,
                                    minLik = fit_mle@min)
      )
    }
  }

  z
}

