

#' Estimating hybrid models
#'
#' @description Estimation of hybrid model using continuous and dichotomous data e.g. EQ-5D data
#'
#' @param formula model `formula`, can be linear or non-linear. For non-linear formulas, variables and parameters
#' must be provided and `formula_type_classic` must be set to `FALSE`. Using `|xg` will include a grouping variable `xg`. see Details.
#' @param data a `data.frame` containing the data. see Details.
#' @param type either the name of the column in `data` containing an indicator of whether an
#'   observation is continuous or dichotomous (as `character`), or a `vector` containing the indicator.
#'   see Details.
#' @param type_cont value of `type` referring to continuous data. see Details.
#' @param type_dich Value of `type` referring to dichotomous data. see Details.
#' @param k `numeric`. Number of latent classes to be estimated via [flexmix::flexmix()].
#' @param control control list for [flexmix::flexmix()].
#' @param stv `named vector` or `list` of named vectors containing start values for all coefficients
#'  formula, including sigma and theta, see Details
#' @param offset offset as in [flexmix::flexmix()]
#' @param optimizer `charachter`,optimizer to be used in [bbmle::mle2()], default `"optim"`
#' @param opt_method `charachter`, optimization method to be used in optimizer, default `"BFGS"`
#' @param lower `numeric`, lower bound for censored data, default `-INF`. If this is used, `opt_method` must be set to `"L-BFGS-B"`,
#' @param upper `numeric`, upper bound for censored data, default `INF`. If this is used, `opt_method` must be set to `"L-BFGS-B"`,
#' @param latent `charachter`,data type to use in component identification, must be one of `"both"`, `"cont"` or `"dich"`,
#'  default  `“both”`,  see Details
#' @param id_col `character`, name of the grouping variable, only needed if `latent != "both”`, see Details
#' @param classes_only `logical`, default `FALSE`, indicates whether the function should perform only
#'  classification, rather than both classification and model estimation, only possible for `latent != "both"`,
#'   see Datails
#' @param variables_both `character vector`; variables to be fitted on both continuous and dichotomous data.
#'  If not specified, all variables from `formula` are used. If provided and not all variables from `formula`
#'  are included, `variables_cont` and `variables_dich` must be provided as well, while one of them can be `NULL`,
#'   see Details.
#' @param variables_cont `character` vector; variables to be fitted only on continous data. If provided,
#' `variables_both` and `variables_dich` must be provided as well.
#' @param variables_dich `character` vector; variables to be fitted only on dichotomous data, if provided,
#'  `variables_both` and `variables_cont` must be provided as well.
#' @param formula_type_classic `logical`; is the provided `formula` a classic R formula containing only variables (`TRUE`)
#'   or does it include both variables and parameters (`FALSE`)? default `TRUE`, see Details
#' @param ... additional arguments for [flexmix::flexmix()] or [bbmle::mle2()]
#'
#' @return model object of type `flexmix` or `list` of model objects of type `flexmix`.
#' Please note, that the estimates for `sigma` and `theta` are on a log-scale and have to be transformed using `exp()`to get the correct estimated values.
#'
#' @details
#' see details of different inputs listed below.
#'
#'
#'@section formula:
#' a classic R formula containing only variables(e.g.`y ~ x1 + x2 + …`) can be provided as well as a formula
#' including variables and parameters (non-classic) e.g. `y ~ x1 * beta1 + x2 * beta2`  or `y ~ 1/exp(x1 * beta1 + x2 * beta2)`,
#'  where `beta` are the parameters to be estimated and the`x`s are column names from the dataset.
#'  Non-linear models and the 8-parameter model for EQ-5D data can only be estimated using a non-classic formula.
#'  If the provided formula is non-classic, `formula_type_classic` must be set to `FALSE`.
#' When estimating an intercept, the formula must explicitly include a parameter named `"INTERCEPT"`(without a corresponding variable from the dataset)
#' Additionally, it is possible to include a grouping variable for repeated measures by using
#' `“| xg”` where `xg` is the column containing the group-memberships. The resulting formula will look
#' like this:  `y ~ x1 + x2 +… | xg`.  In `flexmix`, this is called the concomitant variable specification:
#' the  model is fit conditional on grouping, so that all observations with the same group are treated
#' as belonging together when computing likelihood contributions. One possible grouping variable can be
#' an id number to identify answers by the same participants. We highly recommend using a grouping variable,
#'  since otherwise the algorithm for k = 2 tends to classify all continuous data into one estimated class
#' and all dichotomous data into the other.
#'
#'
#' @section data:
#'  a dataframe having the following columns: all independent variables (x)
#'  and the dependent variable y used in `formula`, one column for the grouping variable xg if grouping
#'  should be used, e.g. id numbers of participants with repeated measurements, one column indicating
#'  if the observations belongs to continuous or dichotomous data with the entries `type_cont`
#'  and `type_dich` (e.g., for a column called `"type"` with the entries "TTO" for continuous datapoints
#'  and "DCE" for dichotomous datapoints, `type_cont` will be "TTO" and `type_dich` will be "DCE").
#'  One row should match one observation (one datapoint).
#'
#'
#' @section start values (stv):
#' if the same start values `stv` are to be used for all latent classes,
#' the given start values must be a `named vector`. Otherwise (if different start values are assumed for
#'  each latent class), a `list` of named vectors should be used . In this case, there must be one entry
#'  in the list for each latent class.  Each start value vector must include start values for sigma and
#'  theta. Currently, it is necessary to use the names `"sigma"` and `"theta"` for these values.
#'  If users are unsure for which variables start values must be provided (in the linear formula case), this can be checked by
#'  calling `colnames(model.matrix(formula,data))`. In this call, the `formula` should not include the
#'  grouping variable.
#'
#'
#' @section latent, id_col, classes_only:
#' in some situations, it can be useful to identify the latent classes on
#' only one `type` of data while estimating the model parameters on both `types` of data. In such cases,
#' the input variable `latent` can be used to specify on which type of data the classification should be done.
#'  If `“cont”` or `“dich”` is used, `formula` must contain a grouping variable and additionally the
#'  input parameter `id_col` must be specified and gives the name,
#'  i.e. a `character string`, of the grouping variable for classification. Some groups may be removed from
#'  the data, since they have only continuous or only dichotomous observations. Then in a first step,
#'  a model is estimated only on the continuous/dichotomous data and the achieved classification is stored.
#'  In a next step, model parameters are estimated separately for each identified class on both `types` of data
#'  using this classification. The output object of `hyreg2` in this case is a `list` of k models.
#'  Additionally, at position k+1 of the list, a data frame containing the corresponding classifications
#'  from the first step is returned. Each element k in the list contains the estimated parameters for one
#'  of the latent classes. When setting the input variable `classes_only` to `TRUE`, the second step is left
#'  out and the estimated classes from step one are given as output.
#'
#'
#' @section variables_both, variables_cont, variables_dich:
#' It is possible to specify partial coefficients,
#'  which are used only on continuous or dichotomous data.
#'*  Example:  Suppose different models should be specified for continuous and dichotomous  data:
#'*  Model continuous data: `y ~  x1 + x3`
#'*  Model dichotomous data: `y ~  x1 + x2`
#'*  The `formula` input to `hyreg2` must then include all parameters that occur in either model:
#'   `y ~ x1 + x2 + x3`
#'* The assignment of parameters to data types is then achieved via the input arguments `variables_both`,
#'   `variables_cont`, and `variables_dich`:
#'*  `variables_both` = `“x1”`,
#'*  `variables_cont` = `“x3”` and
#'* `variables_dich` = `“x2”`.
#'* Every variable included in the provided `formula` (except the grouping variable ) must appear in exactly
#'   one of these vectors. One of the `variables_` vectors can also be `NULL`, if no variables should be used only on this type of the data.

#' @author Svenja Elkenkamp, Kim Rand and John Grosser
#' @examples
#'
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'
#'k <- 2

#'stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control = list(iter.max = 1000, verbose = 4)

#'rm(counter)

#'mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type, # also "type" would work
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "cont",
#'                     id_col = "id"
#')
#'summary_hyreg2(mod)

#' @importFrom flexmix flexmix
#' @importFrom bbmle mle2
#' @export


hyreg2 <-function(formula,
                  data,
                  type,
                  type_cont,
                  type_dich,
                  k = 1,
                  control = NULL,
                  stv = NULL, # has to include starting values for sigma and theta as well
                  offset = NULL,
                  opt_method = "BFGS",
                  optimizer = "optim",
                  lower = -Inf,
                  upper = Inf,
                  latent = "both", # one of "both", "cont", "dich"
                  id_col = NULL, # as character
                  classes_only = FALSE,
                  variables_both = NULL,
                  variables_dich = NULL,
                  variables_cont = NULL,
                  formula_type_classic = TRUE,
                  # additional arguments for flexmix or optimizer ?
                  ...){

  dotarg <- list(...)

  if(is.character(type) & length(type) == 1 ){
    type <- data[,type]
  }

  # assign k to k in package environment to be used during M-step driver
  assign("k", k, envir=the)
  if(exists("counter", envir = the)){
    rm("counter", envir = the)
  }

  # prepare formula handling
  formula_string <- paste(deparse(formula), collapse = "")
  formula_parts <- strsplit(formula_string, "\\|")[[1]]
  formula_short <- as.formula(formula_parts[1])




  # STV CHECK BASED ON FORMULA_TYPE_CLASSIC
  if(isTRUE(formula_type_classic)){

    ### STV Check Classic ###
    # check stv for names in x (model.matrix(formula,data))

    # no stv
    if(!is.list(stv)){
      if(is.null(stv)){
        warning(paste0("Argument stv not provided. Setting all start values from formula to 0.1 and start values for sigma and theta to 1."))
        stv <- setNames(c(rep(0.1,dim(model.matrix(formula_short,data))[2]),1,1), c(colnames(model.matrix(formula_short,data)),c("sigma","theta")))
      }else{

        # one or more stv missing
        if(any(!is.element(colnames(model.matrix(formula_short,data)),names(stv)))){
          miss <- colnames(model.matrix(formula_short,data))[!is.element(colnames(model.matrix(formula_short,data)),names(stv))]
          stop(paste0("Start value(s) missing for ", paste(miss, collapse = ", ") ,". Please provide start values for all relevant formula variables."
          ))
        }

        # theta missing
        if(!is.element("theta",names(stv))){
          stv <- c(stv,setNames(1,"theta"))
          warning(paste0("Start value missing for theta, setting to 1."))
        }

        # sigma missing
        if(!is.element("sigma",names(stv))){
          stv <- c(stv,setNames(1,"sigma"))
          warning(paste0("Start value missing for sigma, setting to 1."))
        }

        # stv for variables not in formula given, d.h. zu viele angegeben
        if(any(!is.element(names(stv), c(colnames(model.matrix(formula_short,data)),"theta","sigma")))){
          much <- names(stv)[!is.element(names(stv), c(colnames(model.matrix(formula_short,data)),"theta","sigma"))]
          stop(paste0("Start values provided for variables not in formula: ",paste(much, collapse = ", ")))
        }
        #  check order in FLXMRhyreg
      }

    }else{ # stv is a list

      # warning(paste0("stv is a list. Please ensure that variables_both, variables_dich and variables_cont are provided"))

      for(i in 1:length(stv)){
        # one or more stv missing
        if(any(!is.element(colnames(model.matrix(formula_short,data)),names(stv[[i]])))){
          miss <- colnames(model.matrix(formula_short,data))[!is.element(colnames(model.matrix(formula_short,data)),names(stv[[i]]))]
          stop(paste0("Start value(s) missing for ", paste(miss, collapse = ", ") ,". Please provide start values for all relevant formula variables."
          ))
        }

        # theta missing
        if(!is.element("theta",names(stv[[i]]))){
          stv[[i]] <- c(stv[[i]],setNames(1,"theta"))
          warning(paste0("Start value missing for theta, setting to 1."))
        }

        # sigma missing
        if(!is.element("sigma",names(stv[[i]]))){
          stv[[i]] <- c(stv[[i]],setNames(1,"sigma"))
          warning(paste0("Start value missing for sigma, setting to 1."))
        }

        # stv for variables not in formula given, d.h. zu viele angegeben
        if(any(!is.element(names(stv[[i]]), c(colnames(model.matrix(formula_short,data)),"theta","sigma")))){
          much <- names(stv[[i]])[!is.element(names(stv[[i]]), c(colnames(model.matrix(formula_short,data)),"theta","sigma"))]
          stop(paste0("Start values provided for variables not in formula: ",paste(much, collapse = ", ")))
        }
      }
    }
  }else{  # close formula_type_classic TRUE

    # ADAPTIONS FOR NON-CLASSIC FORMULAS

    # save data and formula in the environment
    assign("formula_non", formula_short, envir=the)
    assign("data", data, envir=the)

    # transform formulas for flexmix inputs
    formula <- get_data_vars(formula, data)
    formula_short <- get_data_vars(formula_short, data)


    ### STV Check NON-CLASSIC ###

    # check stv for names in formula
    vars <- all.vars(the$formula_non[[3]])[!is.element(all.vars(the$formula_non[[3]]),names(data))]

    # no stv
    if(!is.list(stv)){
      if(is.null(stv)){
        warning(paste0("Argument stv not provided. Setting all start values from formula to 0.1 and start values for sigma and theta to 1."))
        stv <- setNames(c(rep(0.1,length(vars)),1,1), c(vars,c("sigma","theta")))
      }else{

        # one or more stv missing
        if(any(!is.element(vars,names(stv)))){
          miss <- vars[!is.element(vars,names(stv))]
          stop(paste0("Start value(s) missing for ", paste(miss, collapse = ", ") ,". Please provide start values for all relevant formula variables."
          ))
        }

        # theta missing
        if(!is.element("theta",names(stv))){
          stv <- c(stv,setNames(1,"theta"))
          warning(paste0("Start value missing for theta, setting to 1."))
        }

        # sigma missing
        if(!is.element("sigma",names(stv))){
          stv <- c(stv,setNames(1,"sigma"))
          warning(paste0("Start value missing for sigma, setting to 1."))
        }

        # stv for variables not in formula given, d.h. zu viele angegeben
        if(any(!is.element(names(stv), c(vars,"theta","sigma")))){
          much <- names(stv)[!is.element(names(stv), c(vars,"theta","sigma"))]
          stop(paste0("Start values provided for variables not in formula: ",paste(much, collapse = ", ")))
        }
        #  check order in FLXMRhyreg
      }

    }else{ # stv is a list

      # warning(paste0("stv is a list. Please ensure that variables_both, variables_dich and variables_cont are provided"))

      for(i in 1:length(stv)){
        # one or more stv missing
        if(any(!is.element(vars,names(stv[[i]])))){
          miss <- vars[!is.element(vars,names(stv[[i]]))]
          stop(paste0("Start value(s) missing for ", paste(miss, collapse = ", ") ,". Please provide start values for all relevant formula variables."
          ))
        }

        # theta missing
        if(!is.element("theta",names(stv[[i]]))){
          stv[[i]] <- c(stv[[i]],setNames(1,"theta"))
          warning(paste0("Start value missing for theta, setting to 1."))
        }

        # sigma missing
        if(!is.element("sigma",names(stv[[i]]))){
          stv[[i]] <- c(stv[[i]],setNames(1,"sigma"))
          warning(paste0("Start value missing for sigma, setting to 1."))
        }

        # stv for variables not in formula given, d.h. zu viele angegeben
        if(any(!is.element(names(stv[[i]]), c(vars,"theta","sigma")))){
          much <- names(stv[[i]])[!is.element(names(stv[[i]]), c(vars,"theta","sigma"))]
          stop(paste0("Start values provided for variables not in formula: ",paste(much, collapse = ", ")))
        }
      }
    }
  } # close formula_type_classic FALSE check (non-classic formulas)



  ### TYPE Check ###
  if(is.null(type) | is.null(type_dich) | is.null(type_cont)){
    stop(paste0("Argument(s) type, type_dich and/or typ_cont not provided."))
  }else{
    if(!is.element(type_dich,unique(type))){
      warning(paste0("Invalid type_dich: provided name not present in values of type"))
    }
    if(!is.element(type_cont,unique(type))){
      warning(paste0("Invalid type_cont: provided name not present in values of type"))
    }
  }



  ### VARIABALES Check ###
  if(!is.list(stv)){
    if(is.null(variables_both) & is.null(variables_dich) & is.null(variables_cont)){
      variables_both <- names(stv)[!is.element(names(stv),c("sigma","theta"))]
    }else{
      if(any(!is.element(names(stv)[!is.element(names(stv),c("sigma","theta"))],
                         c(variables_both,variables_dich,variables_cont)))){
        # check if all variables are included
        stop(paste0("All variables named in stv must be contained in exactly one of: (1) variables_both, (2) variables_dich or (3) variables_cont."))
        # alternative: do not provide any of the vectors
        # than all relevant variables are set to variables_both automatically
      }
      if(any(table(c(variables_both,variables_cont,variables_dich))>1)){
        stop(paste0("All variables named in stv must be contained in exactly one of: (1) variables_both, (2) variables_dich or (3) variables_cont."))
      }
    }

  }else{

    # if stv is a list, both classes have to depend on the same set of variables,
    # stv of different classes have same names (hence using stv[[1]] is okay here)
    # different set of variables for each class not supported yet!

    if(is.null(variables_both) & is.null(variables_dich) & is.null(variables_cont)){
      variables_both <- names(stv[[1]])[!is.element(names(stv[[1]]),c("sigma","theta"))]
    }else{
      if(any(!is.element(names(stv[[1]])[!is.element(names(stv[[1]]),c("sigma","theta"))],
                         c(variables_both,variables_dich,variables_cont)))){
        # check if all variables are included
        stop(paste0("All variables named in stv must be contained in exactly one of: (1) variables_both, (2) variables_dich or (3) variables_cont."))
        # alternative: do not provide any of the vectors
        # than all relevant variables are set to variables_both automatically
      }
      if(any(table(c(variables_both,variables_cont,variables_dich))>1)){
        stop(paste0("All variables named in stv must be contained in exactly one of: (1) variables_both, (2) variables_dich or (3) variables_cont."))
      }
    }
  }


  ### PREPARE M-STEP FOR NON-CLASSIC FORMULAS ###

  if(isFALSE(formula_type_classic)){

    # hard coding for stv as list because formulas must be the same for all components
    # change this if formulas can differ between components

    assign("stv", stv, envir=the)

    if(!is.list(stv)){
      assign("stv_cont", stv[!is.element(names(stv),c("sigma","theta", variables_dich))], envir = the)
      assign("stv_dich", stv[!is.element(names(stv),c("sigma","theta", variables_cont))], envir = the)
    }else{
      assign("stv_cont", stv[[1]][!is.element(names(stv[[1]]),c("sigma","theta", variables_dich))], envir = the)
      assign("stv_dich", stv[[1]][!is.element(names(stv[[1]]),c("sigma","theta", variables_cont))], envir = the)
    }

    assign("formula_cont", formulacheck_variables(the$formula_non, the$data, the$stv_cont), envir = the)
    assign("formula_dich", formulacheck_variables(the$formula_non, the$data, the$stv_dich), envir = the)
  }


  ### ESTIMATION ###
  if(latent == "both"){

    model <- list(FLXMRhyreg(type= type,
                             stv = stv,
                             type_cont = type_cont,
                             type_dich = type_dich,
                             variables_both = variables_both,
                             variables_cont = variables_cont,
                             variables_dich = variables_dich,
                             opt_method = opt_method,
                             optimizer = optimizer,
                             lower = lower,
                             upper = upper,
                             formula_type_classic = formula_type_classic))



    fit <- flexmix::flexmix(formula = formula, data = data, k = k, model = model, control = control)
    rm(counter, envir = the) # counter will be created during the M-step driver

    return(fit)

  }else{
    # latent = "cont" or "dich"

    # Prepare first step
    if(is.null(id_col)){
      stop("id_col needed")
    }

    idframe <- data.frame(id = data[,id_col],type)
    idcount <- as.data.frame(table(unique(idframe)))
    data <- data[!is.element(as.character(data[,id_col]), as.character(idcount[idcount$Freq == 0,"id"])),]
    type <- idframe[!is.element(as.character(idframe[,"id"]), as.character(idcount[idcount$Freq == 0,"id"])),"type"]

    if(any(idcount$Freq == 0)){
      miss <- idcount[idcount$Freq == 0,"id"]
      warning(paste0(id_col ,paste(miss, collapse = ", "), " were removed, since they had only continuous or only dichotomous observations."))
    }

    # FIRST STEP: GET LATENT CLASSES
    if(latent == "cont"){
      data_cont <- data[type == type_cont,]
      assign("data", data_cont, envir = the)

      model <- list(FLXMRhyreg(type= type[type == type_cont],
                               stv = stv,
                               type_cont = type_cont,
                               type_dich = type_dich,
                               variables_both = variables_both,
                               variables_cont = variables_cont,
                               variables_dich = variables_dich,
                               opt_method = opt_method,
                               optimizer = optimizer,
                               lower = lower,
                               upper = upper,
                               formula_type_classic = formula_type_classic))


      mod <- flexmix::flexmix(formula = formula, data = data_cont, k = k, model = model, control = control)
      rm(counter, envir = the)



      # id_col must be something different than id if each observation should be classified
      data_cont$mod_comp <- mod@cluster

      data$roworder <- 1:nrow(data)
      data <- merge(data, unique(data_cont[,c(id_col,"mod_comp")]), by = id_col)
      data <- data[order(data$roworder), ]

      id_classes <- data_cont[,c(id_col,"mod_comp")]

    }

    if(latent == "dich"){
      data_dich <- data[type == type_dich,]
      assign("data", data_dich, envir = the)

      model <- list(FLXMRhyreg(type= type[type == type_dich],
                               stv = stv,
                               type_cont = type_cont,
                               type_dich = type_dich,
                               variables_both = variables_both,
                               variables_cont = variables_cont,
                               variables_dich = variables_dich,
                               opt_method = opt_method,
                               optimizer = optimizer,
                               lower = lower,
                               upper = upper,
                               formula_type_classic = formula_type_classic))


      mod <- flexmix::flexmix(formula = formula, data = data_dich, k = k, model = model, control = control)
      rm(counter, envir = the) # counter will be created during the M-step driver

      data_dich$mod_comp <- mod@cluster

      data$roworder <- 1:nrow(data)
      data <- merge(data, unique(data_dich[,c(id_col,"mod_comp")]), by = id_col)
      data <- data[order(data$roworder), ]

      # which id belongs to which class
      id_classes <- data_dich[,c(id_col,"mod_comp")]
    }

    # return only estimated classes without model new model coefficients
    if(classes_only == TRUE){
      return(id_classes)
    }


    # SECOND STEP: GET MODEL ESTIMATES

    # assign 1 to k in package environment to be used during M-step driver
    assign("k", 1, envir=the)

    data_list <- list()
    for(i in unique(mod@cluster)){
      data_list[[i]] <- data[data$mod_comp == i,]
    }


    mod_list <- lapply(data_list, function(xy){
      if(is.null(xy)){
        mod <- NULL
        warning( paste("One or more components are empty. Setting mod to NULL."))
      }else{
        model <- list(FLXMRhyreg(type= type[data$mod_comp == unique(xy$mod_comp)],
                                 stv = stv, # stv can be a list
                                 type_cont = type_cont,
                                 type_dich = type_dich,
                                 variables_both = variables_both,
                                 variables_cont = variables_cont,
                                 variables_dich = variables_dich,
                                 opt_method = opt_method,
                                 optimizer = optimizer,
                                 lower = lower,
                                 upper = upper,
                                 formula_type_classic = formula_type_classic))

        assign("data", xy, envir = the)

        mod <- flexmix::flexmix(formula = formula, data = xy, k = 1, model = model, control = control)
        rm(counter, envir = the) # counter will be created during the M-step driver
      }
      return(mod)
    })

    mod_list$id_classes <- unique(id_classes)
    return(mod_list)
  }
  rm(k, envir = the)
}





###############
### SUMMARY ###
###############


#' model summary for hyreg2 objects
#'
#' @description get model parameters of model generated by `hyreg2` or `hyreg2_het`
#'
#' @param object `modelobject` generated with [hyreg2()] or [hyreg2_het()]
#' @return `summary` object of [bbmle::mle2()] model,  Please note
#' that the outputs for `sigma` and `theta` are on a log-scale and have to be transformed using `exp()`to get the correct estimated values.
#'
#'
#' @author Svenja Elkenkamp
#' @examples
#'
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'k <- 1
#'stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control = list(iter.max = 1000, verbose = 4)

#'rm(counter)

#'mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both",
#'                     id_col = "id"
#')
#'summary_hyreg2(mod)
#'
#' @importFrom bbmle mle2
#' @export


summary_hyreg2 <- function(object){
  if(is.list(object)){
    object <- object[1:(length(object) - 1)]
    out <- lapply(object, function(k){
      comp <- k@components
      out_list <- lapply(comp, function(j){bbmle::summary(j[[1]]@parameters[["fit_mle"]])})
    })
  }else{
    comp <- object@components
    out <- lapply(comp, function(j){bbmle::summary(j[[1]]@parameters[["fit_mle"]])})
  }
  return(out)
}


#########################

#' extract parameter estimates as named vector
#'
#' @description function to export coefficient values and names from a `model` fitted by `hyreg2` or `hyreg2_het`
#' These values can be used as `stv` for a new model with `k > 1`
#' @param mod `model`output from [hyreg2] oder [hyreg2_het]. If `latent` was `"cont"` or `"dich"` only one element of the output list can be used.
#' @param comp `charachter`, default `"Comp.1"`. `"Comp.x"` indicating values from which model `component` (`x`) should be exported
#' @return `named vector` of parameter estimates from `mod`. Can be used as `stv` for additional model estimations using
#' `hyreg2` or `hyreg2_het`
#'
#' @author Svenja Elkenkamp
#' @examples
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'
#'k <- 1
#'stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control = list(iter.max = 1000, verbose = 4)

#'rm(counter)

#'mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both",
#'                     id_col = "id"
#')
#'new_stv <- get_stv(mod)
#'
#' # these new_stv can be used in an other estimation using hyreg2 as stv

#' @export

get_stv <- function(mod, comp = "Comp.1"){
  stv <- c(mod@components[[comp]][[1]]@parameters$coef,
           mod@components[[comp]][[1]]@parameters$sigma,
           mod@components[[comp]][[1]]@parameters$theta)
  return(stv)
}


##################################



