


# Helper functions not for user use #



######################
### refit function ###
######################

# internal function, not for direct user use
gendf <- function(list){
  df <- list()
  for(j in 1:length(list)){
    df[[j]] <- data.frame("Estimates" = c(list[[j]]@parameters[["coef"]], #  maybe include an other apply function here?
                                          list[[j]]@parameters[["sigma"]],
                                          list[[j]]@parameters[["theta"]]),
                          "Std Error" = bbmle::summary(list[[j]]@parameters[["fit_mle"]])@coef[,2],
                          "pvalue" =bbmle::summary(list[[j]]@parameters[["fit_mle"]])@coef[,4])
    # give AIC, loglik here as additional arguments?
  }

  return(df)
}


### refit function ###

# get parametervalues of models
# object is outcome of hyreg2
refit <- function(object){
  if(is.list(object)){
    out <- lapply(object, function(k){
      comp <- k@components
      out_list <- lapply(comp, function(j){gendf(j)})
    })
  }else{
    comp <- object@components
    out <- lapply(comp, function(j){gendf(j)})
  }
  return(out)
}






### HELPER FUNCTIONS FOR USE OF NON-CLASSIC FORMULAS ####

################################
### eval non-linear formulas ###
################################

# function to get Xb during M-step for non-classic case
eval_formula_non <- function(formula, data, stv, form_back = F){

  rhs <- formula[[3]]
  int <- is.element("INTERCEPT",all.vars(formula[[3]])) # intercept check

  # if(int & !is.element("(Intercept)", colnames(data))){
  #   data$"(Intercept)" <- rep(1,dim(data)[1])
  # }

  form_wght <- replace_vars(rhs, stv)

  # new env to be used in eval()
  env <- list2env(as.list(data), parent = parent.frame())

  # check for intercept
  if(int){
    intercept_val <- stv["INTERCEPT"]
    if(is.na(intercept_val)){
      intercept_val <- 0
    }
    form_wght <- call("+", intercept_val, form_wght)
  }

  if(form_back == TRUE){ # return formula before it is evaluated
    return(form_wght)
  }

  # get Xb, eval expression using stv (used in form_wght)
  xb <- eval(form_wght, envir = env)

  return(xb)
}




# replace variables by numbers given in stv
# fct used in eval_formula_non
replace_vars <- function(formula, stv){

  if(is.symbol(formula)){

    name <- as.character(formula)

    if(is.element(name, names(stv))){
      return(stv[[name]])
    }else{
      return(formula)
    }

  }else if(is.call(formula)){
    return(as.call(lapply(formula, replace_vars, stv = stv)))
  }else{
    return(formula)
  }
}




#####################################
### get formula for flexmix input ###
#####################################

# which variables from formula are part of data colnames?
# create formula to be uses in flexmix call
# returns a formula which can be used in flexmix

get_data_vars <- function(formula, data) {

  # split up formula
  #formula_parts <- strsplit(deparse(formula), "\\|")[[1]]

  # allow long formulas
  formula_str <- paste(deparse(formula, width.cutoff = 500), collapse = " ")
  formula_parts <- strsplit(formula_str, "\\|")[[1]]

  # find conditional part
  if(length(formula_parts) > 1){
    conditional <- trimws(formula_parts[2])
  }else{
    conditional <- NULL
  }

  # main formula
  main_formula <- as.formula(trimws(formula_parts[1]))

  # check which vars in formula are in data
  vars_data <- intersect(all.vars(main_formula[[3]]), names(data))

  rhs <- paste(vars_data, collapse = " + ")

  # intercept
  int <- is.element("INTERCEPT",all.vars(main_formula[[3]]))
  if(!int){
    rhs <- paste("-1 +", rhs)
  }

  # new main formula
  new_formula_str <- paste(deparse(main_formula[[2]]), "~", rhs)

  # add conditional part
  if (!is.null(conditional)) {
    new_formula_str <- paste(new_formula_str, "|", conditional)
  }


  return(as.formula(new_formula_str))
}




###################################################
### data and formula handling chehcks in hyreg2 ###
###################################################

# simplify formula for only cont or dich variables
simplify_formula <- function(formula){

  if(is.call(formula)){
    fun <- as.character(formula[[1]])
    args <- as.list(formula[-1])
    args <- lapply(args, simplify_formula)

    # Multiply
    if(fun == "*"){
      #if one argument is 0, everything is 0
      if (any(vapply(args, function(a) is.numeric(a) && a == 0, logical(1)))) {
        return(0)
      }

      # 1 * x → x
      args <- Filter(function(a) !(is.numeric(a) && a == 1), args)
      if(length(args) == 0){
        return(1)
      }
      if(length(args) == 1){
        return(args[[1]])
      }
      return(as.call(c(list(as.name("*")), args)))
    }

    # addition
    if (fun == "+") {

      # 0 + x → x
      args <- Filter(function(a) !(is.numeric(a) && a == 0), args)
      if(length(args) == 0){
        return(0)
      }
      if(length(args) == 1){
        return(args[[1]])
      }
      return(as.call(c(list(as.name("+")), args)))
    }

    # return everything else without canges
    return(as.call(c(list(as.name(fun)), args)))
  }else{
    return(formula)
  }
}


# check which vars from formula are undefined here (not part of variables for one type of data/ stv)
formulacheck_variables <- function(formula, data, stv) {
  if (is.symbol(formula)) {
    name <- as.character(formula)
    if (name %in% c(colnames(data), names(stv))) {
      return(formula)
    } else {
      return(0)
    }
  } else if (is.call(formula)) {
    fun <- formula[[1]]
    args <- as.list(formula[-1])
    new_args <- lapply(args, formulacheck_variables, data = data, stv = stv)
    out <- as.call(c(list(fun), new_args))
    return(simplify_formula(out))
  } else {
    return(formula)
  }
}

# check formula for variables in stv,
# fct used in stv check for non-classic fcts
get_stv_vars <- function(formula, stv) {
  return(intersect(all.vars(formula[[3]]), names(stv)))
}
