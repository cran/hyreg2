## ----chunk-1, message=FALSE, warning=FALSE------------------------------------
# simulate dataframe called simulated_data_norm

library(hyreg2)
library(flexmix)


### CLASS 1 ###

# Set parameters
set.seed(42)
n_samples_tto1 <- 200
n_samples_dce1 <- 100
sigma1 <- 1.0
theta1 <- 5
beta1 <- c(0.5, -0.3, 0.8)  # true values for beta

# simulate matrix
x_tto1 <- matrix(rnorm(n_samples_tto1 * length(beta1)), ncol = length(beta1))
x_dce1 <- matrix(rnorm(n_samples_dce1 * length(beta1)), ncol = length(beta1))



#Xb (linear predictor following formula xb = x1*beta1 + x2*beta2 + x3*beta3)
Xb_tto1 <- x_tto1 %*% beta1
Xb_dce1 <- (x_dce1 %*% beta1) * theta1


# generate y
y_tto1 <- rnorm(n_samples_tto1, mean = Xb_tto1, sd = sigma1)  # continuous outcomes
logistic_tmp1 <- 0.5 + 0.5 * tanh(Xb_dce1 / 2)
y_dce1 <- rbinom(n_samples_dce1, size = 1, prob = logistic_tmp1)  # binary outcomes



### dataset ###
data_tto1 <- data.frame(
  type = rep("TTO", n_samples_tto1),
  y = y_tto1,
  x1 = x_tto1[, 1],
  x2 = x_tto1[, 2],
  x3 = x_tto1[, 3]
)

data_dce1 <- data.frame(
  type = rep("DCE_A", n_samples_dce1),
  y = y_dce1,
  x1 = x_dce1[, 1],
  x2 = x_dce1[, 2],
  x3 = x_dce1[, 3]
)

simulated_data1 <- rbind(data_tto1, data_dce1)


### CLASS 2 ###

# Set parameters
set.seed(37)
n_samples_tto2 <- 200
n_samples_dce2 <- 100
sigma2 <- 0.5
theta2 <- 2
beta2 <- c(1.4, 2.3, -0.2)   # true values for beta


# simulate matrix
x_tto2 <- matrix(rnorm(n_samples_tto2 * length(beta2)), ncol = length(beta2))
x_dce2 <- matrix(rnorm(n_samples_dce2 * length(beta2)), ncol = length(beta2))


#  Xb (linear predictor following formula xb = x1*beta1 + x2*beta2 + x3*beta3)
Xb_tto2 <- x_tto2 %*% beta2
Xb_dce2 <- (x_dce2 %*% beta2) * theta2


# generate y
y_tto2 <- rnorm(n_samples_tto2, mean = Xb_tto2, sd = sigma2)  # continuous outcomes
logistic_tmp2 <- 0.5 + 0.5 * tanh(Xb_dce2 / 2)
y_dce2 <- rbinom(n_samples_dce2, size = 1, prob = logistic_tmp2)  # binary outcomes



### dataset ###
data_tto2 <- data.frame(
  type = rep("TTO", n_samples_tto2),
  y = y_tto2,
  x1 = x_tto2[, 1],
  x2 = x_tto2[, 2],
  x3 = x_tto2[, 3]
)

data_dce2 <- data.frame(
  type = rep("DCE_A", n_samples_dce2),
  y = y_dce2,
  x1 = x_dce2[, 1],
  x2 = x_dce2[, 2],
  x3 = x_dce2[, 3]
)

simulated_data2 <- rbind(data_tto2, data_dce2)


### DATASET with two classes ###
simulated_data1$class <- 1
simulated_data2$class <- 2

# set id numerbs
simulated_data1$id <- c(1:100, 1:100, 1:100)
simulated_data2$id <- c(101:200, 101:200, 101:200)


simulated_data_norm<- rbind(simulated_data1,simulated_data2)


# clean environment
rm(list = setdiff(ls(), "simulated_data_norm"))


## ----chunk-2, message=FALSE, warning=FALSE------------------------------------

# visualize true classes from simulated_data_norm
p <- plot_hyreg2(data = simulated_data_norm,
           x = "id",
           y = "y",
           id_col = "id",
           class_df_model = simulated_data_norm[,c("id","class")],
           colors = c("#F8766D","#00BFC4"))
p <- p +  ggplot2::labs(title = "Figure 1")
p


## ----chunk-3, message=FALSE, warning=FALSE, results='hide'--------------------

set.seed(277)
formula <- y ~  -1 + x1 + x2 + x3
k <- 2
stv <- setNames(c(0.2,0,1,1,1),c("x1","x2","x3","sigma","theta"))
control <- list(iter.max = 1000, verbose = 4)

mod <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "both"
)


## ----chunk-4, message=FALSE, warning=FALSE------------------------------------

summary_hyreg2(mod)


## ----chunk-5, message=FALSE, warning=FALSE------------------------------------

# if model was estimated without grouping variable
simulated_data_norm$observation <- rownames(simulated_data_norm)


p <- plot_hyreg2(data = simulated_data_norm,
           x = "id",
           y = "y",
           id_col = "observation",
           class_df_model = give_class(data = simulated_data_norm,
                                 model = mod,
                                 id_col = "observation"))
p <- p +  ggplot2::labs(title = "Figure 2")
p


## ----chunk-6, message=FALSE, warning=FALSE------------------------------------

# only continuous datapoints 
p <- plot_hyreg2(data = simulated_data_norm,
           x = "id",
           y = "y",
           id_col = "observation",
           class_df_model = give_class(data = simulated_data_norm,
                                 model = mod,
                                 id_col = "observation"),
           type_to_plot = list("type","TTO"))
p <- p +  ggplot2::labs(title = "Figure 3")
p



## ----chunk-7, message=FALSE, warning=FALSE------------------------------------

# only dichotomous datapoints
p <- plot_hyreg2(data = simulated_data_norm,
           x = "id",
           y = "y",
           id_col = "observation",
           class_df_model = give_class(data = simulated_data_norm,
                                 model = mod,
                                 id_col = "observation"),
           type_to_plot = list("type","DCE_A"))
p <- p +  ggplot2::labs(title = "Figure 4")
p



## ----chunk-8, message=FALSE, warning=FALSE------------------------------------

# ratio of correct classified datapoints
(sum(mod@cluster == simulated_data_norm$class))/dim(simulated_data_norm)[1]



# cross table of observation types and assigned classes
proof <- merge(simulated_data_norm,
                 give_class(simulated_data_norm,mod,"observation"), by = "observation")
table(proof$type,proof$mod_comp)





## ----chunk-9, message=FALSE, warning=FALSE, results='hide'--------------------
set.seed(227)
formula <- y ~  -1 + x1 + x2 + x3 | id
k <- 2
stv <- setNames(c(0.2,0,1,1,1),c("x1","x2","x3","sigma","theta"))
control <- list(iter.max = 1000, verbose = 4)

mod_id <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "both"
)


## ----chunk-9a, message=FALSE, warning=FALSE-----------------------------------

summary_hyreg2(mod_id)



## ----chunk-10, message=FALSE, warning=FALSE-----------------------------------

p <- plot_hyreg2(data = simulated_data_norm,
           x ="id",
           y = "y",
           id_col ="id",
           class_df_model = unique(give_class(data = simulated_data_norm,
                                        model = mod_id,
                                        id_col = "id")))
p <- p +  ggplot2::labs(title = "Figure 5")
p



## ----chunk-13, message=FALSE, warning=FALSE-----------------------------------

# ratio of correct classified datapoints
(sum(mod_id@cluster == simulated_data_norm$class))/dim(simulated_data_norm)[1]


# cross table of observation types and assigned classes

proof <- merge(simulated_data_norm,
                 unique(give_class(simulated_data_norm,mod_id,"id")), by = "id")
table(proof$type,proof$mod_comp)




## ----chunk-14, message=FALSE, warning=FALSE, results='hide'-------------------
set.seed(227)
formula <- y ~  -1 + x1 + x2 + x3 | id
k <- 2
stv <- setNames(c(0.2,0,1,1,1),c("x1","x2","x3","sigma","theta"))
control <- list(iter.max = 1000, verbose = 4)

mod_cont <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "cont",
                     id_col = "id"
)

## ----chunk-14a, message=FALSE, warning=FALSE----------------------------------

summary_hyreg2(mod_cont)


## ----chunk-15, message=FALSE, warning=FALSE-----------------------------------

p <- plot_hyreg2(data = simulated_data_norm,
           x ="id",
           y = "y",
           id_col ="id",
           class_df_model = give_class(data = simulated_data_norm,
                                model = mod_cont,
                                id_col = "id"))
p <- p +  ggplot2::labs(title = "Figure 6")
p



## ----chunk-16, message=FALSE, warning=FALSE-----------------------------------

# ratio of correct classified datapoints
proof <- merge(unique(simulated_data_norm[,c("id","class")]),mod_cont[["id_classes"]], by = "id")
sum((proof$class == proof$mod_comp)/dim(proof)[1])


# cross table of observation types and assigned classes

proof <- merge(simulated_data_norm,
                mod_cont[[k+1]], by = "id")
table(proof$type,proof$mod_comp)




## ----chunk-18, message=FALSE, warning=FALSE, results='hide'-------------------

# upper bound for simulated_data_norm$y at 3
# new variable y_cens

for(i in 1:length(simulated_data_norm$y)){
  if(simulated_data_norm$y[i] <= 3){
    simulated_data_norm$y_cens[i] <- simulated_data_norm$y[i]
  }else{
    simulated_data_norm$y_cens[i] <- 3
  }
}


set.seed(227)
formula <- y_cens ~  -1 + x1 + x2 + x3 | id
k <- 2
stv <- setNames(c(0.2,0,1,1,1),c("x1","x2","x3","sigma","theta"))
control <- list(iter.max = 1000, verbose = 4)

mod_cens <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     upper = 3,
                     latent = "cont",
                     id_col = "id"
)


## ----chunk-19, message=FALSE, warning=FALSE, results='hide', fig.show='hide'----
summary_hyreg2(mod_cens)


plot_hyreg2(data = simulated_data_norm,
           x ="id",
           y = "y_cens",
           id_col = "id",
           class_df_model = give_class(data = simulated_data_norm,
                                 model = mod_cens,
                                 id_col = "id"))




## ----chunk-22, message=FALSE, warning=FALSE, results='hide'-------------------
set.seed(227)
formula <- y ~  -1 + x1 + x2 + x3 | id
k <- 2
stv1 <- setNames(c(0.2,0,1,1,1),c("x1","x2","x3","sigma","theta"))
stv2 <- setNames(c(1,2,0.5,1,1),c("x1","x2","x3","sigma","theta"))
control <- list(iter.max = 1000, verbose = 4)

mod_stvl<- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = list(stv1,stv2),
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "cont",
                     id_col = "id"
)

summary_hyreg2(mod_stvl)


## ----chunk-24, message=FALSE, warning=FALSE, results='hide'-------------------

set.seed(227)
formula <- y ~  -1 + x1 + x2 + x3 | id
formula_sigma <- y ~ x2 + x3
k <- 2
stv1 <- get_stv(mod_cont[[1]])[-4]  # -4 since "sigma" is on fourth place in return vector
stv2 <- get_stv(mod_cont[[2]])[-4]  # -4 since "sigma" is on fourth place in return vector
stv_sigma <- setNames(rep(0.1,3), c("x2","x3","(Intercept)"))
control <- list(iter.max = 1000, verbose = 4)


mod_het <- hyreg2_het(formula = formula,
                     formula_sigma = formula_sigma,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = list(stv1,stv2),
                     stv_sigma = stv_sigma,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "cont",
                     id_col = "id"
)

## ----chunk-24a, message=FALSE, warning=FALSE----------------------------------
summary_hyreg2(mod_het)


## ----chunk-23, message=FALSE, warning=FALSE, results='hide'-------------------
set.seed(227)
formula <- y ~  -1 + x1 + x2 + x3 | id
k <- 2
stv <- setNames(c(0.2,0,1,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
control <- list(iter.max = 1000, verbose = 4)

mod_partialc <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "cont",
                     id_col = "id",
                     variables_both = "x1",
                     variables_cont = "x3",
                     variables_dich = "x2"
)

summary_hyreg2(mod_partialc)


## ----chunk-26, message=FALSE, warning=FALSE-----------------------------------

# simulate vector y_non and add this to simulated_data_norm

#library(hyreg2)

### CLASS 1 ###

# Set parameters (same as in simulation above)
set.seed(42)
n_samples_tto1 <- 200
n_samples_dce1 <- 100
sigma1 <- 1.0
theta1 <- 5
beta_non1 <- c(0.5, -0.3, 0.8)  # true values for beta

# simulate matrix
x_tto1 <- matrix(rnorm(n_samples_tto1 * length(beta_non1)), ncol = length(beta_non1))
x_dce1 <- matrix(rnorm(n_samples_dce1 * length(beta_non1)), ncol = length(beta_non1))


### NON-classic ###
# x1,x2 and x3 are columnnames of the dataset, beta1, beta2 and beta3 are the parameters to be estimated
formula1 <- y ~ (x1 * beta1 + x2 * beta2 ) * (x1 *beta1 + x3 * beta3)
x_tto_non1 <- data.frame(x1 = x_tto1[,1], x2 = x_tto1[,2], x3 = x_tto1[,3])
x_dce_non1 <- data.frame(x1 = x_dce1[,1], x2 = x_dce1[,2], x3 = x_dce1[,3])
beta_non1 <- setNames(beta_non1[1:3], c("beta1","beta2","beta3"))

Xb_tto_non1 <- hyreg2:::eval_formula_non(formula1, x_tto_non1, beta_non1)
Xb_dce_non1 <- (hyreg2:::eval_formula_non(formula1, x_dce_non1, beta_non1)) * theta1


# generate y_non
y_tto_non1 <- rnorm(n_samples_tto1, mean = Xb_tto_non1, sd = sigma1)  # continuous outcomes
logistic_tmp_non1 <- 0.5 + 0.5 * tanh(Xb_dce_non1 / 2)
y_dce_non1 <- rbinom(n_samples_dce1, size = 1, prob = logistic_tmp_non1)  # binary outcomes



### CLASS 2 ###

# Set parameters (same as in simulation above)
set.seed(37)
n_samples_tto2 <- 200
n_samples_dce2 <- 100
sigma2 <- 0.5
theta2 <- 2
beta_non2 <- c(1.4, 2.3, -0.2)   # true values for beta


# simulate matrix
x_tto2 <- matrix(rnorm(n_samples_tto2 * length(beta_non2)), ncol = length(beta_non2))
x_dce2 <- matrix(rnorm(n_samples_dce2 * length(beta_non2)), ncol = length(beta_non2))



### NON-classic ###
# x1,x2 and x3 are columnnames of the dataset, beta1, beta2 and beta3 are the parameters to be estimated
formula2 <- y  ~ (x1 * beta1 + x2 * beta2 ) * (x1 * beta1 + x3 * beta3)
x_tto_non2 <- data.frame(x1 = x_tto2[,1], x2 = x_tto2[,2], x3 = x_tto2[,3])
x_dce_non2 <- data.frame(x1 = x_dce2[,1], x2 = x_dce2[,2], x3 = x_dce2[,3])
beta_non2 <- setNames(beta_non2[1:3], c("beta1","beta2","beta3"))

Xb_tto_non2 <- hyreg2:::eval_formula_non(formula2, x_tto_non2, beta_non2)
Xb_dce_non2 <- (hyreg2:::eval_formula_non(formula2, x_dce_non2, beta_non2)) * theta2


# generate y_non
y_tto_non2 <- rnorm(n_samples_tto2, mean = Xb_tto_non2, sd = sigma2)  # continuous outcomes
logistic_tmp_non2 <- 0.5 + 0.5 * tanh(Xb_dce_non2 / 2)
y_dce_non2 <- rbinom(n_samples_dce2, size = 1, prob = logistic_tmp_non2)  # binary outcomes


simulated_data_norm$y_non <- c(y_tto_non1,y_dce_non1,y_tto_non2,y_dce_non2)

# clean environment
rm(n_samples_tto1,
   n_samples_dce1,
   sigma1,
   theta1,
   beta_non1,
   n_samples_tto2,
   n_samples_dce2,
   sigma2,
   theta2,
   beta_non2,
   logistic_tmp1,
   logistic_tmp2,
   logistic_tmp_non1,
   logistic_tmp_non2,
   x_tto_non1,
   x_tto_non2,
   x_dce_non1,
   x_dce_non2,
   Xb_tto_non1,
   Xb_tto_non2,
   Xb_dce_non1,
   Xb_dce_non2,
   x_dce1,
   x_dce2,
   x_tto1,
   x_tto2,
   y_dce_non1,
   y_dce_non2,
   y_tto_non1,
   y_tto_non2,
   formula1,
   formula2
)


## ----chunk-25, message=FALSE, warning=FALSE, results='hide'-------------------
set.seed(259)
formula <- y_non ~ (x1 * beta1 + x2 * beta2 ) * (x1 * beta1 + x3 * beta3)| id
k <- 2
stv <- setNames(c(0.2,0,1,1,1),c("beta1","beta2","beta3","sigma","theta"))
control <- list(iter.max = 1000, verbose = 4)

mod_non <- hyreg2(formula = formula,
                     data =  simulated_data_norm,
                     type =  simulated_data_norm$type,
                     stv = stv,
                     k = k,
                     type_cont = "TTO",
                     type_dich = "DCE_A",
                     opt_method = "L-BFGS-B",
                     control = control,
                     latent = "cont",
                     id_col = "id",
                     formula_type_classic = FALSE
)


## ----chunk-21, message=FALSE, warning=FALSE, results='hide'-------------------

summary_hyreg2(mod_non)

p <- plot_hyreg2(data = simulated_data_norm, 
           x ="id",
           y = "y_non",
           id_col = "id",
           class_df_model = give_class(data = simulated_data_norm,
                                 model = mod_non,
                                 id_col = "id"))
p <- p +  ggplot2::labs(title = "Figure 7")
p



