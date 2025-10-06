
##############################
### Desciption of datasets ###
##############################


#' simulated_data_norm
#'
#' @format ## `simulated_data_norm`
#' A simulated data frame with 600 rows and 8 columns, following a combination of normal and binomial distribution
#' \describe{
#'   \item{type}{type of data}
#'   \item{y}{result of the model y = x1 + x2 + x3}
#'   \item{x1, x2, x3}{random numbers from rnorm}
#'   \item{class}{original class of the data point}
#'   \item{id}{id number of observations to simulated different persons}
#'   \item{y_cens}{column y censored at 3}
#'   ...
#' }
#' @source {simulated with true parameter values:
#'  Class 1: sigma  = 1.0, theta = 5 and c(x1,x2,x3) = c(0.5, -0.3, 0.8)
#'  Class 2: sigma  = 0.5, theta = 2 and c(x1,x2,x3) = c(1.4, 2.3, -0.2)}

"simulated_data_norm"


#' simulated_data_mo
#'
#' @format ## `simulated_data_mo`
#' A simulated data frame with 480 rows and 9 columns, following a combination of normal and binomial distribution
#' \describe{
#'   \item{type}{type of data}
#'   \item{y}{result of the model y = -1 + mo2 + mo3 + mo4 + mo5}
#'   \item{mo2, mo3, mo4, mo5}{dummy variables}
#'   \item{class}{original class of the data point}
#'   \item{id}{id number of observations to simulated different persons}
#'   \item{y_cens}{column y censored at 0 (lower boundary)}
#'   ...
#' }
#' @source {simulated with true parameter values:
#'  Class 1: sigma  = XXX, theta = XXX and c(mo2,mo3,mo4,mo5) = c(XXX)
#'  Class 2: sigma  = XXX, theta = XXX and c(mo2,mo3,mo4,mo5) = c(XXX)}

"simulated_data_mo"





#'  simulated_data
#'
#' @format ## `simulated_data`
#' A simulated data frame with 480 rows and 25 columns, following a combination of normal and binomial distribution
#' \describe{
#'   \item{type}{type of data}
#'   \item{y}{result of the model y = -1 + mo2 + mo3 + ... + ad4 + ad5}
#'   \item{mo2, mo3, mo4, mo5,
#'          sc2, sc3, sc4, sc5,
#'          ua2, ua3, ua4, ua5,
#'          pd2, pd3, pd4, pd5,,
#'          ad2, ad3, ad4, ad5,}{dummy variables for EQ5D data simulation}
#'   \item{class}{original class of the data point}
#'   \item{id}{id number of observations to simulated different persons}
#'   \item{y_cens}{column y censored at 2 (upper boundary)}
#'   ...
#' }
#' @source {simulated with true parameter values:
#'  Class 1: sigma  = XXX, theta = XXX and c(mo2,mo3,mo4,mo5, XXX) = c(XXX)
#'  Class 2: sigma  = XXX, theta = XXX and c(mo2,mo3,mo4,mo5, XXX) = c(XXX)}
#'
 "simulated_data"



