
##############################
### Desciption of datasets ###
##############################


#' simulated_data_norm
#'
#' @format ## `simulated_data_norm`
#' A simulated data frame with 600 rows and 9 columns, following a combination of normal and binomial distribution
#' \describe{
#'   \item{type}{type of data, "TTO" indicates normal distribution, "DCE_A" indicates binomial distribution}
#'   \item{y}{result of the formula y ~ x1 + x2 + x3}
#'   \item{x1, x2, x3}{random numbers from rnorm}
#'   \item{class}{original class of the data point}
#'   \item{id}{id number of observations to simulated different persons}
#'   \item{y_non}{result of the formula y ~ (x1 \* beta1 + x2 \* beta3) \* (x1 \* beta1 + x3 \* beta3)}
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
#'   \item{type}{type of data, "TTO" indicates normal distribution, "DCE_A" indicates binomial distribution}
#'   \item{y}{result of the formula y ~ -1 + mo2 + mo3 + mo4 + mo5}
#'   \item{mo2, mo3, mo4, mo5}{dummy variables}
#'   \item{class}{original class of the data point}
#'   \item{id}{id number of observations to simulated different persons}
#'   \item{y_cens}{column y censored at 0 (lower boundary)}
#'   ...
#' }
#' @source {simulated with true parameter values:
#'  Class 1: sigma  = 0.001, theta = 0.2 and c(mo2,mo3,mo4,mo5) = c(0.005, 0.01, 0.08, 0.1)
#'  Class 2: sigma  = 0.1, theta = 2 and c(mo2,mo3,mo4,mo5) = c(0.2, 0.4, 0.6, 0.8)}

"simulated_data_mo"





#'  simulated_data
#'
#' @format ## `simulated_data`
#' A simulated data frame with 480 rows and 25 columns, following a combination of normal and binomial distribution
#' \describe{
#'   \item{type}{type of data, "TTO" indicates normal distribution, "DCE_A" indicates binomial distribution}
#'   \item{y}{result of the formula y ~ -1 + mo2 + mo3 + ... + ad4 + ad5}
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
#'  Class 1: sigma  = 0.02, theta = 2 and
#'   c(mo2,mo3,mo4,mo5) = c(0.001, 0.05, 0.08, 0.1),
#'   c(sc2,sc3,sc4,sc5) = c(0.01, 0.2, 0.36, 0.5),
#'   c(ua2,au3,ua4,ua5) = c(0.015, 0.25, 0.5, 0.8),
#'   c(pd2,pd3,pd4,pd5) = c(0.1, 0.3, 0.4, 0.6),
#'   c(ad2,ad3,ad4,ad5) = c(0.09, 0.19, 0.6, 0.7)
#'
#'  Class 2: sigma  = 0.1, theta = 3 and
#'   c(mo2,mo3,mo4,mo5) = c(0.2, 0.4, 0.6, 0.8),
#'   c(sc2,sc3,sc4,sc5) = c(0.1, 0.3, 0.4, 0.5),
#'   c(ua2,au3,ua4,ua5) = c(0.2, 0.25, 0.6, 0.7),
#'   c(pd2,pd3,pd4,pd5) = c(0.05, 0.2, 0.27, 0.8),
#'   c(ad2,ad3,ad4,ad5) = c(0.15, 0.35, 0.4, 0.65)
#' }
#'
 "simulated_data"



