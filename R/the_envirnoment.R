

### SETTING UP PACKAGE ###

#' creating environment for package internal objects
#' @export
#' @importFrom utils globalVariables


the <- new.env(parent = emptyenv())
the$test <- TRUE


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("counter"))
}
