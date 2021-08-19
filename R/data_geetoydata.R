#' Toydata generated for variable selection under GEE framework without missingness or missing completely at random
#'
#'
#' @docType data
#'
#' @usage data(geetoydata)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{y}{The outcomes generated from Poisson distribution with three repeated measurements from each subject}
#'  \item{x}{A covariate matrix, of which the first column are all ones and rest columns contain normally distributed. Two are time-dependent variables, and one is time-independent variable.}
#' }
#' @references This data set was artificially created for the ELCIC package.
#' @keywords datasets
#' @examples
#'
#' data(geetoydata)
#' geetoydata$y
#'
"geetoydata"
