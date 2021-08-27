#' Toydata generated for variable selection under GEE framework without missingness or missing completely at random
#'
#'
#' @docType data
#'
#' @usage data(wgeesimdata)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{y}{The outcomes generated from Biomial distribution with three repeated measurements from each subject. The data is dropout missing.}
#'  \item{x}{A full covariate matrix. The first column corresponds to the intercept; the second column contains time-independent variable x1; the third column contains a doctor-visit variable x2; the third column contains time-dependent variable x3.}
#'  \item{x_mis}{A covariate matrix for missing data model. The first column corresponds to the intercept; the second column contains continuous variable x_mis1; the third column contains outcome y with lag-1 (x_mis2).}
#'  \item{id}{A vector indicating subject id.}
#'  \item{obs_ind}{A vector indicating missingness: 1 for observed records, and 0 for unobserved records.}
#' }
#' @references This data set was artificially created for the ELCIC package.
#' @keywords datasets
#' @examples
#'
#' data(wgeesimdata)
#' wgeesimdata$y
#'
"wgeesimdata"

