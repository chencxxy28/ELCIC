#' The data are from a clinical trial of patients with respiratory illness, where 111 patients from two different clinics were randomized to receive either placebo or an active treatment. Patients were examined at baseline and at four visits during treatment.
#'
#' @title Data from a clinical trial comparing two treatments for a respiratory illness
#' @docType data
#'
#' @usage data(respiratorydata)
#'
#' @format An object of class \code{"list"}
#' \describe{
#'  \item{y}{Respiratory status at each visit, categorized as 1 = good, 0 = poor.}
#'  \item{x}{A full covariate matrix. It contains intercept, center(1: center 2, 0: center 1), sex (1: male,0: female), treat (1: treatment, 0: placebo), visit: id of each of four visits, baseline: respiratory status at baseline, age: in years at baseline.}
#'  \item{id}{Patient ID}
#'  \item{r}{A vector indicating the observation of data (1: observed; 0: missing). All data are observed.}
#' }
#' @references add here
#' @keywords datasets
#'

#'
"respiratorydata"

