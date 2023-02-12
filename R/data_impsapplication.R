#' The imps data frame has 1544 rows and 8 columns. The data is from National Institute of the Mental Health Schizophrenia Collaborative Study, where the effect of chlorpromazine, fluphenazine, or thioridazine treatment on the overall severity of the schizophrenia disorder is of interest.
#'
#' @title Inpatient Multidimensional Psychiatric Scale (IMPS)
#' @docType data
#'
#' @usage data(impsdata)
#'
#' @format An object of class \code{"list"}
#' \describe{
#'  \item{y}{The binary outcomes indicating whether IMPS >=4, which is longitudinal dropout and missing at random. IMPS describes severity of the schizophrenia disorder (ranges from 0 to 7)}
#'  \item{x}{A full covariate matrix. It contains intercept, sex (1:male,0:female), drug (1: chlorphromazine, fluphenazine, or thioridazine treatment; 0: placebo), time: square root of the week covariate, and their two-way interactions.}
#'  \item{x_mis}{A covariate matrix for missing data model. It contains intercept, drug, time, and sex.}
#'  \item{id}{Patient ID}
#'  \item{r}{An indicator of the missingness (1: observed; 0: missing).}
#' }
#' @references add here
#' @keywords datasets
#' @examples
#'
#' data(impsdata)
#' id<-impsdata$id
#' r<-impsdata$r
#' data.real <- data.frame(id=id,y=impsdata$y, impsdata$x)
#' head(data.real,n=10)
#' # each participant has three records
#' time<-4
#' # the outcome is binary
#' family=binomial()
#' models <- list(y~Time, y~Drug,y~Time+Drug,
#'                y~Time*Drug,y~Time+Sex+Drug,
#'                y~Time+Sex+Drug+Time:Sex+Sex:Drug+Drug:Time)
#'
#' model_mis<-r~Drug+Time+Sex
#'
#' candidate.cor.sets<-c("exchangeable","independence","ar1")
#' #not run
#' #output_ELCIC<-ELCICwgee(models, candidate.cor.sets,data=data.real,model_mis,family,r,id,time)
#' #output_MLIC<-MLICwgee(models, candidate.cor.sets,data=data.real,model_mis,family,r,id,time)
#' #output_QICW<-QICWwgee(models, candidate.cor.sets,data=data.real,model_mis,family,r,id,time)

#'
"impsdata"

