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
#' y<-impsdata$y
#' x<-impsdata$x
#' x_mis<-impsdata$x_mis
#' id<-impsdata$id
#' r<-impsdata$r
#' dist="binomial"
#' time<-4
#' candidate.sets<-list(c(1,2),c(1,4),c(1,2,4),c(1,2,4,7),c(1,2,3,4),c(1,2,3,4,5,6,7))
#' candidate.cor.sets<-c("exchangeable","ar1","independence")
#' #not run
#' #criterion.elcic<-ELCIC.wgee(x,y,x_mis,r,id,time,candidate.sets,
#' #name.var.sets=NULL,dist,candidate.cor.sets,joints=TRUE,lag=2)
#' #criterion.mlic<-MLIC.wgee(x,y,x_mis,r,id,time,candidate.sets,
#' #name.var.sets=NULL,dist,candidate.cor.sets,joints=TRUE,lag=2)
#' #criterion.qicw<-QICW.wgee(x,y,x_mis,r,id,time,candidate.sets,
#' #name.var.sets=NULL,dist,candidate.cor.sets,joints=TRUE,lag=2)

#'
"impsdata"
