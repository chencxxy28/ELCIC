

#'@title Variable selection based on ELCIC under the syntax of GLM (Main function).
#'@description The function \code{\link{ELCICglm}} provides the variable selection under the syntax of GLM package.
#'@usage ELCICglm(models, data, family)
#'@param models A list of formulas. See corresponding documentation to glm.
#'@param data A data frame containing the variables in the model.
#'@param family A description of the error distribution and link function to be used in the model.
#'The details are given under "Details".
#'@return A list with two items: model selection result based on ELCIC, AIC, BIC, and GIC;
#'An object of "glm" based on the selected model.
#'@details Three common distributions are considered: "gaussian", "poisson", "binomial".
#'For the current package, the identical link is consider for "gaussian" distribution;
#'the log link is consider for "poisson" distribution; the logit link is consider for "binomial" distribution;
#'
#'
#'
#'@examples
#'## tests
#'# load data
#'data(glmsimdata)
#'dat <- data.frame(y=glmsimdata$y, glmsimdata$x)
#'models <- list(y~x1, y~x1+x2, y~x1+x2+x3)
#'output<-ELCICglm(models, dat, poisson())
#'output$model.selection
#'output$glm.output
#'
#'@export
#'@importFrom stats poisson binomial gaussian model.response model.frame model.matrix


ELCICglm <- function(models, data, family)
{
    if (!inherits(family, "family"))
        stop("dist must be an object of type 'family'")
    ## dist should be passed directly to glm, but I don't here
    ## because I don't want to change the package functions
    dist <- family$family
    if (!is.list(models))
        models <- list(models)
    y <- model.response(model.frame(models[[1]], data))
    l <- lapply(models, function(li) m <- model.matrix(li, data))
    name.var.sets <- lapply(l, colnames)
    dat <- do.call(cbind, l)
    dat <- dat[,!duplicated(colnames(dat))]
    ci <- ELCIC.glm(x=dat, y=y, name.var.sets = name.var.sets, dist=dist)
    attr(ci, "formula") <- models
    attr(ci, "type") <- "glm"
    class(ci) <- "elcic"

    print.elcic.glm(ci)

    ##fit glm based on ELCIC
    fit.glm<-glm(models[[which.min(ci[1,])]], data, family=family)

    list(model.selection=ci,glm.output=fit.glm)
}




#'@title Model selection based on ELCIC under the syntax of GEE (Main function).
#'@description The function \code{\link{ELCICgee}} provides the model selection under the syntax of geepack package.
#'@usage ELCICgee(models, candidate.cor.sets,data, family,r,id,time)
#'@param models A list of formulas. See corresponding documentation to geeglm.
#'@param candidate.cor.sets A vector containing candidate correlation structures. It can be any subset of c("independence","exchangeable", "ar1").
#'@param data A data frame containing the variables in the model.
#'@param family A description of the error distribution and link function to be used in the model.
#'The details are given under "Details".
#'@param r A vector indicating the observation of data: 1 for observed records (both outcome and covariates are observed for a given subject), and 0 for unobserved records.
#'The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject.
#'@return A list with two items: model selection result based on ELCIC;
#'An object of "geeglm" based on the selected model.
#'@details Three common distributions are considered: "gaussian", "poisson", "binomial".
#'For the current package, the identical link is consider for "gaussian" distribution;
#'the log link is consider for "poisson" distribution; the logit link is consider for "binomial" distribution;
#'
#'
#'
#'@examples
#'## tests
#'# load data
#'data(geesimdata)
#'id<-geesimdata$id
#'r<-rep(1,length(id))
#'time<-3
#'dat <- data.frame(y=geesimdata$y, geesimdata$x,id=id)
#'models <- list(y~x1, y~x1+x2)
#'candidate.cor.sets<-c("exchangeable")
#'family<-poisson()
#'output<-ELCICgee(models, candidate.cor.sets,data=dat,family,r,id,time)
#'output$model.selection
#'output$gee.output
#'@export


ELCICgee <- function(models, candidate.cor.sets,data, family,r,id,time)
{
    if (!inherits(family, "family"))
        stop("dist must be an object of type 'family'")
    ## dist should be passed directly to glm, but I don't here
    ## because I don't want to change the package functions
    dist <- family$family
    if (!is.list(models))
        models <- list(models)
    y <- model.response(model.frame(models[[1]], data))
    l <- lapply(models, function(li) m <- model.matrix(li, data))
    name.var.sets <- lapply(l, colnames)
    dat <- do.call(cbind, l)
    dat <- dat[,!duplicated(colnames(dat))]
    ci <- ELCIC.gee(x=dat,y=y,r=r,id=id,time=time,
                    name.var.sets=name.var.sets,dist=dist,candidate.cor.sets=candidate.cor.sets)

    attr(ci, "formula") <- models
    attr(ci, "type") <- "gee"
    class(ci) <- "elcic"

    print.elcic.gee(ci,candidate.cor.sets)

    ##fit gee based on ELCIC
    fit.gee<-geeglm(models[[which(ci==min(ci),arr.ind = TRUE)[2]]],data=data,family =family,id=id,
           corstr = candidate.cor.sets[which(ci==min(ci),arr.ind = TRUE)[1]])

    list(model.selection=ci,gee.output=fit.gee)
}


#'@title Model selection based on ELCIC under the syntax of WGEE (Main function).
#'@description The function \code{\link{ELCICwgee}} provides the model selection under the syntax of wgeesel package.
#'@usage ELCICwgee(models, candidate.cor.sets, data, model_mis, family,r,id,time)
#'@param models A list of formulas. See corresponding documentation to wgeesel.
#'@param candidate.cor.sets A vector containing candidate correlation structures. It can be any subset of c("independence","exchangeable", "ar1").
#'@param data A data frame containing the variables in both the main model and the missing model.
#'@param model_mis A formulas used in the missing data model.
#'@param family A description of the error distribution and link function to be used in the model.
#'The details are given under "Details".
#'@param r A vector indicating the observation of data: 1 for observed records (both outcome and covariates are observed for a given subject),
#'and 0 for unobserved records.
#'The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject.
#'@return A list with two items: model selection result based on ELCIC;
#'An object of "geeglm" based on the selected model.
#'@details Three common distributions are considered: "gaussian", "poisson", "binomial".
#'For the current package, the identical link is consider for "gaussian" distribution;
#'the log link is consider for "poisson" distribution; the logit link is consider for "binomial" distribution;
#'
#'
#'
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'family<-binomial()
#'r<-wgeesimdata$obs_ind
#'id<-wgeesimdata$id
#'time<-3
#'dat <- data.frame(y=wgeesimdata$y, wgeesimdata$x,x_mis1=wgeesimdata$x_mis[,2])
#'models <- list(y~x1, y~x1+x2)
#'model_mis<-r~x_mis1
#'candidate.cor.sets<-c("exchangeable")
#'output<-ELCICwgee(models, candidate.cor.sets,data=dat,model_mis,family,r,id,time)
#'output$model.selection
#'output$gee.output
#'
#'
#'@export


ELCICwgee <- function(models, candidate.cor.sets,data, model_mis, family,r,id,time)
{
    if (!inherits(family, "family"))
        stop("dist must be an object of type 'family'")
    ## dist should be passed directly to glm, but I don't here
    ## because I don't want to change the package functions
    dist <- family$family
    if (!is.list(models))
        models <- list(models)
    y <- model.response(model.frame(models[[1]], data,na.action=NULL))
    l <- lapply(models, function(li)
        {
        m <- model.frame(li, data,na.action=NULL)
        m <- model.matrix(li, m)
    })
    l_mis <- model.matrix(model_mis, data)
    name.var.sets <- lapply(l, colnames)
    name.var.sets.mis <- colnames(l_mis)

    dat <- do.call(cbind, l)
    dat <- dat[,!duplicated(colnames(dat))]
    dat_mis<-l_mis
    ci <- ELCIC.wgee(x=dat,y,x_mis=dat_mis,r=r,id=id,time=time,name.var.sets=name.var.sets,
                     dist=dist,candidate.cor.sets=candidate.cor.sets)

    attr(ci, "formula") <- models
    attr(ci, "type") <- "wgee"
    class(ci) <- "elcic"

    print.elcic.wgee(ci,candidate.cor.sets)

    ##fit wgee based on ELCIC
    data$r<-r
    fit.wgee<-wgee(models[[which(ci==min(ci),arr.ind = TRUE)[2]]],data=data,id=id,family=dist,
              corstr= candidate.cor.sets[which(ci==min(ci),arr.ind = TRUE)[1]],scale = NULL,mismodel =model_mis)

    list(model.selection=ci,
         wgee.output=fit.wgee
         )
}


