#'@title Variable selection in generalized linear models (GLMs)
#'@description This function will provide values of several model selection criteria including AIC, BIC, GIC, and ELCIC.
#'@usage ELCIC.glm(x, y, index.var=NULL, name.var = NULL, dist)
#'@param x A matrix containing covariates. The first column should contain all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param index.var A vector containing index corresponding to candidate covariates (including the intercept).
#'@param name.var A vector containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'
#'@return A vector containing information criteria including ELCIC, AIC, BIC, and GIC.
#'
#'@examples
#'## tests
#'# load data
#'data(glmtoydata)
#'x=glmtoydata$x
#'y=glmtoydata$y
#'#candidate model index
#'name.var=c("intercept","x1","x2")
#'index.var=c(1,2,3)
#'criteria<-ELCIC.glm(x, y, index.var =index.var, name.var = NULL, dist="poisson")
#'criteria
#'
#'@export
#'@import MASS
#'@importFrom stats as.formula glm rbinom rnbinom rnorm rpois

ELCIC.glm<-function (x,y,index.var=NULL,name.var=NULL,dist)
{
    if(!is.matrix(x))
    {stop("x should be in a matrix format")}else if(!is.vector(y)){stop("y should be in a vector format")}else
        if (!dist %in% c("gaussian","binomial","poisson")){stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")}

    if(!is.null(name.var))
    {
        if(length(name.var)>ncol(x)|length(unique(name.var))<length(name.var)){stop("Invalid candidate model provided")}
        index.var<-match(name.var,colnames(x))
    }
        if(length(index.var)>ncol(x)|max(index.var)>ncol(x)|length(unique(index.var))<length(index.var)){stop("Invalid candidate model provided")}
        samplesize<-nrow(x)
        betahat<-rep(0,ncol(x))
        x_candidate<-x[,index.var]
        reg<-glm(y~x_candidate-1, family = dist)
        betahat[index.var]<-reg$coefficients
        p<-length(reg$coefficients)


    lambda<-lambda.find.glm(x=x,y=y,betahat=betahat, dist=dist)
    Z<-ee.glm(x,y,betahat,dist)

    likelihood<-apply(Z,2,function(x) {
    2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
    ELCIC<-likelihood+p*log(samplesize)

    #get aic, bic, and gic
    AIC<-reg$aic
    tr_matrix<-II(x=x_candidate,y=y,betahat=reg$coefficients,size=length(index.var),samplesize=samplesize, dist=dist)%*%ginv(JJ(x=x_candidate,y=y,betahat=reg$coefficients,size=length(index.var),samplesize=samplesize,dist=dist))
    GIC<-AIC-2*p+2*sum(diag(tr_matrix))
    BIC<-AIC-2*p+p*log(samplesize)
    return(c(ELCIC=ELCIC, AIC=AIC, BIC=BIC, GIC=GIC))
}





#'@title Calculate ELCIC value for a given candidate model under GEE framework without missingness
#'@description This is the function to calculate ELCIC value for a given candidate model. It is able to simultaneously evaluate mean model and working correlation structure.
#'@usage ELCIC.gee(x, y, r, id, time, index.var=NULL, name.var = NULL, dist, corstr, joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param time The number of observations for each subject
#'@param index.var A vector containing index corresponding to candidate covariates.
#'@param name.var A vector containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'@param joints A logic value for joint selection of marginal mean and working correlation structure. The default is TRUE.
#'
#'@return A ELCIC value for a given candidate model.
#'
#'@examples
#'## tests
#'# load data
#'data(geetoydata)
#'x=geetoydata$x
#'y=geetoydata$y
#'id=geetoydata$id
#'r=rep(1,nrow(x))
#'time=3
#'corstr="exchangeable"
#'dist="poisson"
#'criteria<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,index.var=c(1,2,3),
#'            name.var=NULL,dist=dist,corstr=corstr)
#'criteria
#'
#'@export
#'@importFrom geepack geeglm

ELCIC.gee<-function(x,y,r,id,time,index.var=NULL,name.var=NULL,dist,corstr,joints=TRUE)
{
    if(!is.matrix(x))
    {stop("x should be in a matrix format")}else if(!is.vector(y)){stop("y should be in a vector format")}else
        if (!dist %in% c("gaussian","binomial","poisson")){stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")}
    if (!corstr %in% c("ar1","exchangeable","independence")){stop("Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")}

    if(!is.null(name.var))
    {
        if(length(name.var)>ncol(x)|length(unique(name.var))<length(name.var)){stop("Invalid candidate model provided")}
        index.var<-match(name.var,colnames(x))
    }
    if(length(index.var)>ncol(x)|max(index.var)>ncol(x)|length(unique(index.var))<length(index.var)){stop("Invalid candidate model provided")}

    if(joints)
    {
    samplesize<-length(unique(id))
    betahat<-rep(0,ncol(x))
    x_candidate<-(x[,index.var])
    colnames(x_candidate)<-rep(1:ncol(x_candidate))
    y<-as.matrix(y,col=1)
    data<-data.frame(y,x_candidate)

        fit<-geeglm(y~x_candidate-1,data=data,family =dist,id=id,corstr = corstr)
        betahat[index.var]<-fit$coefficients
        pbeta<-length(fit$coefficients)

        if(corstr=="independence")
        {
            ro=0
            p<-pbeta
        }else{
            ro<-unlist(summary(fit)$corr[1]) #correlation coefficients
            p<-pbeta+1
        }
        phihat<-unlist(summary(fit)$dispersion[1])
        Z<-as.matrix(ee.gee(y=y,x=x,r=r,id=id,beta=betahat,ro=ro,phi=phihat,dist=dist,corstr=corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.gee(x=x,y=y,id=id,betahat=betahat,r=r,dist=dist,ro=ro,phi=phihat,corstr=corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)
    return(ELCIC=ELCIC)
    }else{
        samplesize<-length(unique(id))
        betahat<-rep(0,ncol(x))
        x_candidate<-(x[,index.var])
        colnames(x_candidate)<-rep(1:ncol(x_candidate))
        y<-as.matrix(y,col=1)
        data<-data.frame(y,x_candidate)

        fit<-geeglm(y~x_candidate-1,data=data,family =dist,id=id,corstr = corstr)
        betahat[index.var]<-fit$coefficients
        p<-length(fit$coefficients)

        # if(corstr=="independence")
        # {
        #     ro=0
        #     p<-pbeta
        # }else{
        #     ro<-unlist(summary(fit)$corr[1]) #correlation coefficients
        #     p<-pbeta+1
        # }
        phihat<-unlist(summary(fit)$dispersion[1])
        Z<-as.matrix(ee.gee.onlymean(y=y,x=x,r=r,id=id,beta=betahat,ro=ro,phi=phihat,dist=dist,corstr=corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.gee.onlymean(x=x,y=y,id=id,betahat=betahat,r=r,dist=dist,ro=ro,phi=phihat,corstr=corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)
        return(ELCIC=ELCIC)
    }
}




#'@title Calculate ELCIC value for a given candidate model under WGEE framework with data missing at random
#'@description This is the function to calculate ELCIC value for a given candidate model. It is able to simultaneously evaluate mean model and working correlation structure. The data is missing at random.
#'@usage ELCIC.wgee(x,y,x_mis,r,id,time,index.var=NULL,name.var=NULL,dist,corstr,joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records.
#'@param id A vector indicating subject id.
#'@param time The number of observations for each subject
#'@param index.var A vector containing index corresponding to candidate covariates.
#'@param name.var A vector containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'@param joints A logic value for joint selection of marginal mean and working correlation structure. The default is TRUE.

#'@return A matrix containing values of calculated estimating equations.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'corstr="exchangeable"
#'dist="binomial"
#'x=wgeetoydata$x
#'y=wgeetoydata$y
#'x_mis=wgeetoydata$x_mis
#'r=wgeetoydata$obs_ind
#'id=wgeetoydata$id
#'time=3
#'index.var=c(1,2,3)
#'ELCIC_value=ELCIC.wgee(x,y,x_mis,r,id,time,index.var,name.var=NULL,
#'                      dist,corstr,joints=TRUE)
#'ELCIC_value

#'@export
#'@importFrom wgeesel wgee
#'
ELCIC.wgee<-function(x,y,x_mis,r,id,time,index.var=NULL,name.var=NULL,dist,corstr,joints=TRUE)
{
    if(!is.matrix(x))
    {stop("x should be in a matrix format")}else if(!is.vector(y)){stop("y should be in a vector format")}else if(!is.matrix(x_mis)){stop("x_mis should be in a matrix format")}else
        if (!dist %in% c("gaussian","binomial","poisson")){stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")}
    if (!corstr %in% c("ar1","exchangeable","independence")){stop("Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")}

    if(!is.null(name.var))
    {
        if(length(name.var)>ncol(x)|length(unique(name.var))<length(name.var)){stop("Invalid candidate model provided")}
        index.var<-match(name.var,colnames(x))
    }
    if(length(index.var)>ncol(x)|max(index.var)>ncol(x)|length(unique(index.var))<length(index.var)){stop("Invalid candidate model provided")}

    if(joints)
    {
        samplesize<-length(unique(id))
        beta<-rep(0,ncol(x))
        x_candidate<-(x[,index.var])
        colnames(x_candidate)<-rep(1:ncol(x_candidate))
        y<-as.matrix(y,col=1)
        data_wgee<-data.frame(y,x_candidate,r,x_mis)

        mismodel_formula<-paste("r~",colnames(x_mis)[2])
        if(ncol(x_mis)>2)
        {
        for(jj in 3:ncol(x_mis))
        {
        mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
        }
        }
        mismodel_formula<-as.formula(mismodel_formula)
        fit=wgee(y~x_candidate-1,data_wgee,id,family=dist,corstr =corstr,scale = NULL,mismodel =mismodel_formula)

        beta[index.var]<-as.vector(summary(fit)$beta)
        pbeta<-length(summary(fit)$beta)

        if(corstr=="independence")
        {
            ro=0
            p<-pbeta
        }else{
            ro<-summary(fit)$corr
            p<-pbeta+1
        }
        phi=summary(fit)$phi
        gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
        pi=prob.obs(x_mis,gamma)
        Z<-as.matrix(ee.wgee(y,x,r, pi,id,time=3,beta,ro,phi,dist,corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.wgee(y,x,r, pi,id,time=3,beta,ro,phi,dist,corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)
        return(ELCIC=ELCIC)
    }else{
        samplesize<-length(unique(id))
        beta<-rep(0,ncol(x))
        x_candidate<-(x[,index.var])
        colnames(x_candidate)<-rep(1:ncol(x_candidate))
        y<-as.matrix(y,col=1)
        data_wgee<-data.frame(y,x_candidate,r,x_mis)

        mismodel_formula<-paste("r~",colnames(x_mis)[2])
        if(ncol(x_mis)>2)
        {
            for(jj in 3:ncol(x_mis))
            {
                mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
            }
        }
        mismodel_formula<-as.formula(mismodel_formula)
        fit=wgee(y~x_candidate-1,data_wgee,id,family=dist,corstr =corstr,scale = NULL,mismodel =mismodel_formula)

        beta[index.var]<-as.vector(summary(fit)$beta)
        p<-length(summary(fit)$beta)

        # if(corstr=="independence")
        # {
        #     ro=0
        #     p<-pbeta
        # }else{
        #     ro<-summary(fit)$corr
        #     p<-pbeta+1
        # }
        phi=summary(fit)$phi
        gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
        pi=prob.obs(x_mis,gamma)
        Z<-as.matrix(ee.wgee.onlymean(y,x,r, pi,id,time=3,beta,ro,phi,dist,corstr))
        # epi<-1/samplesize
        # model<-function(lambda)
        # {
        #     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
        #     {x/(1+t(lambda)%*%x)}
        #     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,samplesize)
        # }
        lambda<-lambda.find.wgee.onlymean(y,x,r, pi,id,time=3,beta,ro,phi,dist,corstr)
        # lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1))$root
        likelihood<-apply(Z,2,function(x) {
            2*log(1+t(lambda)%*%x)})%*%rep(1,samplesize)
        ELCIC<-likelihood+p*log(samplesize)
        return(ELCIC=ELCIC)
    }
}






#'@title The whole variable selection procedure in GLM
#'@description This function provides the overall procedure for variable selection in GLM.
#'@usage ELCIC.glm.procedure(x,y,candidate.sets,name.var.sets=NULL,dist)
#'@param x A matrix containing covariates. The first column should contain all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param candidate.sets A list containing index corresponding to candidate covariates in each candidate model.
#'@param name.var.sets A list containing names of candidate covariates coresponding to each candidate model. The names should be subset of column names of the x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'
#'@return A matrix with each element containing ELCIC value for each candidate model (in columns) and (in rows)
#'
#'@examples
#'## tests
#'# load data
#'data(glmtoydata)
#'x=glmtoydata$x
#'y=glmtoydata$y
#'#candidate model index
#'candidate.sets<-list(c(1,2),c(1,2,3),c(1,2,3,4))
#'criteria<-ELCIC.glm.procedure(x, y, candidate.sets, name.var.sets = NULL, dist="poisson")
#'criteria
#'
#'@export

ELCIC.glm.procedure<-function(x,y,candidate.sets,name.var.sets=NULL,dist)
{
    if(!is.null(name.var.sets))
    {
    criterion.all<-rep()
    for (i in 1:length(name.var.sets))
    {
        criterion<-ELCIC.glm(x=x,y=y,index.var=NULL,name.var=name.var.sets[[i]],dist=dist)
        criterion.all<-cbind(criterion.all,criterion)
    }
    criterion.all
    }else{
        criterion.all<-rep()
        for (i in 1:length(candidate.sets))
        {
            criterion<-ELCIC.glm(x=x,y=y,index.var=candidate.sets[[i]],name.var=NULL,dist=dist)
            criterion.all<-cbind(criterion.all,criterion)
        }
        criterion.all
    }
}





#'@title The whole procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness
#'@description This function provides the overall procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness.
#'@usage ELCIC.gee.procedure(x,y,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,
#'       candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param time The number of observations for each subject
#'@param candidate.sets A list containing index corresponding to candidate covariates.
#'@param name.var.sets A list containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param candidate.cor.sets A vector containing condidate correlation structures. When joint=TRUE, it is c("independence","exchangeable", "ar1") as default. When joint=FALSE, it should be one of c("independence","exchangeable", "ar1").
#'@param joints A logic value for joint selection of marginal mean and working correlation structure. The default is TRUE.
#'
#'@return A matrix with each element containing ELCIC value for each candidate model.
#'
#'@examples
#'## tests
#'# load data
#'data(geetoydata)
#'x=geetoydata$x
#'y=geetoydata$y
#'id=geetoydata$id
#'r=rep(1,nrow(x))
#'time=3
#'candidate.sets<-list(c(1,2),c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'dist="poisson"
#'criterion.elcic<-ELCIC.gee.procedure(x=x,y=y,r=r,id=id,time=time,candidate.sets=candidate.sets,
#'                 name.var.sets=NULL,dist=dist,candidate.cor.sets=candidate.cor.sets)
#'criterion.elcic
#'
#'@export

ELCIC.gee.procedure<-function(x,y,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence",
                                                                                                   "exchangeable", "ar1"), joints=TRUE)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            criterion.elcic<-rep()
            for(j in 1:length(candidate.cor.sets))
            {
                corstr<-candidate.cor.sets[j]
                criterion.all<-rep()
                for (i in 1:length(name.var.sets))
                {
                    criterion<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,index.var=NULL,name.var=name.var.sets[[i]],dist=dist,corstr=corstr)
                    criterion.all<-c(criterion.all,criterion)
                }
                criterion.elcic<-rbind(criterion.elcic,criterion.all)
                #print(j)
            }
            colnames(criterion.elcic)<-name.var.sets
            rownames(criterion.elcic)<-candidate.cor.sets
            criterion.elcic
        }else{
    criterion.elcic<-rep()
    for(j in 1:length(candidate.cor.sets))
    {
        corstr<-candidate.cor.sets[j]
        criterion.all<-rep()
        for (i in 1:length(candidate.sets))
        {
            criterion<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,index.var=candidate.sets[[i]],name.var=NULL,dist=dist,corstr=corstr)
            criterion.all<-c(criterion.all,criterion)
        }
        criterion.elcic<-rbind(criterion.elcic,criterion.all)
        #print(j)
    }
    colnames(criterion.elcic)<-candidate.sets
    rownames(criterion.elcic)<-candidate.cor.sets
    criterion.elcic
        }
    }else{
        if(!is.null(name.var.sets))
        {
                criterion.elcic<-rep()
                corstr<-candidate.cor.sets
                for (i in 1:length(name.var.sets))
                {
                    criterion<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,index.var=NULL,name.var=name.var.sets[[i]],dist=dist,corstr=corstr)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-name.var.sets
            criterion.elcic
        }else{
                criterion.elcic<-rep()
                corstr<-candidate.cor.sets
                for (i in 1:length(candidate.sets))
                {
                    criterion<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,index.var=candidate.sets[[i]],name.var=NULL,dist=dist,corstr=corstr)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-candidate.sets
            criterion.elcic
        }
    }
}




#'@title The whole procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness
#'@description This function provides the overall procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness.
#'@usage ELCIC.wgee.procedure(x,y,x_mis,r,id,time,candidate.sets=NULL,name.var.sets=NULL,
#'      dist,candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records.
#'@param id A vector indicating subject id.
#'@param time The number of observations for each subject
#'@param candidate.sets A list containing index corresponding to candidate covariates.
#'@param name.var.sets A list containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param candidate.cor.sets A vector containing condidate correlation structures. When joint=TRUE, it is c("independence","exchangeable", "ar1") as default. When joint=FALSE, it should be one of c("independence","exchangeable", "ar1").
#'@param joints A logic value for joint selection of marginal mean and working correlation structure. The default is TRUE.

#'@return A matrix with each element containing ELCIC value for each candidate model.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'dist="binomial"
#'x=wgeetoydata$x
#'y=wgeetoydata$y
#'x_mis=wgeetoydata$x_mis
#'r=wgeetoydata$obs_ind
#'id=wgeetoydata$id
#'time=3
#'candidate.sets<-list(c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'criterion.elcic<-ELCIC.wgee.procedure(x,y,x_mis,r,id,time,candidate.sets,name.var.sets=NULL,
#'                                     dist,candidate.cor.sets,joints=TRUE)
#'criterion.elcic
#'
#'@export

ELCIC.wgee.procedure<-function(x,y,x_mis,r,id,time,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            criterion.elcic<-rep()
            for(j in 1:length(candidate.cor.sets))
            {
                corstr<-candidate.cor.sets[j]
                criterion.all<-rep()
                for (i in 1:length(name.var.sets))
                {
                    criterion<-ELCIC.wgee(x,y,x_mis,r,id,time,index.var=NULL,name.var=name.var.sets[[i]],dist,corstr,joints=TRUE)
                    criterion.all<-c(criterion.all,criterion)
                }
                criterion.elcic<-rbind(criterion.elcic,criterion.all)
                #print(j)
            }
            colnames(criterion.elcic)<-name.var.sets
            rownames(criterion.elcic)<-candidate.cor.sets
            criterion.elcic
        }else{
            criterion.elcic<-rep()
            for(j in 1:length(candidate.cor.sets))
            {
                corstr<-candidate.cor.sets[j]
                criterion.all<-rep()
                for (i in 1:length(candidate.sets))
                {
                    criterion<-ELCIC.wgee(x,y,x_mis,r,id,time,index.var=candidate.sets[[i]],name.var=NULL,dist,corstr,joints=TRUE)
                    criterion.all<-c(criterion.all,criterion)
                }
                criterion.elcic<-rbind(criterion.elcic,criterion.all)
                #print(j)
            }
            colnames(criterion.elcic)<-candidate.sets
            rownames(criterion.elcic)<-candidate.cor.sets
            criterion.elcic
        }
    }else{
        if(!is.null(name.var.sets))
        {
                criterion.elcic<-rep()
                corstr<-candidate.cor.sets
                for (i in 1:length(name.var.sets))
                {
                    criterion<-ELCIC.wgee(x,y,x_mis,r,id,time,index.var=NULL,name.var=name.var.sets[[i]],dist,corstr,joints=FALSE)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-name.var.sets
            criterion.elcic
        }else{
            criterion.elcic<-rep()
                corstr<-candidate.cor.sets
                for (i in 1:length(candidate.sets))
                {
                    criterion<-ELCIC.wgee(x,y,x_mis,r,id,time,index.var=candidate.sets[[i]],name.var=NULL,dist,corstr,joints=FALSE)
                    criterion.elcic<-c(criterion.elcic,criterion)
                }

            names(criterion.elcic)<-candidate.sets
            criterion.elcic
        }
    }
}



#'@title Joint selection of marginal mean and correlation structures in longitudinal data based on QIC
#'@description This function provides the Joint selection of marginal mean and correlation structures in longitudinal data based on QIC.
#'@usage QICc.procedure(x,y,id,dist,candidate.sets=NULL, name.var.sets=NULL,
#'    candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param id A vector indicating subject id.
#'@param candidate.sets A list containing index corresponding to candidate covariates.
#'@param name.var.sets A list containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param candidate.cor.sets A vector containing condidate correlation structures. When joint=TRUE, it is c("independence","exchangeable", "ar1") as default. When joint=FALSE, it should be one of c("independence","exchangeable", "ar1").
#'@param joints A logic value for joint selection of marginal mean and working correlation structure. The default is TRUE.
#'
#'@return A vector with each element containing QIC value for each candidate model. The row name of this vector is the selected correlation structure.
#'
#'@examples
#'## tests
#'# load data
#'data(geetoydata)
#'x=geetoydata$x
#'y=geetoydata$y
#'id=geetoydata$id
#'r=rep(1,nrow(x))
#'time=3
#'candidate.sets<-list(c(1,2),c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'dist="poisson"
#'criterion.qic<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=candidate.sets,
#'                     name.var.sets=NULL,candidate.cor.sets=candidate.cor.sets)
#'criterion.qic
#'
#'@export
#'@importFrom geepack QIC

QICc.procedure<-function (x,y,id,dist,candidate.sets,name.var.sets=NULL,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
    data<-data.frame(y,x)
    modelQ<-rep()
    for(i in 1:length(candidate.cor.sets))
    {
        fit<-geeglm(y~x-1,data=data,family =dist,id=id,corstr = candidate.cor.sets[i])
        modelQ[i]<-as.numeric(QIC(fit)[4])##look at QIC and QICu
    }
    QICcorstr<-candidate.cor.sets[which.min(modelQ)]

    criterion.mean<-rep()
    for (i in 1:length(name.var.sets))
    {
        x.candidate<-x[,name.var.sets[[i]]]
        data<-data.frame(y,x.candidate)
        fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
        criterion.mean<-c(criterion.mean,QIC(fit)[1])
    }
    criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
    rownames(criterion.mean)<-QICcorstr
    colnames(criterion.mean)<-name.var.sets
    criterion.mean
        }else{
            data<-data.frame(y,x)
            modelQ<-rep()
            for(i in 1:length(candidate.cor.sets))
            {
                fit<-geeglm(y~x-1,data=data,family =dist,id=id,corstr = candidate.cor.sets[i])
                modelQ[i]<-as.numeric(QIC(fit)[4])##look at QIC and QICu
            }
            QICcorstr<-candidate.cor.sets[which.min(modelQ)]

            criterion.mean<-rep()
            for (i in 1:length(candidate.sets))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                data<-data.frame(y,x.candidate)
                fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
                criterion.mean<-c(criterion.mean,QIC(fit)[1])
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }else{
        if(!is.null(name.var.sets))
        {
            data<-data.frame(y,x)
            QICcorstr<-candidate.cor.sets

            criterion.mean<-rep()
            for (i in 1:length(name.var.sets))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                data<-data.frame(y,x.candidate)
                fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
                criterion.mean<-c(criterion.mean,QIC(fit)[1])
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            data<-data.frame(y,x)
            QICcorstr<-candidate.cor.sets

            criterion.mean<-rep()
            for (i in 1:length(candidate.sets))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                data<-data.frame(y,x.candidate)
                fit<-geeglm(y~x.candidate-1,data=data,family =dist,id=id,corstr = QICcorstr)
                criterion.mean<-c(criterion.mean,QIC(fit)[1])
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }
}




#'@title The whole MLIC procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness
#'@description This function provides the overall MLIC procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness.
#'@usage MLIC.wgee.procedure(x,y,x_mis,r,id,candidate.sets=NULL, name.var.sets=NULL,dist,
#'       candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records.
#'@param id A vector indicating subject id.
#'@param candidate.sets A list containing index corresponding to candidate covariates.
#'@param name.var.sets A list containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param candidate.cor.sets A vector containing condidate correlation structures. When joint=TRUE, it is c("independence","exchangeable", "ar1") as default. When joint=FALSE, it should be one of c("independence","exchangeable", "ar1").
#'@param joints A logic value for joint selection of marginal mean and working correlation structure. The default is TRUE.

#'@return A vector with each element containing MLIC value for each candidate model. The row name of this vector is the selected correlation structure.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'dist="binomial"
#'x=wgeetoydata$x
#'y=wgeetoydata$y
#'x_mis=wgeetoydata$x_mis
#'r=wgeetoydata$obs_ind
#'id=wgeetoydata$id
#'candidate.sets<-list(c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'criterion.mlic<-MLIC.wgee.procedure(x,y,x_mis,r,id,candidate.sets,
#'             name.var.sets=NULL,dist,candidate.cor.sets,joints=FALSE)
#'criterion.mlic

#'@export
#'@importFrom wgeesel wgee data_sim MLIC.gee QICW.gee

MLIC.wgee.procedure<-function(x,y,x_mis,r,id,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in 1:length(candidate.cor.sets))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-MLIC.gee(fit,fit)$MLICc##look at MLICc
            }
            MLICcorstr<-candidate.cor.sets[which.min(modelQ)]

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in 1:length(name.var.sets))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-MLICcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in 1:length(candidate.cor.sets))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-MLIC.gee(fit,fit)$MLICc##look at MLICc
            }
            MLICcorstr<-candidate.cor.sets[which.min(modelQ)]

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in 1:length(candidate.sets))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-MLICcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }else{
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)

            MLICcorstr<-candidate.cor.sets

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in 1:length(name.var.sets))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            colnames(criterion.mean)<-name.var.sets
            rownames(criterion.mean)<-candidate.cor.sets
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)

            MLICcorstr<-candidate.cor.sets

            #fit a full model
            fitf<-wgee(y~x-1,data_wgee,id,family=dist,corstr = MLICcorstr,scale = NULL,mismodel =mismodel_formula)

            criterion.mean<-rep()
            for (i in 1:length(candidate.sets))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =MLICcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,MLIC.gee(fitm,fitf)$MLIC)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            colnames(criterion.mean)<-candidate.sets
            rownames(criterion.mean)<-candidate.cor.sets
            criterion.mean
        }
    }
}




#'@title The whole QICW procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness
#'@description This function provides the overall QICW procedure for joint selection of mean structure and correlation structure in longitudinal data without missingness.
#'@usage QICW.wgee.procedure(x,y,x_mis,r,id,candidate.sets,name.var.sets=NULL,
#'      dist,candidate.cor.sets=c("independence","exchangeable", "ar1"), joints=TRUE)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records.
#'@param id A vector indicating subject id.
#'@param candidate.sets A list containing index corresponding to candidate covariates.
#'@param name.var.sets A list containing names of candidate covariates. The names should be subset of column names of x matrix.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param candidate.cor.sets A vector containing condidate correlation structures. When joint=TRUE, it is c("independence","exchangeable", "ar1") as default. When joint=FALSE, it should be one of c("independence","exchangeable", "ar1").
#'@param joints A logic value for joint selection of marginal mean and working correlation structure. The default is TRUE.

#'@return A vector with each element containing QICW value for each candidate model. The row name of this vector is the selected correlation structure.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'dist="binomial"
#'x=wgeetoydata$x
#'y=wgeetoydata$y
#'x_mis=wgeetoydata$x_mis
#'r=wgeetoydata$obs_ind
#'id=wgeetoydata$id
#'candidate.sets<-list(c(1,2,3))
#'candidate.cor.sets<-c("exchangeable")
#'criterion.qicw<-QICW.wgee.procedure(x,y,x_mis,r,id,candidate.sets,
#'           name.var.sets=NULL,dist,candidate.cor.sets,joints=FALSE)
#'criterion.qicw
#'
#'@export
#'@importFrom wgeesel wgee data_sim MLIC.gee QICW.gee

QICW.wgee.procedure<-function(x,y,x_mis,r,id,candidate.sets=NULL,name.var.sets=NULL,dist,candidate.cor.sets=c("independence","exchangeable", "ar1"),joints=TRUE)
{
    if(joints)
    {
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in 1:length(candidate.cor.sets))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-QICW.gee(fit)$QICWr##look at QICWr
            }
            QICWcorstr<-candidate.cor.sets[which.min(modelQ)]

            criterion.mean<-rep()
            for (i in 1:length(name.var.sets))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            modelQ<-rep()
            for(i in 1:length(candidate.cor.sets))
            {
                fit<-wgee(y~x-1,data_wgee,id,family=dist,corstr =candidate.cor.sets[[i]],scale = NULL,mismodel =mismodel_formula)
                modelQ[i]<-QICW.gee(fit)$QICWr##look at QICWr
            }
            QICWcorstr<-candidate.cor.sets[which.min(modelQ)]

            criterion.mean<-rep()
            for (i in 1:length(candidate.sets))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }else{
        if(!is.null(name.var.sets))
        {
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            QICWcorstr<-candidate.cor.sets

            criterion.mean<-rep()
            for (i in 1:length(name.var.sets))
            {
                x.candidate<-x[,name.var.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-name.var.sets
            criterion.mean
        }else{
            y<-as.matrix(y,col=1)
            data_wgee<-data.frame(y,x,r,x_mis)

            mismodel_formula<-paste("r~",colnames(x_mis)[2])
            if(ncol(x_mis)>2)
            {
                for(jj in 3:ncol(x_mis))
                {
                    mismodel_formula<-paste(mismodel_formula,"+",colnames(x_mis)[jj])
                }
            }
            mismodel_formula<-as.formula(mismodel_formula)
            QICWcorstr<-candidate.cor.sets

            criterion.mean<-rep()
            for (i in 1:length(candidate.sets))
            {
                x.candidate<-x[,candidate.sets[[i]]]
                fitm<-wgee(y~x.candidate-1,data_wgee,id,family=dist,corstr =QICWcorstr,scale = NULL,mismodel =mismodel_formula)
                criterion.mean<-c(criterion.mean,QICW.gee(fitm)$QICWr)
            }
            criterion.mean<-t(as.matrix(criterion.mean,nrow=1))
            rownames(criterion.mean)<-QICWcorstr
            colnames(criterion.mean)<-candidate.sets
            criterion.mean
        }
    }
}


