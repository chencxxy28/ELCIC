

#'@title Estimating equation for ELCIC under GLM
#'@description A specified estimating equation for ELCIC under GLM. This estimating equation is used for marginal mean selection.
#'@usage ee.glm(x, y, beta, dist)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept. See more details in
#'@param y A vector containing outcomes.
#'@param beta A plug-in estimator solved by an external estimating procedure.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@note "x" and "y" should be all observed.
#'
#'@examples
#'## tests
#'# load data
#'data(glmsimdata)
#'x<-glmsimdata$x
#'y<-glmsimdata$y
#'# obtain the estimates. Note that x matrix already contains intercept.
#'fit<-glm(y~x-1,family="poisson")
#'beta<-fit$coefficients
#'ee.matrix<-ee.glm(x, y, beta, dist="poisson")
#'apply(ee.matrix,1,mean)
#'
#'@export

######for cross-sectional glm
ee.glm<-function (x,y,beta,dist)
{
    full<-ncol(x)
    samplesize<-nrow(x)
    ee<-matrix(0,nrow=full,ncol=samplesize)

    switch(dist,
           "gaussian"={    for (i in 1:samplesize)
           {
               ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(t(as.matrix(x[i,],ncol=1))%*%beta))
           }
           },
           "binomial"={ for (i in 1:samplesize)
           {
               ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(1+exp(-t(as.matrix(x[i,],ncol=1))%*%beta))^(-1))
           }
           },
           "poisson"={for (i in 1:samplesize)
           {
               ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-exp(t(as.matrix(x[i,],ncol=1))%*%beta))
           }
           },
           stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")
    )

    ee

    # if(dist=="gaussian")
    # {
    # for (i in 1:samplesize)
    # {
    #     ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(t(as.matrix(x[i,],ncol=1))%*%beta))
    # }
    # }else if (dist=="poisson")
    # {
    #     for (i in 1:samplesize)
    #     {
    #     ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-exp(t(as.matrix(x[i,],ncol=1))%*%beta))
    #     }
    # }else if (dist=="binomial")
    # {
    #     for (i in 1:samplesize)
    #     {
    #     ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(1+exp(-t(as.matrix(x[i,],ncol=1))%*%beta))^(-1))
    #     }
    # }

}




#######################################################
###### for joint selection of marginal mean and working correlation

#'@title Estimating equation for GEE without missingness or with data missing completely at random.
#'@description Calculate estimating equation from GEE in ELCIC without missingness or missing completely at random. This estimating equation is used for joint selection of marginal mean and working correlation structure.
#'@usage ee.gee(y,x,r,id,beta,rho,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones the represents the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating the observation of outcomes: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed. See more in details section.
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure.
#'@param rho A correlation coefficients obtained from an external estimation procedure, such as GEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A candidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@details If the element in argument "r" equals zero, the corresponding rows of "x" and "y" should be all zeros.
#'
#'@examples
#'## tests
#'# load data
#'data(geesimdata)
#'x<-geesimdata$x
#'y<-geesimdata$y
#'id<-geesimdata$id
#'corstr<-"exchangeable"
#'dist<-"poisson"
#'# obtain the estimates
#'library(geepack)
#'# x matrix already include the intercept column.
#'fit<-geeglm(y~x-1,data=geesimdata,family =dist,id=id,corstr = "ar1")
#'beta<-fit$coefficients
#'rho<-unlist(summary(fit)$corr[1])
#'phi<-unlist(summary(fit)$dispersion[1])
#'r<-rep(1,nrow(x))
#'ee.matrix<-ee.gee(y,x,r,id,beta,rho,phi,dist,corstr)
#'apply(ee.matrix,1,mean)
#'
#'@export
#'
######for longitudinal data gee
ee.gee<-function(y,x,r,id,beta,rho,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    time<-length(id)/n
    A<-diag(1,time,time)
    rho<-roo(rho,time,corstr)
    R<-R(rho,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    #z.col<-ncol(z)
    m_num_vector<-rep()
     for(m in 1:(time-1))
     {
         m_num<-0
         for (i in 1:n) #formula for rho
         {
         WW<-W*r[((i-1)*time+1):(i*time)]
         for (j in 1:(time-m))
         {
             m_num<-m_num+WW[j+m,j+m]*WW[j,j]
         }
         }
         m_num_vector<-c(m_num_vector,m_num)
     }

    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        #piii<-pii(pi[((i-1)*time+1):(i*time)])  #used when we have ipw
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        #WW<-W*r[((i-1)*time+1):(i*time)]*piii  #used when we have ipw
        WW<-W*r[((i-1)*time+1):(i*time)]
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        for (m in 1:(time-1)) #formula for rho
        {error<-0
        for (j in 1:(time-m))
        {
            error<-error+e[j]*e[j+m]*WW[j+m,j+m]*WW[j,j]
        }
        error<-error-rho[m]*phi*(m_num_vector[m]-p)/n #directly use phi might be more efficient
        wgeei<-rbind(wgeei,error)
        }

        wgeef<-cbind(wgeef,wgeei)
    }

    return(wgeef)
}




#'@title Estimating equation for weighted GEE (WGEE) for missing longitudinal data under the mechanism of missing at random and drop-out
#'@description Calculate estimating equation from WGEE for missing longitudinal data under the mechanism of missing at random and drop-out. This estimating equation is used for joint selection of marginal mean and "working" correlation structure.
#'@usage ee.wgee(y,x,r,pi,id,time,beta,rho,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param r A vector indicating the observation of outcomes: 1 for observed records, and 0 for unobserved records.
#'@param pi A vector containing observing probabilities across all observations.
#'@param time The number of observations for each subject.
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as WGEE.
#'@param rho A correlation coefficients obtained from an external estimation procedure, such as WGEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A candidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'library(wgeesel)
#'data_wgee<-data.frame(do.call(cbind,wgeesimdata))
#'corstr<-"exchangeable"
#'dist<-"binomial"
#'id<-data_wgee$id
#'# obtain the estimates.
#'# Note that "obs_ind" is an indicator of observations in the missing data model.
#'fit<-wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,
#'      scale = NULL,mismodel =obs_ind~x_mis1)
#'beta<-as.vector(summary(fit)$beta)
#'rho<-summary(fit)$corr
#'phi<-summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis<-wgeesimdata$x_mis
#'pi<-cond.prob(x_mis,gamma,id,time=3)
#'wgee.matrix<-ee.wgee(y=wgeesimdata$y,x=wgeesimdata$x,r=wgeesimdata$obs_ind,
#'pi=pi,id=wgeesimdata$id,time=3,beta=beta,rho=rho,phi=phi,dist=dist,corstr=corstr)
#'apply(wgee.matrix,1,mean)
#'
#'@export

ee.wgee<-function(y,x,r,pi,id,time,beta,rho,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    rho<-roo(rho,time,corstr)
    A<-diag(1,time,time)
    R<-R(rho,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    # z.col<-ncol(z)
    # rlag<-ylag(id,r,1)
    # rlag[is.na(rlag)]<-1
    # S<-matrix(rep(0),nrow = z.col, ncol = n )


    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        piii<-marg.prob(pi[((i-1)*time+1):(i*time)])
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        WW<-W*r[((i-1)*time+1):(i*time)]*piii
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        for (m in 1:(time-1)) #formula for rho
        {error<-0
        for (j in 1:(time-m))
        {
            error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        }
        error<-error-rho[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        wgeei<-rbind(wgeei,error)
        }

        wgeef<-cbind(wgeef,wgeei)
        # ss<-t(rlag[id==i]*(r[id==i]-lambda[id==i])*z[id==i,])
        # S[,i]<-as.matrix(apply(ss,1,sum))
    }
    return(wgeef)
}




#######################################################
###### only for marginal mean selection

#'@title Estimating equation of marginal mean for GEE without missingness or missing completely at random
#'@description Calculate estimating equation from GEE in ELCIC. This estimating equation is used for marginal mean selection.
#'@usage ee.gee.mean(y,x,r,id,beta,rho,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating the observation of outcomes: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed. See more in details section.
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as GEE.
#'@param rho A correlation coefficients obtained from an external estimation procedure, such as GEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A candidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@details If the element in argument "r" equals zero, the corresponding rows of "x" and "y" should be all zeros.
#'
#'@note corstr should be prespecified.
#'
#'@examples
#'## tests
#'# load data
#'data(geesimdata)
#'x<-geesimdata$x
#'y<-geesimdata$y
#'id<-geesimdata$id
#'corstr<-"exchangeable"
#'dist<-"poisson"
#'# obtain the estimates
#'library(geepack)
#'fit<-geeglm(y~x-1,data=geesimdata,family =dist,id=id,corstr = corstr)
#'beta<-fit$coefficients
#'rho<-unlist(summary(fit)$corr[1])
#'phi<-unlist(summary(fit)$dispersion[1])
#'r<-rep(1,nrow(x))
#'ee.matrix<-ee.gee.mean(y,x,r,id,beta,rho,phi,dist,corstr)
#'apply(ee.matrix,1,mean)
#'
#'@export
ee.gee.mean<-function(y,x,r,id,beta,rho,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    time<-length(id)/n
    A<-diag(1,time,time)
    rho<-roo(rho,time,corstr)
    R<-R(rho,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    #z.col<-ncol(z)


    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        #piii<-pii(pi[((i-1)*time+1):(i*time)])  #used when we have ipw
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        #WW<-W*r[((i-1)*time+1):(i*time)]*piii  #used when we have ipw
        WW<-W*r[((i-1)*time+1):(i*time)]
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        # e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        # for (m in 1:(time-1)) #formula for rho
        # {error<-0
        # for (j in 1:(time-m))
        # {
        #     error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        # }
        # error<-error-rho[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        # wgeei<-rbind(wgeei,error)
        # }
        #
         wgeef<-cbind(wgeef,wgeei)
    }

    return(wgeef)
}




#'@title Estimating equation for marginal mean under WGEE for missing longitudinal data under the mechanism of missing at random and drop-out
#'@description Calculate estimating function from WGEE. This estimating function is used for marginal mean selection.
#'@usage ee.wgee.mean(y,x,r,pi,id,time,beta,rho,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param r A vector indicating the observation of outcomes: 1 for observed records, and 0 for unobserved records.
#'@param pi A vector containing observing probabilities across all observations.
#'@param time The number of observations for each subject
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as WGEE.
#'@param rho A correlation coefficients obtained from an external estimation procedure, such as WGEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A candidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@note corstr should be prespecified.
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'library(wgeesel)
#'data_wgee<-data.frame(do.call(cbind,wgeesimdata))
#'corstr<-"exchangeable"
#'dist<-"binomial"
#'id<-data_wgee$id
#'# obtain the estimates.
#'# Note that "obs_ind" is an indicator of observations in the missing data model.
#'fit<-wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,
#'      scale = NULL,mismodel =obs_ind~x_mis1)
#'beta<-as.vector(summary(fit)$beta)
#'rho<-summary(fit)$corr
#'phi<-summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis<-wgeesimdata$x_mis
#'pi<-cond.prob(x_mis,gamma,id,time=3)
#'wgee.matrix<-ee.wgee.mean(y=wgeesimdata$y,x=wgeesimdata$x,r=wgeesimdata$obs_ind,
#'pi=pi,id=wgeesimdata$id,time=3,beta=beta,rho=rho,phi=phi,dist=dist,corstr=corstr)
#'apply(wgee.matrix,1,mean)
#'@export

ee.wgee.mean<-function(y,x,r,pi,id,time,beta,rho,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    rho<-roo(rho,time,corstr)
    A<-diag(1,time,time)
    R<-R(rho,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    # z.col<-ncol(z)
    # rlag<-ylag(id,r,1)
    # rlag[is.na(rlag)]<-1
    # S<-matrix(rep(0),nrow = z.col, ncol = n )


    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        piii<-marg.prob(pi[((i-1)*time+1):(i*time)])
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        WW<-W*r[((i-1)*time+1):(i*time)]*piii
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        # e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        # for (m in 1:(time-1)) #formula for rho
        # {error<-0
        # for (j in 1:(time-m))
        # {
        #     error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        # }
        # error<-error-rho[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        # wgeei<-rbind(wgeei,error)
        # }

        wgeef<-cbind(wgeef,wgeei)
        # ss<-t(rlag[id==i]*(r[id==i]-lambda[id==i])*z[id==i,])
        # S[,i]<-as.matrix(apply(ss,1,sum))
    }
    return(wgeef)
}

