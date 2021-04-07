######for cross-sectional glm
ee.glm<-function (x,y,betahat,dist)
{
    full=ncol(x)
    samplesize<-nrow(x)
    ee<-matrix(0,nrow=full,ncol=samplesize)
    if(dist=="gaussian")
    {
    for (i in 1:samplesize)
    {
        ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(t(as.matrix(x[i,],ncol=1))%*%betahat))
    }
    }else if (dist=="poisson")
    {
        for (i in 1:samplesize)
        {
        ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-exp(t(as.matrix(x[i,],ncol=1))%*%betahat))
        }
    }else if (dist=="binomial")
    {
        for (i in 1:samplesize)
        {
        ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(1+exp(-t(as.matrix(x[i,],ncol=1))%*%betahat))^(-1))
        }
    }
    ee
}


######for longitudinal data gee
ee.gee<-function(y,x,r,id,beta,ro,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    time<-length(id)/n
    A<-diag(1,time,time)
    ro<-roo(ro,time,corstr)
    R<-R(ro,id)
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
        e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        for (m in 1:(time-1)) #formula for ro
        {error<-0
        for (j in 1:(time-m))
        {
            error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        }
        error<-error-ro[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        wgeei<-rbind(wgeei,error)
        }

        wgeef<-cbind(wgeef,wgeei)
    }

    return(wgeef)
}

