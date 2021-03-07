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
