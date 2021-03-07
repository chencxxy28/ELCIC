##first derivative of -log EL
R1der<-function(lambda,ZZ)
{
    apply(ZZ,2,function(xx)
    {as.matrix(xx,ncol=1)/as.vector((1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}

##second derivative of -log EL
R2der<-function(lambda,ZZ)
{
    r2der<-0
    for(i in 1:ncol(ZZ))
    {
        r2der_i<--as.matrix(ZZ[,i],ncol=1)%*%t(as.matrix(ZZ[,i],ncol=1))/as.vector(1+t(lambda)%*%as.matrix(ZZ[,i],ncol=1))^2
        r2der<-r2der+r2der_i
    }
    r2der
}

##-log EL
R0der<-function(lambda,ZZ)
{
    apply(ZZ,2, function (xx) {log(as.vector(1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}

#function to find lambda, given beta and an
lambda.find.glm<-function(x,y,betahat,dist)
{
    samplesize<-nrow(x)
    ZZ<-ee.glm(x=x,y=y,betahat=betahat,dist=dist)[,1:(samplesize)]
    #ZZ<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-8

    repeat{
        rl<-R1der(lambda,ZZ)
        rll<-R2der(lambda,ZZ)
        Delta<--ginv(rll)%*%rl
        if(mean(abs(Delta))<tol)
        {break}else{
            repeat{
                mm<-0
                repeat{
                    delta<-gamma*Delta
                    index_1<-apply(ZZ,2,function (xx)
                    {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
                    )
                    if (sum(index_1)>0)
                    {gamma<-gamma/2
                    mm<-mm+1}else{break}}
                index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
                if (index_2==1)
                {gamma<-gamma/2}else{break}
            }
        }
        lambda<-lambda+delta
        c<-c+1
        gamma<-(c)^(-0.5)
    }
    lambda
}
