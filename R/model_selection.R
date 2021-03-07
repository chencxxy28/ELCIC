ELCIC.glm<-function (x,y,index.var,name.var=NULL,dist)
{
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
