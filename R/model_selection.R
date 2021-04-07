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


ELCIC.gee<-function(x,y,r,id,time,index.var,name.var=NULL,dist,corstr)
{
    samplesize<-nrow(x)/time
    betahat<-rep(0,ncol(x))
    x_candidate<-(x[,index.var])
    colnames(x_candidate)<-rep(1:ncol(x_candidate))
    y<-as.matrix(y,col=1)
    data<-data.frame(y,x_candidate)

        fit<-geeglm(y~x_candidate-1,data=data,family =dist,id=id,corstr = corstr)
        betahat[index.var]<-fit$coefficients
        p<-length(fit$coefficients)

        if(corstr=="independence")
        {
            ro=0
        }else{
            ro<-unlist(summary(fit)$corr[1]) #correlation coefficients
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
}

ELCIC.glm.procedure<-function(x,y,candidate.sets,name.var=NULL,dist)
{
    criterion.all<-rep()
    for (i in 1:length(candidate.sets))
    {
        criterion<-ELCIC.glm(x=x,y=y,index.var=candidate.sets[[i]],name.var=NULL,dist=dist)
        criterion.all<-cbind(criterion.all,criterion)
    }
    criterion.all
}


ELCIC.gee.procedure<-function(x,y,r,id,time,candidate.sets,name.var=NULL,dist,candidate.cor.sets)
{
    criterion.elcic<-rep()
    for(j in 1:length(candidate.cor.sets))
    {
        corstr<-candidate.cor.sets[j]
        criterion.all<-rep()
        for (i in 1:length(candidate.sets))
        {
            criterion<-ELCIC.gee(x=x,y=y,r=r,id=id,time=time,index.var=candidate.sets[[i]],name.var=NULL,dist="poisson",corstr=corstr)
            criterion.all<-c(criterion.all,criterion)
        }
        criterion.elcic<-rbind(criterion.elcic,criterion.all)
        #print(j)
    }
    colnames(criterion.elcic)<-candidate.sets
    rownames(criterion.elcic)<-candidate.cor.sets
    criterion.elcic
}



QIC.procedure<-function (x,y,id,dist,candidate.sets,name.var=NULL,candidate.cor.sets)
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
