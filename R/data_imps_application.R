setwd("/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/phd_in_psu/research/Prof_Wang/MEL/coding/newest/reviewer_code_09_02/upload codes/real data/table6/")
# rm(list=ls(all=TRUE))
# library(rootSolve)
# library(PoisNor)
# library(bindata)
# #library(geepack)
# #library(glmnet)
# library(MASS)
# library(regpro)
# #library(SparseM)
# #library(quantreg)
# #library(rqPen)
# library(denpro)
# library(regpro)
library(wgeesel)


# v<-function(x,beta,dist) #function used to calculate variance, mean, and derivative of random variables followed by binomial, poisson, and gaussian
# {
#     n<-length(x[,1])
#     reg<-x%*%beta
#     pi<-exp(reg)/(1+exp(reg))
#     if (dist=="gaussian")
#     {
#         v<-rep(1,n)
#         der<-x
#         mu<-reg
#     }
#     else if (dist=="binomial")
#     {
#         v<-pi*(1-pi)
#         der<-x*as.vector(v)
#         mu<-as.vector(pi)
#     }
#     else if (dist=="poisson")
#     {
#         v<-exp(reg)
#         der<-x*as.vector(exp(reg))
#         mu<-as.vector(exp(reg))
#     }
#     list(v=v,der=der,mu=mu)
# }
#
#
# R<-function(ro,id) #working correlation matrix
# {
#     time<-length(id)/length(unique(id))
#     R<-diag(1,time,time)
#     ro<-c(1,ro)
#     for (i in 1:time)
#     {
#         for (j in 1:time)
#         {
#             R[i,j]<-ro[abs(i-j)+1]
#         }
#     }
#     return(R)
# }
#
#
# pii<-function(pi) #missing prob
# { pi[which(pi<0.00001)]<-0.00001
# pii<-1/pi
# time<-length(pi)
# for (i in 2:time)
# {
#     pii[i]<-pii[i-1]*pii[i]
# }
# return(pii)
# }
#
#
#
# wgeef<-function(y,x,xf,r,lambda,z,pi,id,beta,betaf,ro,phi,dist,corstr)
# { #full wgee
#     n<-length(unique(id))
#     p<-length(beta)
#     pf<-length(betaf)
#     time<-length(id)/n
#     A<-diag(1,time,time)
#     R<-R(ro,id)
#     W<-diag(1,time,time)
#     V<-v(x,beta,dist)
#     wgeef<-rep()
#
#     z.col<-ncol(z)
#     rlag<-ylag(id,r,1)
#     rlag[is.na(rlag)]<-1
#     S<-matrix(rep(0),nrow = z.col, ncol = n )
#
#
#     y[which(is.na(y))]<-0
#     for (i in 1:n)
#     {
#         piii<-pii(pi[((i-1)*time+1):(i*time)])
#         AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
#         WW<-W*r[((i-1)*time+1):(i*time)]*piii
#         e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
#         wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
#         e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
#         for (m in 1:(time-1)) #formula for ro
#         {error<-0
#         for (j in 1:(time-m))
#         {
#             error<-error+e[j]*e[j+m]*WW[j+m,j+m]
#         }
#         error<-error-ro[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
#         wgeei<-rbind(wgeei,error)
#         }
#
#         wgeef<-cbind(wgeef,wgeei)
#
#         ss<-t(rlag[id==i]*(r[id==i]-lambda[id==i])*z[id==i,])
#         S[,i]<-as.matrix(apply(ss,1,sum))
#     }
#
#     wgeef<-rbind(wgeef,S)
#
#     return(wgeef)
# }
#
#
#
#
# roo<-function(ro,time,corstr)
# {
#     if (ro>0.99)
#     {ro<-0.99}
#     roo<-rep(0,(time-1))
#     if (corstr=="ar1")
#     {
#         for (i in 1:(time-1))
#         {roo[i]<-ro^i}
#     }
#     else if (corstr=="exchangeable")
#     {
#         roo<-rep(ro,(time-1))
#     }
#     return(roo)
# }




######################################
time<-4
# countaic<-matrix(0,nrow=3,ncol=6)
# countbic<-matrix(0,nrow=3,ncol=6)
# structure<-c("ar1","exchangeable","independence")

# criteriaaic<-rep()
# criteriabic<-rep()
#correlation<-c("ar1","exchangeable","independence")
dist<-"binomial" #"binomial","poisson"
# WEAIC<-matrix(0,nrow=3,ncol=6)
# WEBIC<-matrix(0,nrow=3,ncol=6)

load("imps.rdata")
imps$subject <- imps$ID
id<-rep(1:(length(imps$subject)/time),each=time)
# lag1y <- ylag(imps$ID,imps$Y,1) ###create lagged y(t-1)##
# lag2y <- ylag(imps$ID,imps$Y,2,na=F) ###create lagged y(t-2)I(t>2)##
# lag3y <- ylag(imps$ID,imps$Y,3,na=F) ###create lagged y(t-3)I(t>3)##
# imps_new <- cbind(imps,lag1y,lag2y,lag3y)
imps_new$Week<-sqrt(imps_new$Week)
y<-imps_new$Y
r<-imps_new$R
# n<-length(id)/time
# epi<-1/n
# mis<-cbind(1,imps_new$Drug,imps_new$Time,imps_new$Sex,imps_new$lag1y,imps_new$lag2y,imps_new$lag3y)
x_mis<-cbind(Intercept=1,drug=imps_new$Drug,time=imps_new$Time,sex=imps_new$Sex)
x<-cbind(1,imps_new$Time,imps_new$Sex,imps_new$Drug,imps_new$Time*imps_new$Sex,
          imps_new$Sex*imps_new$Drug,imps_new$Drug*imps_new$Time)
colnames(x)<-c("Intercept","Time","Sex","Drug","Time*Sex","Sex*Drug","Drug*Time")
# zz<-z
# zz[is.na(zz)]<-1
impsdata<-list(y=y,x=x,r=r,id=id,x_mis=x_mis)
save(impsdata, file = "/Users/chixiang.chen/OneDrive - University of Maryland School of Medicine/packages/ELCIC/data/impsdata.RData")


candidate.sets<-list(c(1:2),c(1,4),c(1:2,4),c(1:2,4,7),c(1:4),c(1:7))
candidate.cor.sets<-c("exchangeable","ar1","independence")
criterion.elcic<-ELCIC.wgee(x,y,x_mis,r,id,time,candidate.sets,name.var.sets=NULL,
                                     dist,candidate.cor.sets,joints=TRUE,lag=2)
criterion.mlic<-MLIC.wgee(x,y,x_mis,r,id,time,candidate.sets,
                                     name.var.sets=NULL,dist,candidate.cor.sets,joints=TRUE,lag=2)
criterion.qicw<-QICW.wgee(x,y,x_mis,r,id,time,candidate.sets,
                                     name.var.sets=NULL,dist,candidate.cor.sets,joints=TRUE,lag=2)

criterion.elcic
criterion.mlic
criterion.qicw









#
#
# ########################
# corstr<-"ar1"
# ########################
#
# #time
# x<-xf[,c(1,2)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1,2)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[1,1]<-likelihood+2*(p+1)
# WEBIC[1,1]<-likelihood+(p+1)*log(n)
#
#
# #drug
# x<-xf[,c(1,4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1,4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[1,2]<-likelihood+2*(p+1)
# WEBIC[1,2]<-likelihood+(p+1)*log(n)
#
#
# #time drug
# x<-xf[,c(1:2,4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:2,4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[1,3]<-likelihood+2*(p+1)
# WEBIC[1,3]<-likelihood+(p+1)*log(n)
#
# #time drug time*drug
# x<-xf[,c(1:2,4,7)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:2,4,7)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[1,4]<-likelihood+2*(p+1)
# WEBIC[1,4]<-likelihood+(p+1)*log(n)
#
# #time sex drug
# x<-xf[,c(1:4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[1,5]<-likelihood+2*(p+1)
# WEBIC[1,5]<-likelihood+(p+1)*log(n)
#
# #time,sex,drug, time*sex,sex*drug, drug*time
# x<-xf
# colnames(xf)<-c("Intercept","Time","Sex","Drug","Time*Sex","Sex*Drug","Drug*Time")
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[1:7]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[1,6]<-likelihood+2*(p+1)
# WEBIC[1,6]<-likelihood+(p+1)*log(n)
#
#
# ########################
# corstr<-"exchangeable"
# ########################
#
# #time
# x<-xf[,c(1,2)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1,2)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[2,1]<-likelihood+2*(p+1)
# WEBIC[2,1]<-likelihood+(p+1)*log(n)
#
#
# #drug
# x<-xf[,c(1,4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1,4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[2,2]<-likelihood+2*(p+1)
# WEBIC[2,2]<-likelihood+(p+1)*log(n)
#
#
# #time drug
# x<-xf[,c(1:2,4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:2,4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[2,3]<-likelihood+2*(p+1)
# WEBIC[2,3]<-likelihood+(p+1)*log(n)
#
# #time drug time*drug
# x<-xf[,c(1:2,4,7)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:2,4,7)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[2,4]<-likelihood+2*(p+1)
# WEBIC[2,4]<-likelihood+(p+1)*log(n)
#
# #time sex drug
# x<-xf[,c(1:4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[2,5]<-likelihood+2*(p+1)
# WEBIC[2,5]<-likelihood+(p+1)*log(n)
#
# #time,sex,drug, time*sex,sex*drug, drug*time
# x<-xf
# colnames(xf)<-c("Intercept","Time","Sex","Drug","Time*Sex","Sex*Drug","Drug*Time")
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[1:7]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[2,6]<-likelihood+2*(p+1)
# WEBIC[2,6]<-likelihood+(p+1)*log(n)
#
# ########################
# corstr<-"independence"
# ########################
#
# #time
# x<-xf[,c(1,2)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1,2)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[3,1]<-likelihood+2*(p)
# WEBIC[3,1]<-likelihood+(p)*log(n)
#
#
# #drug
# x<-xf[,c(1,4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1,4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[3,2]<-likelihood+2*(p)
# WEBIC[3,2]<-likelihood+(p)*log(n)
#
#
# #time drug
# x<-xf[,c(1:2,4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:2,4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[3,3]<-likelihood+2*(p)
# WEBIC[3,3]<-likelihood+(p)*log(n)
#
# #time drug time*drug
# x<-xf[,c(1:2,4,7)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:2,4,7)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[3,4]<-likelihood+2*(p)
# WEBIC[3,4]<-likelihood+(p)*log(n)
#
# #time sex drug
# x<-xf[,c(1:4)]
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[c(1:4)]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[3,5]<-likelihood+2*(p)
# WEBIC[3,5]<-likelihood+(p)*log(n)
#
# #time,sex,drug, time*sex,sex*drug, drug*time
# x<-xf
# colnames(xf)<-c("Intercept","Time","Sex","Drug","Time*Sex","Sex*Drug","Drug*Time")
# fit <- wgee(Y~x-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# beta<-as.vector(summary(fit)$beta)  #c(1.0624,0.9779)
# p<-length(beta)
# gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1] )    #c(0.7929,1.0142,1.0698)
# ro<-summary(fit)$corr
# ro<-roo(ro,time,corstr)
# phihat<-summary(fit)$phi
# pi<-exp(z%*%gamma)/(1+exp(z%*%gamma))
# pi[which(is.na(pi))]<-1
# betahat<-rep(0,7)
# betahat[1:7]<-beta
#
# Z<-as.matrix(wgeef(y,xf,xf,r,pi,zz,pi,id,betahat,betahat,ro,phihat,dist,corstr))
# #Z%*%rep(1,n)
#
#
# model<-function(lambda)
# {
#     apply(Z,2,function(x) if (1+t(lambda)%*%x>=epi)
#     {x/(1+t(lambda)%*%x)}
#     else {2/epi*x-(1+t(lambda)%*%x)/(epi^2)*x})%*%rep(1,n)
# }
#
# lambda<-multiroot(f = model, start = rep(0,length(betahat)+time-1+ncol(zz)))$root
#
# likelihood<-apply(Z,2,function(x) {
#
#     2*log(1+t(lambda)%*%x)})%*%rep(1,n)
#
# WEAIC[3,6]<-likelihood+2*(p)
# WEBIC[3,6]<-likelihood+(p)*log(n)
#
# ########################
# #MLIC and QICW
# ########################
# #########################################begin correlation selection based on full model
# corstr<-"ar1"
#
# fit <- wgee(Y~xf-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
#
#
# MLIC1<-MLIC.gee(fit,fit)$MLICc
# QICW1<-QICW.gee(fit)$QICWr
#
#
#
# #########################################
# corstr<-"exchangeable"
#
# fit <- wgee(Y~xf-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
# MLIC2<-MLIC.gee(fit,fit)$MLICc
# QICW2<-QICW.gee(fit)$QICWr
#
#
#
#
# #########################################
# corstr<-"independence"
#
# fit <- wgee(Y~xf-1,imps_new,
#             imps_new$ID,family=dist,corstr =corstr,scale = NULL,
#             mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
#
# MLIC3<-MLIC.gee(fit,fit)$MLICc
# QICW3<-QICW.gee(fit)$QICWr
#
# corindexm<-which.min(c(MLIC1,MLIC2,MLIC3))
# corstrm<-structure[corindexm]
#
# corindexq<-which.min(c(QICW1,QICW2,QICW3))
# corstrq<-structure[corindexq]
#
# #########now to do model selection
#
#
# fitf <- wgee(Y~xf-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrm,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
# #1:time
# #########################################
# x<-xf[,1:2]
#
# fitm <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrm,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# fitq <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrq,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
# MLIC1<-MLIC.gee(fitm,fitf)$MLIC
# QICW1<-QICW.gee(fitq)$QICWr
#
# #2:drug
# x<-xf[,c(1,4)]
# fitm <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrm,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# fitq <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrq,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
# MLIC2<-MLIC.gee(fitm,fitf)$MLIC
# QICW2<-QICW.gee(fitq)$QICWr
#
# #3:time drug
# x<-xf[,c(1:2,4)]
# fitm <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrm,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# fitq <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrq,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
# MLIC3<-MLIC.gee(fitm,fitf)$MLIC
# QICW3<-QICW.gee(fitq)$QICWr
#
# #4:time drug time*drug
# x<-xf[,c(1:2,4,7)]
# fitm <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrm,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# fitq <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrq,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
# MLIC4<-MLIC.gee(fitm,fitf)$MLIC
# QICW4<-QICW.gee(fitq)$QICWr
#
# #5:time sex drug
# x<-xf[,c(1:4)]
# fitm <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrm,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# fitq <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrq,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
# MLIC5<-MLIC.gee(fitm,fitf)$MLIC
# QICW5<-QICW.gee(fitq)$QICWr
#
# #6:time,sex,drug, time*sex,sex*drug, drug*time
# x<-xf
# fitm <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrm,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
# fitq <- wgee(Y~x-1,imps_new,
#              imps_new$ID,family=dist,corstr =corstrq,scale = NULL,
#              mismodel =R~Drug+Time+Sex+lag1y+lag2y+lag3y)
#
#
# MLIC6<-MLIC.gee(fitm,fitf)$MLIC
# QICW6<-QICW.gee(fitq)$QICWr
#
# which(WEAIC == min(WEAIC), arr.ind = TRUE)
# which(WEBIC == min(WEBIC), arr.ind = TRUE)
# which.min(c(MLIC1,MLIC2,MLIC3,MLIC4,MLIC5,MLIC6))
# which.min(c(QICW1,QICW2,QICW3,QICW4,QICW5,QICW6))
#
# WEAIC
# WEBIC
# c(MLIC1,MLIC2,MLIC3,MLIC4,MLIC5,MLIC6)
# c(QICW1,QICW2,QICW3,QICW4,QICW5,QICW6)



