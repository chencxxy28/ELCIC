data("wgeesimdata")

test_that("output error: y is not a vector",{expect_error(ELCIC.wgee.single(x=wgeesimdata$x,y=data.frame(wgeesimdata$y)
                                                                     ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"y should be in a vector format")})

test_that("output error: x is not a matrix",{expect_error(ELCIC.wgee.single(x=data.frame(wgeesimdata$x),y=(wgeesimdata$y)
                                                                     ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"x should be in a matrix format")})

test_that("output error: x_mis is not a matrix",{expect_error(ELCIC.wgee.single(x=(wgeesimdata$x),y=(wgeesimdata$y)
                                                                     ,x_mis=data.frame(wgeesimdata$x_mis),r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"x_mis should be in a matrix format")})

test_that("output error: dist is not undefined",{expect_error(ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                                                                         ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="gamma",corstr="exchangeable",joints=T),"Invalid type of dist. It should be one of gaussian,binomial,poisson")})


test_that("output error: non-unique index",{expect_error(ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                                                                    ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,1,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"Invalid candidate model provided")})

test_that("output error: undefined candidate model",{expect_error(ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                                                                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2","x3","x4"),dist="binomial",corstr="exchangeable",joints=T),"Invalid candidate model provided")})

test_that("output error: non-unique variable name",{expect_error(ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                                                                            ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x1"),dist="binomial",corstr="exchangeable",joints=T),"Invalid candidate model provided")})

test_that("output error: invalid correlation structure for outcomes in gee without missing",{expect_error(ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                                                                                                                     ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="binomial",corstr="ar2",joints=T),"Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")})

test_that("output error: incorrect lag input",{expect_error(ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                                                                                ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="poisson",corstr="exchangeable",joints=T,lag=3),"Invalid type of lag. It should be less than 3 and time")})

test_that("output error: incorrect lag input",{expect_warning(ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                                                                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="poisson",corstr="exchangeable",joints=T,lag=0),"No lag of response added may indicate missing completely at random. GEE may be used")})


output1<-ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                    ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=T)
output2<-ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                    ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=T)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                    ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=FALSE)
output2<-ELCIC.wgee.single(x=wgeesimdata$x,y=(wgeesimdata$y)
                    ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                    ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                    ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets=c("exchangeable","ar1"),joints=FALSE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets=c("independence","ar1"),joints=FALSE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="independence",joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})

#try different distributions
set.seed(515413)
library(wgeesel)
n=200
time=3
rho=0.5
id=rep(1:n,each=time)
truecorstr<-"exchangeable" #"ar1","exchangeable"
truedist<-"poisson" #"binary","poisson"
xf <- cbind(1, rep(runif(n), each=time), rep(0:(time-1),times=n), rnorm(time*n))
xT<-xf[,1:3]
betaT=c(-1,1,0.4)
x_mis=cbind(matrix(1,nrow=length(id),ncol=1),runif(n*time,-0.5,0.5))
para=c(1.74,0.5,-0.8)
data_gernated<-data_sim(id,rho,phi,xT,betaT,x_mis,para,truecorstr,truedist,lag_level = 1)

x<-xf
colnames(x)<-c("intercept","x1","x2","x3")
x_mis<-as.matrix(data.frame(intercept=data_gernated$data$V1,x_mis1=data_gernated$data$V2))
wgeesimdata<-list(y=data_gernated$data$response_mis,x=x,x_mis=x_mis,id=id,obs_ind=data_gernated$data$ind)


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets=c("exchangeable","ar1"),joints=FALSE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="gaussian",candidate.cor.sets=c("exchangeable","ar1"),joints=FALSE)
output2<-ELCIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="gaussian",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})

