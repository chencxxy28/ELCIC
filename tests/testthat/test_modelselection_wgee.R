data("wgeetoydata")

test_that("output error: y is not a vector",{expect_error(ELCIC.wgee(x=wgeetoydata$x,y=data.frame(wgeetoydata$y)
                                                                     ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"y should be in a vector format")})

test_that("output error: x is not a matrix",{expect_error(ELCIC.wgee(x=data.frame(wgeetoydata$x),y=(wgeetoydata$y)
                                                                     ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"x should be in a matrix format")})

test_that("output error: x_mis is not a matrix",{expect_error(ELCIC.wgee(x=(wgeetoydata$x),y=(wgeetoydata$y)
                                                                     ,x_mis=data.frame(wgeetoydata$x_mis),r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"x_mis should be in a matrix format")})

test_that("output error: dist is not undefined",{expect_error(ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                                                                         ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="gamma",corstr="exchangeable",joints=T),"Invalid type of dist. It should be one of gaussian,binomial,poisson")})


test_that("output error: non-unique index",{expect_error(ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                                                                    ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=c(1,1,3),name.var=NULL,dist="binomial",corstr="exchangeable",joints=T),"Invalid candidate model provided")})

test_that("output error: undefined candidate model",{expect_error(ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                                                                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2","x3","x4"),dist="binomial",corstr="exchangeable",joints=T),"Invalid candidate model provided")})

test_that("output error: non-unique variable name",{expect_error(ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                                                                            ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x1"),dist="binomial",corstr="exchangeable",joints=T),"Invalid candidate model provided")})

test_that("output error: invalid correlation structure for outcomes in gee without missing",{expect_error(ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                                                                                                                     ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="binomial",corstr="ar2",joints=T),"Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")})

output1<-ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                    ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=T)
output2<-ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                    ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=T)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                    ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=FALSE)
output2<-ELCIC.wgee(x=wgeetoydata$x,y=(wgeetoydata$y)
                    ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="binomial",corstr="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-ELCIC.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                    ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-ELCIC.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                    ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-ELCIC.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                              ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-ELCIC.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                              ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})

output1<-ELCIC.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                              ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets=c("exchangeable","ar1"),joints=FALSE)
output2<-ELCIC.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                              ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})



