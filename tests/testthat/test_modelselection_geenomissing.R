data("geetoydata")

x<-geetoydata$x

test_that("output error: y is not a vector",{expect_error(ELCIC.gee(x=geetoydata$x,y=as.data.frame(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="poisson",corstr="exchangeable"),"y should be in a vector format")})

test_that("output error: x is not a matrix",{expect_error(ELCIC.gee(x=as.data.frame(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="poisson",corstr="exchangeable"),"x should be in a matrix format")})

test_that("output error: dist is not undefined",{expect_error(ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=c(1,2,3),name.var=NULL,dist="gamma",corstr="exchangeable"),"Invalid type of dist. It should be one of gaussian,binomial,poisson")})


test_that("output error: non-unique index",{expect_error(ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=c(1,1,3),name.var=NULL,dist="poisson",corstr="exchangeable"),"Invalid candidate model provided")})

test_that("output error: undefined candidate model",{expect_error(ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2","x3","x4"),dist="poisson",corstr="exchangeable"),"Invalid candidate model provided")})

test_that("output error: non-unique variable name",{expect_error(ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x1"),dist="poisson",corstr="exchangeable"),"Invalid candidate model provided")})

test_that("output error: invalid correlation structure for outcomes in gee without missing",{expect_error(ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar2"),"Invalid type of correlation structure for outcomes. It should be one of ar1,exchangeable,independence")})

output1<-ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=TRUE)
output2<-ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=TRUE)
test_that("output equal: equal values given both index and name, when joint=true",{expect_equal(output1,output2)})

output1<-ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=NULL,name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=FALSE)
output2<-ELCIC.gee(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),dist="poisson",corstr="ar1",joint=FALSE)
test_that("output equal: equal values given both index and name, when joint=false",{expect_equal(output1,output2)})

output1<-ELCIC.gee.procedure(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=TRUE)
output2<-ELCIC.gee.procedure(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=TRUE)
test_that("output equal: equal values given both index and name, when joint=true",{expect_equal(output1,output2)})

output1<-ELCIC.gee.procedure(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
output2<-ELCIC.gee.procedure(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
test_that("output equal: equal values given both index and name, when joint=false",{expect_equal(output1,output2)})

output1<-ELCIC.gee.procedure(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets="ar1",joint=FALSE)
output2<-ELCIC.gee.procedure(x=(geetoydata$x),y=(geetoydata$y),r=rep(1,nrow(x)),id=geetoydata$id,time=3,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="poisson",candidate.cor.sets=c("ar1","exchangeable"),joint=FALSE)
test_that("output equal: equal values given more than 1 correlation structures, when joint=false",{expect_equal(output1,output2)})

