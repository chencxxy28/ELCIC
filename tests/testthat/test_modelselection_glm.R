data("glmsimdata")

test_that("output error: y is not a vector",{expect_error(ELCIC.glm.single(x=glmsimdata$x,y=as.data.frame(glmsimdata$y),index.var=c(1,2,3),name.var=NULL,
                                              dist = "gaussian"),"y should be in a vector format")})

test_that("output error: x is not a matrix",{expect_error(ELCIC.glm.single(x=data.frame(glmsimdata$x),y=as.data.frame(glmsimdata$y),index.var=c(1,2,3),name.var=NULL,
                                                                    dist = "gaussian"),"x should be in a matrix format")})

test_that("output error: dist is not undefined",{expect_error(ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=c(1,2,3),name.var=NULL,
                                                 dist = "gamma"),"Invalid type of dist. It should be one of gaussian,binomial,poisson")})


test_that("output error: non-unique index",{expect_error(ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=c(1,1,3),name.var=NULL,
                                                 dist = "gaussian"),"Invalid candidate model provided")})

test_that("output error: undefined candidate model",{expect_error(ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=NULL,name.var=c("intercept","x1","x2","x3","x4"),
                                                 dist = "gaussian"),"Invalid candidate model provided")})

test_that("output error: non-unique variable name",{expect_error(ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=NULL,name.var=c("intercept","x1","x1"),
                                                 dist = "gaussian"),"Invalid candidate model provided")})

output1<-ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=NULL,name.var=c("intercept","x1","x2"),
          dist = "gaussian")
output2<-ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),
                   dist = "gaussian")
test_that("output equal: same output given both index and var.names",{expect_equal(output1,output2)})


output1<-ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=NULL,name.var=c("intercept","x1","x2"),
                   dist = "poisson")
output2<-ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),
                   dist = "poisson")
test_that("output equal: same output given both index and var.names",{expect_equal(output1,output2)})


output1<-ELCIC.glm(x=glmsimdata$x,y=glmsimdata$y,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),
                   dist = "poisson")
output2<-ELCIC.glm(x=glmsimdata$x,y=glmsimdata$y,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),
                   dist = "poisson")
test_that("output equal: same output given both index and var.names",{expect_equal(output1,output2)})


output1<-as.vector(ELCIC.glm(x=glmsimdata$x,y=glmsimdata$y,name.var.sets=list(c("intercept","x1","x2")),
                             dist = "poisson"))
output2<-as.vector(ELCIC.glm(x=glmsimdata$x,y=glmsimdata$y,candidate.sets=list(c(1,2,3)),
                             dist = "poisson"))
test_that("output equal: same output given both index and var.names",{expect_equal(output1,output2)})


#generate binomial
set.seed(515413)
samplesize<-100
glmsimdata<-glm.generator(beta=c(0.5,0.5,0.5,0),samplesize=samplesize,rho=0.5,dist="binomial")
rownames(glmsimdata$x)<-seq_len(nrow(glmsimdata$x))
colnames(glmsimdata$x)<-c("intercept","x1","x2","x3")

output1<-ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=NULL,name.var=c("intercept","x1","x2"),
                   dist = "binomial")
output2<-ELCIC.glm.single(x=glmsimdata$x,y=glmsimdata$y,index.var=c(1,2,3),name.var=c("intercept","x1","x2"),
                   dist = "binomial")
test_that("output equal: same output given both index and var.names",{expect_equal(output1,output2)})







# name.var.set <- list(c("intercept","x1"),c("intercept","x1","x2"))
# transform.index<-list(c(1,2),c(1:3))
# test_that("var.name transforms to index",{expect_equal(match.index(x=glmsimdata$x,name.var.set),transform.index)})
