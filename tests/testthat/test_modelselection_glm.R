data("glmtoydata")

test_that("output error: y is not a vector",{expect_error(ELCIC.glm(x=glmtoydata$x,y=as.data.frame(glmtoydata$y),index.var=c(1,2,3),name.var=NULL,
                                              dist = "gaussian"),"y should be in a vector format")})

test_that("output error: x is not a matrix",{expect_error(ELCIC.glm(x=data.frame(glmtoydata$x),y=as.data.frame(glmtoydata$y),index.var=c(1,2,3),name.var=NULL,
                                                                    dist = "gaussian"),"x should be in a matrix format")})

test_that("output error: dist is not undefined",{expect_error(ELCIC.glm(x=glmtoydata$x,y=glmtoydata$y,index.var=c(1,2,3),name.var=NULL,
                                                 dist = "gamma"),"Invalid type of dist. It should be one of gaussian,binomial,poisson")})


test_that("output error: non-unique index",{expect_error(ELCIC.glm(x=glmtoydata$x,y=glmtoydata$y,index.var=c(1,1,3),name.var=NULL,
                                                 dist = "gaussian"),"Invalid candidate model provided")})

test_that("output error: undefined candidate model",{expect_error(ELCIC.glm(x=glmtoydata$x,y=glmtoydata$y,index.var=NULL,name.var=c("intercept","x1","x2","x3","x4"),
                                                 dist = "gaussian"),"Invalid candidate model provided")})

test_that("output error: non-unique variable name",{expect_error(ELCIC.glm(x=glmtoydata$x,y=glmtoydata$y,index.var=NULL,name.var=c("intercept","x1","x1"),
                                                 dist = "gaussian"),"Invalid candidate model provided")})



# name.var.set <- list(c("intercept","x1"),c("intercept","x1","x2"))
# transform.index<-list(c(1,2),c(1:3))
# test_that("var.name transforms to index",{expect_equal(match.index(x=glmtoydata$x,name.var.set),transform.index)})
