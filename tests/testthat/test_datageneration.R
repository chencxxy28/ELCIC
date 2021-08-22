
samplesize<-100
glmtoydata1<-glm.generator(beta=c(0.5,0.5,0.5,0),samplesize=samplesize,rho=0.5,dist="NB",ov=2)
test_that("output equal: same length",{expect_equal(length(glmtoydata1),2)})


samplesize<-200
time=3
geetoydata<-gee.generator(beta=c(-1,1,0.5,0),samplesize=samplesize,time=time,num.time.dep=2,num.time.indep=1,rho=0.4,x.rho=0.2,dist="gaussian",cor.str="exchangeable",x.cor.str="exchangeable")
test_that("output equal: same length",{expect_equal(length(geetoydata),3)})

