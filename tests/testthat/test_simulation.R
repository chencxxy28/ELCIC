
test_that("output error: invalid dist type in glm",{expect_error(glm.generator(beta=c(0.5,0.5,0.5,0),samplesize=100,rho=0.5,dist="gamma",ov=2),"Invalid type of dist. It should be one of gaussian,binomial,poisson,NB")})

test_that("output error: invalid dist type in gee without missing",{expect_error(gee.generator(beta=c(-1,1,0.5,0),samplesize=200,time=3,num.time.dep=2,num.time.indep=1,rho=0.4,x.rho=0.2,dist="gamma",cor.str="exchangeable",x.cor.str="exchangeable"),"Invalid type of dist. It should be one of gaussian,binomial,poisson")})

test_that("output error: invalid correlation type for outcomes in gee without missing",{expect_error(gee.generator(beta=c(-1,1,0.5,0),samplesize=200,time=3,num.time.dep=2,num.time.indep=1,rho=0.4,x.rho=0.2,dist="gaussian",cor.str="ar2",x.cor.str="exchangeable"),"Invalid type of correlation.struture for outcomes. It should be one of exchangeable, ar1, and independence")})

test_that("output error: invalid correlation type for outcomes in gee without missing",{expect_error(gee.generator(beta=c(-1,1,0.5,0),samplesize=200,time=3,num.time.dep=2,num.time.indep=1,rho=0.4,x.rho=0.2,dist="gaussian",cor.str="ar1",x.cor.str="ar2"),"Invalid type of correlation.struture for covariates. It should be one of exchangeable, ar1, and independence")})

test_that("output error: invalid input for num.time.dep and num.time.indep",{expect_error(gee.generator(beta=c(-1,1,0.5,0),samplesize=200,time=3,num.time.dep=5,num.time.indep=1,rho=0.4,x.rho=0.2,dist="gaussian",cor.str="ar1",x.cor.str="ar1"),"Invalid input for number of time independent covariate and dependent covariate. Sum should be equal to the length of beta minus one")})

