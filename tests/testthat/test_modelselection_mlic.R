data("wgeesimdata")


output1<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="ar1",joints=FALSE)
output2<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets=c("ar1","exchangeable"),joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})


output1<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-MLIC.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})
