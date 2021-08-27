data("wgeesimdata")


output1<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                              ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})



output1<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets=c("exchangeable","ar1"),joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})


output1<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-QICW.wgee(x=wgeesimdata$x,y=(wgeesimdata$y)
                             ,x_mis=wgeesimdata$x_mis,r=wgeesimdata$obs_ind,id=wgeesimdata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=FALSE",{expect_equal(output1,output2)})
