data("wgeetoydata")


output1<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                              ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                              ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})



output1<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
output2<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=TRUE)
test_that("output equal: same output given both index and var.names, given joints=true",{expect_equal(output1,output2)})


output1<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=list(c(1,2,3)),name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=false",{expect_equal(output1,output2)})


output1<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=NULL,name.var.sets=list(c("intercept","x1","x2")),dist="binomial",candidate.cor.sets=c("exchangeable","ar1"),joints=FALSE)
test_that("output equal: same output given multiple correlation structures, given joints=false",{expect_equal(output1,output2)})


output1<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
output2<-QICW.wgee.procedure(x=wgeetoydata$x,y=(wgeetoydata$y)
                             ,x_mis=wgeetoydata$x_mis,r=wgeetoydata$obs_ind,id=wgeetoydata$id,candidate.sets=list(c(1,2,3)),name.var.sets=NULL,dist="binomial",candidate.cor.sets="exchangeable",joints=FALSE)
test_that("output equal: same output given both index and var.names, given joints=FALSE",{expect_equal(output1,output2)})
