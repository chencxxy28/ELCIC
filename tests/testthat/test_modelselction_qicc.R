data(geetoydata)
x<-geetoydata$x
y<-geetoydata$y
id<-geetoydata$id
r<-rep(1,nrow(x))
time<-3
candidate.sets<-list(c(1,2))
name.var.sets<-list(c("intercept","x1"))
candidate.cor.sets<-c("exchangeable")
dist="poisson"
output1<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=candidate.sets,
                     name.var.sets=NULL,candidate.cor.sets=candidate.cor.sets,joint=TRUE)
output2<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=candidate.sets,
                        name.var.sets=NULL,candidate.cor.sets=candidate.cor.sets,joint=TRUE)
test_that("output equa, given joints=true",{expect_equal(output1,output2)})


output1<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=NULL,
                        name.var.sets=name.var.sets,candidate.cor.sets=candidate.cor.sets,joint=TRUE)
output2<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=NULL,
                        name.var.sets=name.var.sets,candidate.cor.sets=candidate.cor.sets,joint=TRUE)
test_that("output equa, given joints=true",{expect_equal(output1,output2)})


output1<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=candidate.sets,
                        name.var.sets=NULL,candidate.cor.sets=candidate.cor.sets,joint=FALSE)
output2<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=candidate.sets,
                        name.var.sets=NULL,candidate.cor.sets=candidate.cor.sets,joint=FALSE)
test_that("output equa, given joints=FALSE",{expect_equal(output1,output2)})


output1<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=NULL,
                        name.var.sets=name.var.sets,candidate.cor.sets=candidate.cor.sets,joint=FALSE)
output2<-QICc.procedure(x=x,y=y,id=id,dist=dist,candidate.sets=NULL,
                        name.var.sets=name.var.sets,candidate.cor.sets=candidate.cor.sets,joint=FALSE)
test_that("output equa, given joints=FALSE",{expect_equal(output1,output2)})



