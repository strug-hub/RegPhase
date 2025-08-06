test_that("RoP returns proper p-values", {
  N=100
  Hap1<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
  Hap2<-matrix(rbinom(N*2,size = 1,prob=0.2),ncol=2)
  Y<-rnorm(n=N)
  out<-RoP(Y,Hap1,Hap2,family="gaussian")

  expect_named(out,c("p_cis","p_trans"))
  expect_true(all(out>=0 & out<=1))
})
