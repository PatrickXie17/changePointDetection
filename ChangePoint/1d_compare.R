library(Rcpp)
library(metRology)
#library(shallot)
library(fOptions)
library(ecp)

tms = Sys.time()

source("R/cp_funs2.R")
sourceCpp("src/cp2.cpp")

require("foreach")
require("parallel")
require("doParallel")
require("doMC")
registerDoMC(detectCores())

N1 = 1000
N2 = 1000
d = 300
true_label = c(rep(1,N1), rep(2,N2), rep(3, N1))

res = foreach(z=1:2000, .combine = 'rbind', .errorhandling = "remove") %dopar%{
  set.seed(z)
  X = rnorm(N1,0,1)
  Y1 = rnorm(N2,1,1)
  Y2 = rnorm(N2,0,2)
  Y3 = rt.scaled(N2,df = 2)
  Y4 = rt.scaled(N2,df = 16)
  Z1 = c(X,Y1,X)
  Z2 = c(X,Y2,X)
  Z3 = c(X,Y3,X)
  Z4 = c(X,Y4,X)
  ## Z1
  q1 = rollingWindowGenerator(Z1, d)

  s1 = SOM_R2(q1, 2000, 3)
  
  ze1 = fine_CP_finder(Z1, s1, d)
  
  MJ1 = e.divisive(as.matrix(Z1),k = 2)$estimates[2:3]

  estimated_label_1a = c(rep(1,ze1[1]), rep(2, (ze1[2] - ze1[1])),  rep(3, (length(Z1)-ze1[2])))
  estimated_label_1b = c(rep(1,MJ1[1]), rep(2, (MJ1[2] - MJ1[1])),  rep(3, (length(Z1)-MJ1[2])))
  
  w1a = adj.rand.index(true_label, estimated_label_1a)
  w1b = adj.rand.index(true_label, estimated_label_1b)
  
  ## Z2
  q2 = rollingWindowGenerator(Z2, d)
  
  s2 = SOM_R2(q2, 2000, 3)
  
  ze2 = fine_CP_finder(Z2, s2, d)
  
  MJ2 = e.divisive(as.matrix(Z2),k = 2)$estimates[2:3]
  
  estimated_label_2a = c(rep(1,ze2[1]), rep(2, (ze2[2] - ze2[1])),  rep(3, (length(Z2)-ze2[2])))
  estimated_label_2b = c(rep(1,MJ2[1]), rep(2, (MJ2[2] - MJ2[1])),  rep(3, (length(Z2)-MJ2[2])))
  
  w2a = adj.rand.index(true_label, estimated_label_2a)
  w2b = adj.rand.index(true_label, estimated_label_2b)
  
  ## Z3
  q3 = rollingWindowGenerator(Z3, d)
  
  s3 = SOM_R2(q3, 2000, 3)
  
  ze3 = fine_CP_finder(Z3, s3, d)
  
  MJ3 = e.divisive(as.matrix(Z3),k = 2)$estimates[2:3]
  
  estimated_label_3a = c(rep(1,ze3[1]), rep(2, (ze3[2] - ze3[1])),  rep(3, (length(Z3)-ze3[2])))
  estimated_label_3b = c(rep(1,MJ3[1]), rep(2, (MJ3[2] - MJ3[1])),  rep(3, (length(Z3)-MJ3[2])))
  
  w3a = adj.rand.index(true_label, estimated_label_3a)
  w3b = adj.rand.index(true_label, estimated_label_3b)
  
  
  ## Z4
  q4 = rollingWindowGenerator(Z4, d)
  
  s4 = SOM_R2(q4, 2000, 3)
  
  ze4 = fine_CP_finder(Z4, s4, d)
  
  MJ4 = e.divisive(as.matrix(Z4),k = 2)$estimates[2:3]
  
  estimated_label_4a = c(rep(1,ze4[1]), rep(2, (ze4[2] - ze4[1])),  rep(3, (length(Z4)-ze4[2])))
  estimated_label_4b = c(rep(1,MJ4[1]), rep(2, (MJ4[2] - MJ4[1])),  rep(3, (length(Z4)-MJ4[2])))
  
  w4a = adj.rand.index(true_label, estimated_label_4a)
  w4b = adj.rand.index(true_label, estimated_label_4b)
  
  c(w1a, w1b,w2a, w2b, w3a, w3b, w4a, w4b)
}
#print(res)

#write.csv(res,'compare_1d_a.csv')
saveRDS(res, file="compare_1d_2cp.Rds")

colMeans(res)
for(i in 1:6){
  print(mean(res[,i]))
  print(sd(res[,i]))
}

Sys.time() - tms