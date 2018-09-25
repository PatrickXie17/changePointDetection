library(Rcpp)
library(metRology)
#library(shallot)
library(fOptions)
library(ecp)

tms = Sys.time()

source("R/cp_funs.R")
sourceCpp("src/cp_rcpp_funs.cpp")

require("foreach")
require("parallel")
require("doParallel")
require("doMC")
registerDoMC(detectCores())

set.seed(3)
N1 = 500
N2 = 1000
d = 100

set.seed(2)
X = rnorm(N1,0,1)
Y1 = rnorm(N2,1,1)
Y2 = rnorm(N2,0,sqrt(2))
#Y = rt.scaled(N2,df = 2)
#Y = rt(N2, df = 16)
#Y = rgamma(N2,2,1)
Z1 = c(X,Y1)
Z2 = c(X,Y2)

true_label = c(rep(1,N1), rep(2,N2))
MJ2 = e.divisive(as.matrix(Z2, ncol = 1),k = 1)$estimates
estimated_label_2b = c(rep(1,MJ2), rep(2,length(true_label)-MJ2))

print(MJ2)

q2 = rollingWindowGenerator(Z2, d)
dm2 = rollingWindowDistance(q2)

qs = sample(1:nrow(q2),2)
c1 = q2[qs[1],]
c2 = q2[qs[2],]

s2 = KMeans_R(c1,c2,q2,dm2)


estimated_label_2a = rep(0,(N1+N2))

if(min(s2) == 1){
  ind2 = min(s2):(max(s2)+d)
  estimated_label_2a[ind2] = 1
  estimated_label_2a[-ind2] = 2
  MJ1 = max(s2)
} else {
  ind2 = (min(s2)):(max(s2)+(d-1))
  estimated_label_2a[-ind2] = 1
  estimated_label_2a[ind2] = 2
  MJ1 = min(s2)
}

w2a = adj.rand.index(true_label, estimated_label_2a)

w2b = adj.rand.index(true_label, estimated_label_2b)

c(MJ1, MJ2)
