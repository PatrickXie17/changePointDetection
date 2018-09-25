library(Rcpp)
library(metRology)
library(shallot)
library(fOptions)
library(ecp)

source("R/cp_funs2.R")
sourceCpp("src/cp2.cpp")

require("foreach")
require("parallel")
require("doParallel")
require("doMC")
registerDoMC(detectCores())

set.seed(3)
N1 = 500
N2 = 500
N3 = 500
d = 200

#D = rt.scaled(N1,df = 2)
X = rnorm(N1,0,1)
#Y = rnorm(N2,1,1)
Y = rnorm(N2,1,1)
#Y = rt.scaled(N2,df = 16)
W = rnorm(N3, 0,1)
#Y = rt(N2, df = 16)
#Y = rgamma(N2,2,1)
Z = c(X,Y,W)
#c1 = e.cp3o(as.matrix(Z, ncol = 1))

q = rollingWindowGenerator(Z, d)
#dm = rollingWindowDistance2(q)

set.seed(123)

res = SOM_R2(q, 2000, 3)
res
td = res[1,2]

tdat = Z[(td-d):(td+d)]
ze = e.divisive(as.matrix(tdat),k = 1)$estimates[2] + (td-d)
e.divisive(as.matrix(Z),k = 2)$estimates

f1 = function(){
  s1 = SOM_R(Z, q, 100)
  td1 = s1[1,2]
  tdat1 = Z[(td1-d):(td1+d)]
  ze1 = e.divisive(as.matrix(tdat1),k = 1)$estimates[2] + (td1-d)
}

f2 = function(){
  e.divisive(as.matrix(Z),k = 1)$estimates[2]
}

microbenchmark::microbenchmark(f1(),f2(), times = 10)



N1 = 1000
N2 = 1000
X = rnorm(N1,0,1)
Y = rt(N2,df = 2)
Y = (Y - mean(Y))/sd(Y)
Z = c(X,Y)
d = 200
q = rollingWindowGenerator(Z, d)
s1 = SOM_R(Z, q, 100)
td1 = s1[1,2]
tdat1 = Z[(td1-d):(td1+d)]
ze1 = e.divisive(as.matrix(tdat1),k = 1)$estimates[2] + (td1-d)








w = dm[,cs[order(cs)]]
#s = KMedoids_R(kdr,q,dm)


true_label = c(rep(1,N1), rep(2,N2))
estimated_label = rep(0,(N1+N2))

if(min(s) == 1){
  ind = min(s):(max(s)+d)
  estimated_label[ind] = 1
  estimated_label[-ind] = 2
} else {
  ind = (min(s)):(max(s)+(d-1))
  estimated_label[-ind] = 1
  estimated_label[ind] = 2
}

adj.rand.index(true_label, estimated_label)

plot(Z, main = 'N(0,1) vs t(16)')
lines(ind, Z[ind], col = 'red', type = 'p')
abline(v = 600.5)

ks.cp3o(as.matrix(Z, ncol = 1),K = 1)

newvecs = KMeans_Update(w,q)
nv2 = KMeans_Update_R(q,w)

tms = Sys.time()
d1 = rep(0,1000)
d2 = rep(0,1000)
for(i in 1:1000){
  set.seed(i)
  x = rnorm(100)
  y1 = rnorm(100)
  #y2 = rt.scaled(50,df = 2)
  y2 = rt(100, df = 2)
  #y2 = rnorm(100,0,sqrt(2))
  #y2 = rnorm(50, mean = 0, sd = sqrt(10))
  
  d1[i] = Disco(x,y1)
  d2[i] = Disco(x,y2)
}
Sys.time() - tms

sum(d1 < d2)
