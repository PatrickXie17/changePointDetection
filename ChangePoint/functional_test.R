require('zipfR')
require("MASS")
require('fOptions')
require('fda')
require("foreach")
require("tseries")
require("parallel")
require("doParallel")
require("doMC")
library(ecp)
library(fChange)

# require("Matrix")
# require("pso")
# require("corpcor")

tms = Sys.time()
registerDoMC(detectCores())
#registerDoMC(50)
source("R/high_dim_funs.R")
source("R/cp_funs.R")
sourceCpp("src/cp_rcpp_funs.cpp")

tms = Sys.time()

S = gen.Sphere(1e4,7)

N = 150
nfbs = 7
n_ghhk = 3
K = 11
n_select = 7
np = 75 # 75 observations
times = seq(0,1,len=np) # t \in [0,1]

dgp_bs = create.bspline.basis(rangeval = c(0,1),nbasis = K)
#dgp_bs = create.fourier.basis(rangeval = c(0,1),nbasis = K)
test_bs = create.fourier.basis(rangeval = c(0,1),nbasis = nfbs)

I1 = as.matrix(Diagonal(K,seq(1,0.1,length.out=K)))
I2 = as.matrix(Diagonal(K,seq(1,0.1,length.out=K)))
I2[5,5] = I2[5,5]/2
set.seed(2)
dlt1 = mvrnorm(N, mu = rep(0, K), I1)
dlt2 = mvrnorm(N, mu = c(rep(0, 10),rep(0,1)), I2)
# tz = rt(N,2)
# tz = tz / sd(tz)
# dlt2[,5] = tz
dlts = t(rbind(dlt1, dlt2))
GP_fd = fd(dlts, dgp_bs)
k1 = change_fPCA(GP_fd,d = 3)$change
k2 = change_FF(GP_fd)$change
values = t(eval.fd(times, GP_fd))

res = foreach(z=1:nrow(S), .combine = 'rbind', .errorhandling = "remove") %dopar%{
  ts = S[z,]
  tfd = fd(ts, test_bs)
  tv = eval.fd(times, tfd)
  tprj = values %*% tv
  cp_test = e.divisive(tprj,k = 1)
  pvs = cp_test$estimates[2]
  Q = Disco(tprj[1:pvs], tprj[-c(1:pvs)])
  c(pvs,Q)
}
res[order(res[,2], decreasing = T)[1:20],]
write.csv(res, 'f_res_1.csv')
Sys.time()-tms