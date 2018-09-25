### Brownian kernel
cts = function(t,s){
  return(min(t,s))
}

tfs = function(i,j){
  out_put = function(t,s){
    cts(t,s) * phi[[i]](t) * phi[[j]](s)
  }
  return(out_put)
}

y1 = function(x){
  return(cbind(cos(2*pi*x), sin(2*pi*x)))
}

y2 = function(x,y){
  p1 = sqrt(1 - (1 - 2*x)^2) * y
  p2 = 1 - 2 * x
  return(cbind(p1, p2))
}

yd = function(x, y, d){
  p1 = sqrt(1 - (1 - 2 * Rbeta.inv(x, d/2, d/2))^2) * y
  p2 = 1 - 2 * Rbeta.inv(x, d/2, d/2)
  return(cbind(p1, p2))
}

gen.Sphere = function(N, n){
  rds = runif.sobol(N, n)
  prj1 = y1(rds[,1])
  prj2 = y2(rds[,2], prj1)
  tp = prj2
  for(i in 3:(n-1)){
    np = yd(rds[,i], tp, i)
    tp = np
  }
  return(tp)
}

to.Fourier.Func = function(tcof,np){
  t_dim = length(tcof)
  times = seq(0,1,len = np)
  temp_bs = create.fourier.basis(rangeval = c(0,1),nbasis = t_dim)
  temp_coefs = as.matrix(Diagonal(t_dim))
  temp_fd = fd(temp_coefs, temp_bs)
  temp_vs = t(eval.fd(times, temp_fd))
  phi_t = t(tcof %*% temp_vs)
  return(phi_t)
}

to.Function.Adaptive = function(tcof,temp_vs){
  phi_t = t(tcof %*% temp_vs)
  return(phi_t)
}

dat.Proj.Dt = function(X,S,temp_vs){
  uts = S %*% temp_vs
  projs = uts %*% t(X)
  
  m1 = rowMeans(projs)
  m2 = rowMeans(sweep(projs,1,m1)^2)
  m3 = rowMeans(sweep(projs,1,m1)^3)
  m4 = rowMeans(sweep(projs,1,m1)^4)
  
  b1 = (m3/m2^(3/2))^2
  b2 = (m4/m2^2)
  
  jbs = nrow(X) * b1/6 + nrow(X) * (b2 - 3)^2/24
  
  sks = cbind(jbs, seq(1:nrow(S)))
  
  return(sks)
}

med.From.Sphere = function(theta, i){
  res = 1
  if(i <= length(theta)){
    temp_theta = theta[1:i]
    if(i == 1){
      res = cos(temp_theta)
    } else {
      for(p in 1:(i-1)){
        res = res * sin(temp_theta[p])
      }
        res = res * cos(temp_theta[i])
    }
  } else {
    temp_theta = theta
    for(p in 1:(i-2)){
      res = res * sin(temp_theta[p])
    }
    res = res * sin(temp_theta[i-1])
  }
  return(res)
}

from.Sphere = function(theta){
  s = rep(0,(length(theta)+1))
  for(i in 1:length(s)){
    s[i] = med.From.Sphere(theta, i)
  }
  return(s)
}

to.Sphere = function(s){
  theta = rep(0, length(s)-1)
  theta[1] = acos(s[1])
  for(i in 2:length(theta)){
    theta[i] = acos(s[i] /sqrt(sum(s[i:length(s)]^2)))
  }
  if(s[length(s)] < 0){
    theta[length(theta)] = 2*pi - theta[length(theta)] 
  }
  return(theta)
}

sphere.Proj = function(theta,X,np){
  s = from.Sphere(theta)
  tf = to.Fourier.Func(s, np)
  projs = X %*% tf / np
  return(jarque.bera.test(projs)$statistic)
}

hybrid.Proj.Max.Dt = function(S, X, np, vs){
  reso = dat.Proj.Dt(X, S, np)
  w = reso[order(reso[,1],decreasing = T),]
  s = S[w[1:10,2],]
  cand_dirs = apply(s, MARGIN = 1, FUN = meta.Pso.Emb, tX = X, tnp = np, tvs = vs)
  return(max(cand_dirs[1,]))
}

hybrid.Proj.Dir.Dt = function(S, X, np, vs){
  reso = dat.Proj.Dt(X, S, vs)
  w = reso[order(reso[,1],decreasing = T),]
  s = S[w[1:10,2],]
  cand_dirs = apply(s, MARGIN = 1, FUN = meta.Pso.Emb, tX = X, tnp = np, tvs = vs)
  return(cand_dirs[2:nrow(cand_dirs),which.max(cand_dirs[1,])])
}

meta.Pso.Emb = function(ws,tX,tnp,tvs){
  theta = to.Sphere(ws)
  low_bd = theta - rep(0.1*pi, length(theta))
  up_bd = theta + rep(0.1*pi, length(theta))
  k = pso::psoptim(theta,sphere.Proj.Emb, X = tX, np = tnp, vs = tvs,
                   lower = low_bd, upper = up_bd,
                   control = list(fnscale = -1, vectorize = T, maxit = 1e3))
  return(c(abs(k$value), k$par))
}

sphere.Proj.Emb = function(theta,X,np,vs){
  s = from.Sphere(theta)
  tf = s %*% vs
  projs = X %*% t(tf) / np
  return(jarque.bera.test(projs)$statistic)
}

hybrid.Proj.Max.Emb = function(reso, S, X, np, vs){
  w = reso[order(reso[,1],decreasing = T),]
  s = S[w[1:10,2],]
  cand_dirs = apply(s, MARGIN = 1, FUN = meta.Pso.Emb, tX = X, tnp = np, tvs = vs)
  return(max(cand_dirs[1,]))
}

hybrid.Proj.Dir.Emb = function(reso, S, X, np){
  w = reso[order(reso[,1],decreasing = T),]
  s = S[w[1:10,2],]
  cand_dirs = apply(s, MARGIN = 1, FUN = meta.Pso.Emb, tX = X, tnp = np, tvs = vs)
  return(cand_dirs[2:nrow(cand_dirs),which.max(cand_dirs[1,])])
}

multi.std = function(mat,nfbs){
  x = scale(mat,scale = F)
  #q = cov.shrink(x, verbose = F)
  nq = nearPD(cov(x))
  q = matrix(nq$mat,ncol = nfbs)
  ql = t(chol(q))
  iql = solve(ql)
  zn = t(iql %*% t(x))
  #zr = matrix(zn@x, ncol = np)
  return(zn)
}


####### GHHK method for comparison ##########
GHHK = function(fd,k){
  proj_pca = pca.fd(fd, nharm = k)
  coefs = proj_pca$scores
  st_GHHK = rep(0,k)
  for(j in 1:k){
    prd = coefs[,j]
    tau = (sum((prd - mean(prd))^3) / length(prd)) / (sd(prd))^3
    kappa = (sum((prd - mean(prd))^4) / length(prd)) / (sd(prd))^4
    st_GHHK[j] = length(prd) * (tau^2 / 6 + (kappa - 3)^2 / 24)
  }
  return(sum(st_GHHK))
}

basis.Select = function(dat,stv,nsbs,nfbs){
  select_proj = dat %*% t(stv)
  select_jbs = rep(0, nsbs)
  for(i in 1:nsbs){
    select_jbs[i] = jarque.bera.test(select_proj[,i])$statistic
  }
  return(order(select_jbs,decreasing = T)[1:nfbs])
}

generate.test.bs = function(dat,stv,nsbs,nfbs, select_bs){
  select_proj = dat %*% t(stv)
  select_jbs = rep(0, nsbs)
  for(i in 1:nsbs){
    select_jbs[i] = jarque.bera.test(select_proj[,i])$statistic
  }
  selected_coefs = order(select_jbs,decreasing = T)[1:nfbs]
  m = max(selected_coefs)
  if(m%%2==0){m = m+1}
  TI = rep(1,m)
  TI[selected_coefs] = 0
  dps = which(TI == 1)
  temp_bs = create.fourier.basis(rangeval = c(0,1),dropind = dps,nbasis = m)
  return(temp_bs)
}