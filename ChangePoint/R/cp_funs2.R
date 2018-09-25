adj.rand.index = function (c1, c2) 
{
  n <- length(c1)
  if (length(c2) != n) 
    stop("Clusterings must be the same length.")
  t1 <- table(c1)
  t2 <- table(c2)
  t12 <- table(c1, c2)
  expected <- sum(choose(t1, 2)) * sum(choose(t2, 2))/choose(n, 
                                                             2)
  numerator <- sum(choose(t12, 2)) - expected
  denominator <- 0.5 * (sum(choose(t1, 2)) + sum(choose(t2, 
                                                        2))) - expected
  numerator/denominator
}


SOM_Update_R = function(Z, dat, vs){
  d = ncol(dat) # window size
  # randomly select one window from all possible rolling windows
  si = sample(1:nrow(dat),1)
  tc = dat[si,]
  # calculate the distance between the selected window and two nodes
  d1 = Disco(tc,vs[1,])
  d2 = Disco(tc,vs[2,])
  # select the winning node
  tw = vs[which.min(c(d1,d2)),]
  
  # select the losing node
  tl = vs[which.max(c(d1,d2)),]
  
  # weighted randomly sampling from the selected data
  
  # for winning node, 80% from the winning node
  k1 = sample(1:d, 0.2*d)
  k2 = sample(1:d, 0.8*d)
  nd1 = c(tc[k1], tw[k2])
  
  # for losing node, 100% from the losing node
  nd2 = tl
  
  res = rbind(nd1, nd2)
  return(res)
}

SOM_R = function(Z, dat, iters = 10){
  new_vs = dat[sample(1:nrow(dat),2),]
  for(i in 1:iters){
    new_vs = SOM_Update_R(Z, dat, new_vs)
  }
  dis = SOM_Diff(dat, new_vs)
  if(dis[1,1] > dis[1,2]){
    dis = cbind(dis[,2], dis[,1])
  }
  p = which(dis[,1] < dis[,2])
  if(p[1] != 1){
    p = c(1:nrow(dat))[-p]
  }
  qd = diff(p)
  td = which(qd>1)
  if(length(td) != 0){
    cp = td[which(td>100)][1]
  } else {
    cp = max(p)
  }
  res = matrix(ncol = 2, nrow = 2)
  res[1,1] = 1
  res[1,2] = cp
  res[2,1] = cp+1
  res[2,2] = nrow(dat)
  return(res)
}


######### for multiple change points ##############

SOM_Update_R2 = function(dat, vs, K){
  d = ncol(dat) # window size
  # randomly select one window from all possible rolling windows
  si = sample(1:nrow(dat),1)
  tc = dat[si,]
  # calculate the distance between the selected window and nodes
  ds = rep(0, K)
  for(i in 1:K){
    ds[i] = Disco(tc,vs[i,])
  }
  # select the winning node
  wi = which.min(ds)
  tw = vs[wi,]
  
  # weighted randomly sampling from the selected data
  
  # for winning node, 80% from the winning node
  k1 = sample(1:d, 0.2*d)
  k2 = sample(1:d, 0.8*d)
  nd1 = c(tc[k1], tw[k2])
  
  # for losing node, 100% from the losing node
  vs[wi,] = nd1
  
  res = vs
  return(res)
}

# SOM_Locator = function(dat, nw, tds){
#   tdat = dat[nw,]
#   k = ncol(tds)
#   lc = rep(0,k)
#   for(p in 1:k){
#     lc[k] = Disco(tdat, tds[p,])
#   }
# }

SOM_Selector = function(dis, K){
  tds = dis
  mms = max(colMeans(dis))
  tres = matrix(1, ncol = 2, nrow = K-1)
  ind = 1
  while(ind < K){
    ti = which.min(colMeans(tds[1:50,]))
    tc = tds[,ti]
    nind = SOM_Selector_min_R(tc, tds,mms)
    qd = diff(nind)
    td = which(qd>1)
    if(length(td) != 0){
      cp = td[which(td>100)][1]
    } else {
      cp = max(nind)
    }
    tres[ind,] = c(1, cp)
    tds = tds[-c(tres[ind,1]:tres[ind,2]),-ti]
    ind = ind+1
  }
  res = matrix(ncol = 2, nrow = K)
  res[1,] = tres[1,]
  for(i in 2:(K-1)){
    res[i,] = tres[i,] + res[(i-1),2]
  }
  res[K,] = c((res[(K-1),2]+1), nrow(dis))
  return(res)
}

SOM_R2 = function(dat, iters = 10, K){
  new_vs = dat[sample(1:nrow(dat),K),]
  for(i in 1:iters){
    new_vs = SOM_Update_R2(dat, new_vs, K)
  }
  dis = SOM_Diff(dat, new_vs,K)
  
  ans = SOM_Selector(dis, K)

  return(ans)
}

SOM_Selector_min_R = function(tc, tdis,mms){
  ind = c()
  for(i in 1:nrow(tdis)){
    if(tc[i] < mms){
      ind = c(ind,i)
    }
  }
  return(ind)
}

fine_CP_finder = function(Z, labs, d){
  k = nrow(labs) - 1
  res = rep(0,k)
  for(i in 1:k){
    td = labs[i,2]
    tdat = Z[(td-d):(td+d)]
    res[i] = e.divisive(as.matrix(tdat),k = 1)$estimates[2] + (td-d)
  }
  return(res)
}