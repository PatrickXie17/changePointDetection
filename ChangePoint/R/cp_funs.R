rw = function(dat, ws){
  tempdat = matrix(0, ncol = ws, nrow = (length(dat) - ws +1))
  for(i in 1:nrow(tempdat)){
    tempdat[i,] = dat[i:(i+ws-1)]
  }
  res = rep(0, nrow(tempdat)-1)
  for(j in 1:length(res)){
    res[j] = Discor(tempdat[j,], tempdat[(j+1),])
  }
  return(res)
}

Discor = function(x,y){
  n = length(x)
  m = length(y)
  p1 = 0
  p2 = 0
  p3 = 0
  for(i in 1:n){
    for(j in 1:m){
      p1 = p1 + abs(x[i] - y[j])
    }
  }
  
  for(i in 1:(n-1)){
    for(k in (i+1):n){
      p2 = p2 + abs(x[i] - x[k])
    }
  }
  
  for(j in 1:(m-1)){
    for(l in (j+1):m){
      p3 = p3 + abs(y[j] - y[l])
    }
  }
  
  p1 = p1 * 2 / (m*n)
  p2 = p2 * 2 / (n*(n-1))
  p3 = p3 * 2 / (m*(m-1))
  return(p1 - p2 - p3)
}

KMeans_Distance_R = function(d1, d2, dat, dm){
  i1 = which(dat[,1] == d1[1])
  i2 = which(dat[,1] == d2[1])
  return(dm[,c(i1,i2)])
}


# K-medoids
KMD_Distance_R = function(kdr, dat, dm){
  return(dm[,match(kdr[,1], dat[,1])])
}

KMeans_Update_R = function(dat, dists, dm){
  tcq = which(dists[,1] < dists[,2])
  d1 = min(tcq)
  d2 = nrow(dat) - max(tcq)
  if(d1 < d2){
    dm1 = dm[1:max(tcq),1:max(tcq)]
    dm2 = dm[(max(tcq)+1):nrow(dat),(max(tcq)+1):nrow(dat)]
    g1 = dat[1:max(tcq),]
    g2 = dat[(max(tcq)+1):nrow(dat),]
  } else {
    dm1 = dm[(d1+1):nrow(dat),(d1+1):nrow(dat)]
    dm2 = dm[1:d1,1:d1]
    g1 = dat[(d1+1):nrow(dat),]
    g2 = dat[1:d1,]
  }
  # if(d1 < d2){
  #   g1 = dat[1:max(tcq),]
  #   g2 = dat[(max(tcq)+1):nrow(dat),]
  # } else {
  #   g1 = dat[(d1+1):nrow(dat),]
  #   g2 = dat[1:d1,]
  # }
  #ind2 = which(dists[,1] >= dists[,2])
  # ind2 = seq(1:nrow(dat))[-ind1]
  # g1 = dat[ind1,]
  # g2 = dat[ind2,]
  # dm1 = rollingWindowDistance(g1)
  # dm2 = rollingWindowDistance(g2)
  
  k1 = which.min(colSums(dm1))
  k2 = which.min(colSums(dm2))
  
  res = matrix(0, ncol = 2, nrow = ncol(dat))
  res[,1] = g1[k1,]
  res[,2] = g2[k2,]
  return(res)
}

KMD_Get_cq = function(dat, dists, dm){
  n = ncol(dists)
  flag = TRUE
  fi = 1
  while(flag){
    ind = match(dists[fi,], dm[,fi])
    flag = (length(unique(ind)) != n)
    fi= fi + 1
  }
  ind = ind[order(ind)]
  cq = matrix(nrow = n, ncol = 2)
  td = dists
  for (i in 1:n){
    tcq = range(is_Min(td[,i], td))
    if(i == 1){
      cq[1,] = c(1, min((ind[i+1]-1), tcq[2]))  
    } else {
      if(i == n){
        cq[i,] = c((cq[(i-1),2]+1), nrow(dat))
      } else {
        cq[i,] = c((cq[(i-1),2]+1), min((ind[i+1]-1), tcq[2]))
      }
    }
  }
  return(cq)
}

KMD_Update_R = function(dat, dists, dm){
  n = ncol(dists)
  flag = TRUE
  fi = 1
  while(flag){
    ind = match(dists[fi,], dm[,fi])
    flag = (length(unique(ind)) != n)
    fi= fi + 1
  }
  ind = ind[order(ind)]
  cq = matrix(nrow = n, ncol = 2)
  td = dists
  for (i in 1:n){
    tcq = range(is_Min(td[,i], td))
    if(i == 1){
      cq[1,] = c(1, min((ind[i+1]-1), tcq[2]))  
    } else {
      if(i == n){
        cq[i,] = c((cq[(i-1),2]+1), nrow(dat))
      } else {
        cq[i,] = c((cq[(i-1),2]+1), min((ind[i+1]-1), tcq[2]))
      }
    }
  }
  res = matrix(0, ncol = n, nrow = ncol(dat))
  for(j in 1:n){
    tm = dm[cq[j,1]:cq[j,2], cq[j,1]:cq[j,2]]
    td = dat[cq[j,1]:cq[j,2], ]
    res[,j] = td[which.min(colSums(tm)),]
  }
  return(res)
}

KMeans_R = function(d1, d2, dat, dm){
  dists = KMeans_Distance_R(d1,d2,dat,dm)
  ov1 = dists[,1]
  ov2 = dists[,2]
  diff1 = 1
  diff2 = 1
  while((diff1 != 0) & (diff2 != 0 )){
    vsk = KMeans_Update_R(dat, dists,dm)
    v1 = vsk[,1]
    v2 = vsk[,2]
    dists = KMeans_Distance_R(v1, v2,dat,dm)
    diff1 = Disco(ov1, v1)
    diff2 = Disco(ov2, v2)
    ov1 = v1
    ov2 = v2
  }
  # for(i in 1:reps){
  #   vsk = KMeans_Update_R(dat, dists)
  #   v1 = vsk[,1]
  #   v2 = vsk[,2]
  #   dists = KMeans_Distance(v1, v2,dat)
  # }
  m2 = which(dat[,1] == v2[1])
  cq = which(dists[,1] < dists[,2])
  d1 = min(cq)
  d2 = nrow(dat) - max(cq)
  if(d1 < d2){
    res = c(1:max(cq))
  } else {
    res = c(d1:nrow(dat))
  }
  if(m2 %in% seq(min(cq),max(cq))){
    res = seq(1,nrow(dat))[-res]
  }
  return(res)
}

KMedoids_R = function(kdr, dat, dm){
  dists = KMD_Distance_R(kdr,dat,dm)
  ods = dists
  diff = 1
  dd = 0
  while(diff > 0){
    vsk = t(KMD_Update_R(dat, dists,dm))
    dists = KMD_Distance_R(vsk,dat,dm)
    diff = cal_diff(ods, dists)
    ods = dists
  }
  cq = KMD_Get_cq(dat, dists, dm)
  return(cq)
}

ODC_Update_R = function(cs, dat, dm){
  n = length(cs)
  qs = cs[order(cs)]
  res = ODC_RG_R(qs, dat, dm)
  mds = rep(0, nrow(res))
  for(i in 1:nrow(res)){
    mds[i] = res[i,1] + which.min(colSums(dm[res[i,1]:res[i,2],res[i,1]:res[i,2]])) - 1#floor(median(res[i,1]:res[i,2]))
  }
  return(mds)
}

ODC_RG_R = function(qs, dat, dm){
  dists = dm[,qs]
  res = matrix(ncol= 2, nrow = n)
  for (i in 1:n){
    if(i == 1){
      res[i,1] = 1
      ts = qs[i]
      for(j in qs[i]:qs[i+1]){
        if(dists[j,i] <= dists[j,(i+1)]){
          ts = j
        }
      }
      res[i,2] = ts
    }
    if(i > 1 & i < n){
      res[i,1] = res[(i-1), 2] + 1
      ts = qs[i]
      for(j in qs[i]:qs[i+1]){
        if(dists[j,i] <= dists[j,(i+1)]){
          ts = j
        }
      }
      res[i,2] = ts
    }
    if(i == n){
      res[i,1] = res[(i-1), 2] + 1
      res[i,2] = nrow(dat)
    }
  }
  return(res)
}

ODC_R = function(cs, dat, dm){
  diff = 10
  ndiff = 5
  ocs = cs[order(cs)]
  while(ndiff < diff){
    diff = ndiff
    ncs = ODC_Update_R(ocs, dat, dm)
    ndiff = max(abs(ncs - ocs))
    ocs = ncs
  }
  res = ODC_RG_R(ocs, dat, dm)
  return(res)
}

zero_in_window = function(x){
  flag = sum(x==0)
  return(flag)
}

PKM_RG_R = function(ctrds, dat){
  dists = PKM_Distance(ctrds, dat)
  n = ncol(ctrds)
  res = matrix(ncol= 2, nrow = n)
  for (i in 1:n){
    if(i == 1){
      res[i,1] = 1
      ts = res[i,1] + 1
      lb = mean(dists[,i+1])
      while((min(dists[ts:(ts+ncol(dat)),i])) < lb & ((ts+ncol(dat)) < nrow(dat))){
        ts = ts+1
      }
      res[i,2] = ts
    }
    if(i > 1 & i < n){
      res[i,1] = res[(i-1), 2] + 1
      ts = res[i,1] + 1
      lb = mean(dists[,i+1])
      while((min(dists[ts:(ts+ncol(dat)),i])) < lb & ((ts+ncol(dat)) < nrow(dat))){
        ts = ts+1
      }
      res[i,2] = ts
    }
    if(i == n){
      res[i,1] = res[(i-1), 2] + 1
      res[i,2] = nrow(dat)
    }
  }
  return(res)
}

PKM_Update_R = function(ctrds, dat, Z){
  res = PKM_RG_R(ctrds, dat)
  mds = matrix(ncol = ncol(ctrds), nrow = nrow(ctrds))
  for(i in 1:ncol(mds)){
    tz = Z[res[i,1]:(res[i,2]+ncol(dat)-1)]
    mds[,i] = sample(tz, nrow(ctrds))
  }
  return(mds)
}

PKM_R = function(ctrds, dat, Z){
  octrds = ctrds
  count = 0
  while(count < 10){
    nctrds = PKM_Update_R(octrds, dat, Z)
    count = count + 1
  }
  res = PKM_RG_R(nctrds, dat)
  return(res)
}

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
  #sw = which(dat[,1] == tw[1])
  
  # select the losing node
  tl = vs[which.max(c(d1,d2)),]
  #sl = which(dat[,1] == tl[1])
  
  # # back to the original data
  # zi = Z[si:(si+d-1)]
  # zw = Z[sw:(sw+d-1)]
  # zl = Z[sl:(sl+d-1)]
  
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