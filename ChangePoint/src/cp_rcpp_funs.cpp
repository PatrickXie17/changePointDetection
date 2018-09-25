#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

template <class T>
inline double do_mean( T& x ) {
  return mean(x);
}

NumericVector col_means( NumericMatrix& X ) {
  
  int nCols = X.ncol();
  NumericVector out = no_init(nCols);
  
  for( int j=0; j < nCols; j++ ) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = do_mean( tmp );
  }
  
  return out;
  
}

NumericVector row_means( NumericMatrix& X ) {
  
  int nRows = X.nrow();
  NumericVector out = no_init(nRows);
  
  for( int i=0; i < nRows; i++ ) {
    NumericMatrix::Row tmp = X(i, _);
    out[i] = do_mean( tmp );
  }
  
  return out;
  
}
// [[Rcpp::export]]
NumericMatrix row_erase (NumericMatrix& x, IntegerVector& rowID) {
  rowID = rowID.sort();
  
  NumericMatrix x2(Dimension(x.nrow()- rowID.size(), x.ncol()));
  int iter = 0; 
  int del = 1; // to count deleted elements
  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID[del - 1]) {
      x2.row(iter) = x.row(i);
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}


// NumericMatrix col_erase (NumericMatrix& x, IntegerVector& colID) {
//   colID = colID.sort();
//   
//   NumericMatrix x2(Dimension(x.nrow(), x.ncol()- colID.size()));
//   int iter = 0; 
//   int del = 1; 
//   for (int i = 0; i < x.ncol(); i++) {
//     if (i != colID[del - 1]) {
//       x2.col(iter) = x.column(i);
//       iter++;
//     } else {
//       del++;
//     }
//   }
//   return x2;
// }


// [[Rcpp::export]]
double Disco(NumericVector x, NumericVector y){
  int n = x.size();
  int m = y.size();
  double p1 = 0;
  double p2 = 0;
  double p3 = 0;
  for(int i=0; i<n; ++i){
  	for(int j=0; j<m; ++j){
  	  //p1 += pow(std::abs(x[i] - y[j]),1);
  	  p1 += std::abs(x[i] - y[j]);
  	}
  }

  for(int i=0; i<(n-1);++i){
    for(int k=(i+1); k < n;++k){
      //p2 += pow(std::abs(x[i] - x[k]),1);
      p2 += std::abs(x[i] - x[k]);
    }
  }
  
  for(int j =0; j<(m-1); ++j){
    for(int k=(j+1); k < m;++k){
      //p3 += pow(std::abs(y[j] - y[k]),1);
      p3 += std::abs(y[j] - y[k]);
    }
  }
  
  p1 = p1 * 2 / (n*m);
  p2 = p2 * 2 / (n*(n-1));
  p3 = p3 * 2 / (m*(m-1));
  return(std::max(0.0,p1-p2-p3));
}

// [[Rcpp::export]]
double Disco2(NumericVector x, NumericVector y){
  int n = x.size();
  int m = y.size();
  double p1 = 0;
  double p2 = 0;
  double p3 = 0;
  for(int i=0; i<n; ++i){
    for(int j=0; j<n; ++j){
      p1 += pow(std::abs(x[i] - y[j]),1);
      //p1 += std::abs(x[i] - y[j]);
    }
  }
  
  for(int i=0; i<(n-1);++i){
    for(int k=(i+1); k < n;++k){
      p2 += pow(std::abs(x[i] - x[k]),1);
      //p2 += std::abs(x[i] - x[k]);
    }
  }
  
  for(int j =0; j<(m-1); ++j){
    for(int k=(j+1); k < m;++k){
      p3 += pow(std::abs(y[j] - y[k]),1);
      //p3 += std::abs(y[j] - y[k]);
    }
  }
  
  p1 = p1 * 2 / (n*n);
  p2 = p2 * 2 / (n*(n-1));
  p3 = p3 * 2 / (m*(m-1));
  return(std::max(0.0,p1-p2-p3));
}


// [[Rcpp::export]]
NumericMatrix rollingWindowGenerator(NumericVector dat, int window_size){
  int nrows = dat.size()-window_size+1;
  NumericMatrix tempdat(nrows, window_size);
  for(int i=0; i<nrows; ++i){
    tempdat(i,_) = dat[Range(i, (i+window_size-1))];
  }
  return(tempdat);
}

// [[Rcpp::export]]
NumericMatrix rollingWindowDistance(NumericMatrix x){
  int dr = x.nrow();
  NumericMatrix res(dr,dr);
  for(int i=0; i< dr; ++i){
    for(int j =0; j< dr; ++j){
      res(i,j) = Disco(x(i,_), x(j,_));
    }
  }
  
  return(res);
}

// [[Rcpp::export]]
NumericMatrix rollingWindowDistance2(NumericMatrix x){
  int dr = x.nrow();
  NumericMatrix res(dr,dr);
  for(int i=0; i< dr; ++i){
    for(int j =0; j< dr; ++j){
      if(j >= i){
        res(i,j) = Disco(x(i,_), x(j,_));
      } else {
        res(i,j) = res(j,i);
      }
    }
  }
  return(res);
}

// [[Rcpp::export]]
NumericMatrix PKM_Distance(NumericMatrix ctrds, NumericMatrix dat){
  int n = ctrds.ncol();
  int m = dat.nrow();
  NumericMatrix res(m,n);
  for(int i=0; i<m; ++i){
    for(int j=0; j<n; ++j){
      res(i,j) = Disco(ctrds(_,j), dat(i,_));
    }
  }
  return(res);
}

// // [[Rcpp::export]]
// double KMedoids_Cost(NumericVector A, NumericVector B, NumericMatrix dat){
//   int n = dat.nrow();
//   double res = 0;
//   for(int i=0; i<n; ++i){
//     res = res + Disco(A,dat(i,_)) + Disco(B,dat(i,_));
//   }
//   return(res);
// }

// // [[Rcpp::export]]
// NumericVector KMedoids_Update(int i1, int i2, double cost, NumericMatrix dat){
//   int n = dat.nrow();
//   int new_i1 = 0;
//   for(int i=0; i<n; i++){
//     if(i == i1 | i == i2){
//       break;
//     } else {
//       int temp_c = i;
//       double c1 = KMedoids_Cost()
//     }
//   }
// }


// [[Rcpp::export]]
NumericMatrix KMeans_Update(NumericMatrix dists, NumericMatrix dat){
  int n = dists.nrow();
  int m = dat.ncol();
  NumericMatrix res(m,2);
  NumericVector ind1 = NumericVector::create(-1);
  NumericVector ind2 = NumericVector::create(-1);
  for(int i=0; i<n;++i){
    if(dists(i,0) < dists(i,1)){
      ind1.push_back(i);
    } else {
      ind2.push_back(i);
    }
  }
  ind1.erase(0);
  ind2.erase(0);
  int na = ind1.size();
  int nb = ind2.size();
  NumericMatrix g1(na,m);
  for(int i=0; i<na; ++i){
    g1(i,_) = dat(ind1[i],_);
  }
  NumericMatrix g2(nb,m);
  for(int j=0; j<nb; ++j){
    g2(j,_) = dat(ind2[j],_);
  }
  NumericMatrix d1(na,na);
  d1 = rollingWindowDistance(g1);
  NumericMatrix d2(nb,nb);
  d2 = rollingWindowDistance(g2);
  res(_,0) = col_means(g1);
  res(_,1) = col_means(g2);
  return(res);
}

// // [[Rcpp::export]]
// NumericVector KMeans(NumericVector A, NumericVector B, NumericMatrix dat, int reps){
//   int n = A.size();
//   int m = dat.nrow();
//   NumericMatrix dists(m,2);
//   dists = KMeans_Distance(A,B,dat);
//   NumericVector v1(n);
//   NumericVector v2(n);
//   NumericMatrix vsk(n,2);
//   for(int i = 0; i<reps; ++i){
//      vsk = KMeans_Update(dists, dat);
//      v1 = vsk(_,0);
//      v2 = vsk(_,1);
//      dists = KMeans_Distance(v1, v2,dat);
//   }
//   NumericVector ind = NumericVector::create(-1);
//   for(int i=0; i<m;++i){
//     if(dists(i,0) < dists(i,1)){
//       ind.push_back(i+1);
//     }
//   }
//   ind.erase(0);
//   return(ind);
// }


// [[Rcpp::export]]
NumericVector col_Square_Sums(NumericMatrix dat){
  int n = dat.nrow();
  NumericVector res(n);
  for(int i=0; i<n; ++i){
    NumericVector tempdat = dat(i,_);
    for(int j=0; j<n; ++j){
      res[i] += pow(tempdat[j],2);
    }
  }
  return(res);
}

// [[Rcpp::export]]
NumericVector is_Min(NumericVector k, NumericMatrix td){
  int n = td.nrow();
  NumericVector res = NumericVector::create(-1);
  for(int i=0; i<n; ++i){
    if(k[i] == min(td(i,_))){
      res.push_back(i+1);
    }
  }
  res.erase(0);
  return(res);
}

// [[Rcpp::export]]
double cal_diff(NumericMatrix v1, NumericMatrix v2){
  int n = v1.ncol();
  double res = 0.0;
  for(int i=0; i<n; ++i){
    res += Disco(v1(_,i), v2(_,i));
  }
  return(res);
}

// [[Rcpp::export]]
NumericMatrix SOM_Diff(NumericMatrix dat, NumericMatrix vs){
  int n = dat.nrow();
  NumericMatrix res(n,2);
  for(int i=0; i< n; ++i){
    res(i,0) = Disco(dat(i,_), vs(0,_));
    res(i,1) = Disco(dat(i,_), vs(1,_));
  }
  return(res);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

# /*** R
# x = rnorm(10)
# x
# test(x)
# */