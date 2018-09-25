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
NumericMatrix rollingWindowGenerator(NumericVector dat, int window_size){
  int nrows = dat.size()-window_size+1;
  NumericMatrix tempdat(nrows, window_size);
  for(int i=0; i<nrows; ++i){
    tempdat(i,_) = dat[Range(i, (i+window_size-1))];
  }
  return(tempdat);
}

// [[Rcpp::export]]
NumericMatrix SOM_Diff(NumericMatrix dat, NumericMatrix vs, int K){
  int n = dat.nrow();
  NumericMatrix res(n,K);
  for(int i=0; i< n; ++i){
    for(int j=0; j<K; ++j){
      res(i,j) = Disco(dat(i,_), vs(j,_));
    }
  }
  return(res);
}

// [[Rcpp::export]]
NumericVector SOM_Selector_min(NumericVector tc, NumericMatrix tdis){
  int n = tdis.nrow();
  int k = tdis.ncol();
  NumericVector ind = NumericVector::create(-1);
  NumericVector ms(k);
  for(int l=0; l<k; ++l){
    ms(l) = mean(tdis(l,_));
  }
  
  double mms = min(ms);
  for(int i = 0; i<n; ++i){
    if(tc[i] < mms){
      ind.push_back(i+1);
    }
  }
  ind.erase(0);
  return(ind);
}

