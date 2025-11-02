#include <Rcpp.h>

using namespace Rcpp;

// RCPP implementation of random walk sampling
//
// [[Rcpp::export]]
IntegerVector crwsample(IntegerVector b, IntegerVector c, int m) {
  int n = c.size() - 1;
  IntegerVector out(m);
  int tic = 0; 
  while(tic < m) {
    int pi = sample(n, 1)[0];
    IntegerVector s = seq(c[pi - 1], c[pi] - 1);
    int temp = b[sample(s, 1)[0]] + 1;
    if(std::find(out.begin(), out.end(), temp) == out.end()) {
      out[tic] = temp;
      ++tic;
    } 
  }
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

