#include <Rcpp.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
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

