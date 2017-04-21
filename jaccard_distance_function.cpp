#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
using namespace std;

using namespace arma;
using namespace Rcpp;


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

NumericVector intersect(uvec x, uvec y) {
  std::vector<double> res;
  std::unordered_set<double> s(y.begin(), y.end());
  for (int i=0; i < x.size(); ++i) {
    auto f = s.find(x[i]);
    if (f != s.end()) {
      res.push_back(x[i]);
      s.erase(f);
    }
  }
  return Rcpp::wrap(res);
}

//' Calculate the Jaccard distance between two sets
//' @param x and y two vectors of features (1xn)
//' @param nFeat number of features in each sets
//' @return result distance between the two vectors
// [[Rcpp::export]]
double distJaccard(vec x, vec y, unsigned int nFeat) {
  
  uvec xSort = sort_index(x, "ascend");
  uvec ySort = sort_index(y, "ascend");
  
  xSort.print("xSort");
  ySort.print("ySort");
  
  uvec A = sort(xSort.head(nFeat));
  uvec B = sort(ySort.head(nFeat));
  uvec C = sort(xSort.tail(nFeat));
  uvec D = sort(ySort.tail(nFeat)); 
  
  A.print("A");
  B.print("B");
  C.print("C");
  D.print("D");
  
  vec Avec = conv_to<vec>::from(A);
  vec Bvec = conv_to<vec>::from(B);
  vec Cvec = conv_to<vec>::from(C);
  vec Dvec = conv_to<vec>::from(D);

  std::vector<int> ABunion;
  std::set_union(A.begin(), A.end(), 
                 B.begin(), B.end(), 
                 std::back_inserter(ABunion));
  std::vector<int> CDunion;
  std::set_union(C.begin(), C.end(), 
                 D.begin(), D.end(), 
                 std::back_inserter(CDunion));
  std::vector<int> ABinterect;
  std::set_intersection(A.begin(), A.end(),
                        B.begin(), B.end(),
                        std::back_inserter(ABinterect));
  std::vector<int> CDinterect;
  std::set_intersection(C.begin(), C.end(),
                        D.begin(), D.end(),
                        std::back_inserter(CDinterect));
  
  double dTop = 1 - ABinterect.size()/ABunion.size();
  double dBottom = 1 - CDinterect.size()/CDunion.size();
  
  double distance = (dTop + dBottom)/2.0;

  return distance;
}


/*** R
distJaccard(c(2, 1, 3, 4, 5, 9), c(9, 6, 1, 8, 7, 56), 5)
*/
