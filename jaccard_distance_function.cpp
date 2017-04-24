#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
using namespace std;

using namespace arma;
using namespace Rcpp;

//' Calculate the Jaccard distance between two sets
//' @param x and y two vectors of features (1xn)
//' @param nFeat number of features in each sets
//' @return result distance between the two vectors
double distJaccard(vec x, vec y, unsigned int nFeat) {
  
  uvec xSort = sort_index(x, "ascend");
  uvec ySort = sort_index(y, "ascend");
  
  uvec A = sort(xSort.head(nFeat));
  uvec B = sort(ySort.head(nFeat));
  uvec C = sort(xSort.tail(nFeat));
  uvec D = sort(ySort.tail(nFeat)); 

  std::vector<int> ABunion;
  std::set_union(A.begin(), A.end(), 
                 B.begin(), B.end(), 
                 std::back_inserter(ABunion));
  std::vector<int> CDunion;
  std::set_union(C.begin(), C.end(), 
                 D.begin(), D.end(), 
                 std::back_inserter(CDunion));
  std::vector<int> ABintersect;
  std::set_intersection(A.begin(), A.end(),
                        B.begin(), B.end(),
                        std::back_inserter(ABintersect));
  std::vector<int> CDintersect;
  std::set_intersection(C.begin(), C.end(),
                        D.begin(), D.end(),
                        std::back_inserter(CDintersect));
  
  double dTop = 1.0 - double(ABintersect.size())/ABunion.size();
  double dBottom = 1.0 - double(CDintersect.size())/CDunion.size();
  
  double distance = (dTop + dBottom)/2.0;

  return distance;
}

//' Calculate the Jaccard distances between all pairs of vector (different than themself)
//' @param A matrix of replicate mxn
//' @param nFeat number of features in each sets
//' @return vecotr of Jaccard distance
// [[Rcpp::export]]
NumericVector vecJaccardDistance(mat A, unsigned int nFeat){

  NumericVector distVec(A.n_rows*(A.n_rows-1)/2);
  
  unsigned int idx(0);
  
  for(unsigned int i = 0; i < A.n_rows-1; i++){
    vec x = (A.row(i)).t();
    for(unsigned int j = i+1; j < A.n_rows; j++){
      vec y = (A.row(j)).t();
      distVec[idx] = distJaccard(x, y, nFeat);
      idx++;
    }
    
  }
  return distVec;
}