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


/*** R
#distJaccard(c(2, 1, 3, 4, 5, 9), c(9, 6, 1, 8, 7, 56), 5)

B <- matrix(c(2, 1, 3, 4, 5, 9, 9, 6, 1, 8, 7, 56, 59, 3, 2, 5, 1, 98, 38, 29, 2, 5, 23, 6), 
           nrow=4, 
           ncol=6) 

temp <- vecJaccardDistance(B, 3)

*/
