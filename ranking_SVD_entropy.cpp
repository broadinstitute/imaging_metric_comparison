#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
using namespace std;

using namespace arma;
using namespace Rcpp;

//' Calculate the entropy of a matrix based on SVD
//' @param A dataset mxn (m features and n observations) 
//' @return entropy
double calculate_entropy(mat A) {
	
  // calculate the svd
	vec s = svd(A);
  // calculate eigenvalue
	s = pow(s, 2); 
	// normalize relative values
	vec v = s/sum(s);
	
	// calculate the entropy for values bigger than 0 (to avoid error)
	arma::vec vNoZero = v.elem(find(v > 0));
	double E = -sum(vNoZero % log(vNoZero))/log(vNoZero.size()); 
	
	return E;
}

/* SR method
Simple Ranking:
 select mc features according to the highest ranking order of their CE values
 */
//' Calculate the entropy contribution of a matrix based on SVD
//' @param A dataset mxn (m features and n observations) 
//' @return entropy contribution vector
// [[Rcpp::export]]
NumericVector CE_entropy_SR(NumericMatrix A){
  // convert into matrix (armadillo)
  mat Amat(A.begin(), A.nrow(), A.ncol(), false);
  // total entropy
	double E = calculate_entropy(Amat);
  // Contributio vector to the entropy by a leave-one-out comparison
  NumericVector CE(A.nrow());
  // for each feature calculate the contribution to the entropy by a leave-one-out comparison
	for(unsigned int i = 0; i < A.nrow(); i++){
	  mat Ai = Amat;
	  // remove the row i
	  Ai.shed_row( i ) ;
	  //cout << "Ai size: " << Ai.n_cols << "columns and " << Ai.n_rows << "rows" << endl;
	  double Ei = calculate_entropy(Ai);
	  CE[i] = E - Ei;
	}
	
	return CE;
}

/* FS2 method
Forward selection
*/
//' Calculate the entropy of a matrix based on SVD with FS method
//' @param A dataset mxn (m features and n observations) 
//' @return index of features that were selected
// [[Rcpp::export]]
NumericVector CE_entropy_FS2(NumericMatrix A, unsigned int mc){
  // convert into matrix (armadillo)
  mat Amat(A.begin(), A.nrow(), A.ncol(), false);
  
  // vector of index of the features
  vec idxFeat(linspace<vec>(1, Amat.n_rows, Amat.n_rows));
  
  // best features index
  NumericVector idxBest(mc);
  
  for(unsigned int j = 0; j < mc; j++){
    // total entropy
    double E = calculate_entropy(Amat);
    
    // Contribution vector to the entropy by a leave-one-out comparison
    vec CE(Amat.n_rows);
    
    // for each feature calculate the contribution to the entropy by a leave-one-out comparison
    for(unsigned int i = 0; i < Amat.n_rows; i++){
      mat Ai = Amat;
      // remove row i
      Ai.shed_row(i);
      double Ei = calculate_entropy(Ai);
      CE[i] = E - Ei;
    }
    // find the index of the highest entropy contribution
    int idx = CE.index_max();
    idxBest[j] = idxFeat[idx];
    
    // remove the best feature
    Amat.shed_row(idx);
    idxFeat.shed_row(idx);
  }
  cout << "Best features: " << idxBest << endl;
  
  return idxBest;
}

/* BE method
 Backward Elimination
 Super expensive method...
*/
//' Calculate the entropy of a matrix based on SVD with BE method
//' @param A dataset mxn (m features and n observations) 
//' @return index of features that were selected
// [[Rcpp::export]]
NumericVector CE_entropy_BE(NumericMatrix A, unsigned int mc){
  // convert into matrix (armadillo)
  mat Amat(A.begin(), A.nrow(), A.ncol(), false);
  
  // vector of index of the features
  vec idxFeat(linspace<vec>(1, Amat.n_rows, Amat.n_rows));
  
  // best features index
  //NumericVector idxBest(mc);
  
  for(unsigned int j = 0; j < Amat.n_rows-mc; j++){ // TODO: test that this is correct
    // total entropy
    double E = calculate_entropy(Amat);
    
    // Contribution vector to the entropy by a leave-one-out comparison
    vec CE(Amat.n_rows);
    
    // for each feature calculate the contribution to the entropy by a leave-one-out comparison
    for(unsigned int i = 0; i < Amat.n_rows; i++){
      mat Ai = Amat;
      // remove row i
      Ai.shed_row(i);
      double Ei = calculate_entropy(Ai);
      CE[i] = E - Ei;
    }
    // find the index of the highest entropy contribution
    int idx = CE.index_min();
    //idxBest[j] = idxFeat[idx];
    
    // remove the best feature
    Amat.shed_row(idx);
    idxFeat.shed_row(idx);
  }
  cout << "Best features: " << idxFeat << endl;
  
  return wrap(idxFeat);
}
