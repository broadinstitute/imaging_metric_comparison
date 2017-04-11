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

// [[Rcpp::export]]
NumericVector CE_entropy(NumericMatrix A){
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
  /*mat test(3,3);
  test << 1 << 1 << 1 << endr
       << 2 << 2 << 2 << endr
       << 3 << 3 << 3 << endr;
  //<< 4 << 4 << 4 << endr;
  
  mat testt  = test * test.t();
  
  testt.print("Testt: ");
  
  vec s = svd(test);
  
  s.print("s: ");
  
  vec v = s/sum(s); 
  v.print("v: ");
  uvec idxZer0 = find(v < 0.5);
  idxZer0.print("v: ");
  int l = idxZer0.size();
  cout << l << endl;
  
  
  vec t(4);
  t << 2 << 3 << 4 << 5;
  arma::vec vNoZero = t.elem(find(t >= 0));
  double E = -sum(vNoZero % log(vNoZero))/log(vNoZero.size()); 
  cout << E << endl;*/

/*
A <- profiles %>% select(one_of(variables)) %>% as.matrix() %>% t(.)
*/
