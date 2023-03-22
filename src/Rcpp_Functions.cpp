//' Multiply two matrices using Eigen library
//'
//' This function multiplies two matrices using the Eigen library, which provides
//' fast linear algebra operations.
//'
//' @param A a matrix
//' @param B a matrix
//' @return the product of A and B
//' @export
//'
//' @examples
//' A <- matrix(1:6, nrow=2)
//' B <- matrix(7:12, nrow=2)
//' eigenMapMatMult(A, B)
//'
//' @importFrom Rcpp sourceCpp
//' @useDynLib SPACO


// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Eigen/Dense>


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}
