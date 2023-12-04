#include <RcppEigen.h>
#include <Eigen/Dense>

//' Multiply two matrices using Eigen library
//'
//' This function multiplies two matrices using the Eigen library, which provides
//' fast linear algebra operations.
//'
//' @param A a matrix
//' @param B a matrix
//' @return the product of A and B
//'
//' @importFrom Rcpp sourceCpp
//' @useDynLib SPACO
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  return A * B;
}
