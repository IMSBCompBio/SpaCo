#include <RcppEigen.h>
#include <Eigen/Dense>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
//' Multiply two matrices using Eigen library
//'
//' This function multiplies two matrices using the Eigen library, which provides
//' fast linear algebra operations.
//'
//' @param A a matrix
//' @param B a matrix
//' @return the product of A and B
//' @useDynLib SPACO
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd eigenMapMatMult(const Eigen::MatrixXd& A,
                                const Eigen::MatrixXd& B){
  return A * B;
}
