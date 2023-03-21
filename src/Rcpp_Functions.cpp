// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppEigen.h>
#include <Eigen/Dense>


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}
