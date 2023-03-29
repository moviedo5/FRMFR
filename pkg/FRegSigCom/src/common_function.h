#ifndef COMMON_FUNCTION_H
#define COMMON_FUNCTION_H

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace std;
using namespace Eigen;

void  power_eig(Eigen::Map<Eigen::MatrixXd> Pi_eig, Eigen::Map<Eigen::MatrixXd> M, Eigen::Map<Eigen::MatrixXd> MM_inv, Eigen::MatrixXd J, double & max_value, Eigen::Map<Eigen::VectorXd>& beta_current);
void  power_eig_0(Eigen::Map<Eigen::MatrixXd> Pi_eig, Eigen::Map<Eigen::MatrixXd> M, Eigen::Map<Eigen::MatrixXd> MM_inv, double & max_value, Eigen::Map<Eigen::VectorXd>& beta_current);
List cal_comp_without_max(Map<MatrixXd> Pi_eig, Map<MatrixXd> M, int upper_comp, double thresh);
List cal_comp_with_max(Map<MatrixXd> Pi_eig, Map<MatrixXd> M, int max_comp);
VectorXi extract(VectorXi x, VectorXi ind);
List basic_max(VectorXd a, double lambda);
VectorXd eval_wave_basis(ArrayXd t, VectorXd ref_vec, VectorXd ref_t, int upper);
int pow_int(int base, int power);
Eigen::MatrixXd find_orth_basis(Eigen::MatrixXd X);

#endif
