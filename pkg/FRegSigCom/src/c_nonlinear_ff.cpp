#include <RcppEigen.h>
#include <Rcpp.h>
#include "common_function.h"


// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;

/////////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::MatrixXd c_prod(Eigen::MatrixXd A_ini, Eigen::MatrixXd B_ini)
{
  Map<MatrixXd> A=Map<MatrixXd>(A_ini.data(), A_ini.rows(), A_ini.cols());
  Map<MatrixXd> B=Map<MatrixXd>(B_ini.data(), B_ini.rows(), B_ini.cols());
  return A*B;
}
///////////////////////////////////////////////////////
MatrixXd cal_R_inv(Map<MatrixXd> v, Map<MatrixXd> R_inv, int n_curves)
{
  int q=v.cols();
  MatrixXd out_ini(v.rows(), v.cols());
  Map<MatrixXd> out=Map<MatrixXd>(out_ini.data(), out_ini.rows(), out_ini.cols());
  
  int nbasis=R_inv.cols();
  
  for(int i=0; i<n_curves; i++)
  {
    out.block(i*nbasis,0, nbasis, q)=R_inv*v.block(i*nbasis,0, nbasis, q);
  }
  return out_ini;
}

/////////////////////////////////////////////////
MatrixXd cal_R_trans_inv(Map<MatrixXd> v, Map<MatrixXd> R_inv, int n_curves)
{
  int q=v.cols();
  MatrixXd out_ini(v.rows(), v.cols());
  Map<MatrixXd> out=Map<MatrixXd>(out_ini.data(), out_ini.rows(), out_ini.cols());
  
  int nbasis=R_inv.cols();
  
  for(int i=0; i<n_curves; i++)
  {
    out.block(i*nbasis,0, nbasis, q)=R_inv.transpose()*v.block(i*nbasis,0, nbasis, q);
  }
  return out_ini;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
////////// 
/////////////////////////////////////////////////

////////////////////////////////////////

MatrixXd get_cv_error_nonlinear(MatrixXd T_train, MatrixXd T_valid, MatrixXd Y_valid, List G_Y_Phi_train_list)
{
  int q=G_Y_Phi_train_list.size(), ncol=T_train.cols();
  MatrixXd t_train_mtx(T_train.rows(), ncol+1);
  t_train_mtx.col(0)=MatrixXd::Constant(T_train.rows(),1, 1/sqrt(T_train.rows()));
  t_train_mtx.rightCols(ncol)=T_train;
  MatrixXd t_valid_mtx(T_valid.rows(), ncol+1);
  t_valid_mtx.col(0)=MatrixXd::Constant(T_valid.rows(),1, 1/sqrt(T_train.rows()));
  t_valid_mtx.rightCols(ncol)=T_valid;
  MatrixXd error(q,ncol);
  for(int k=0; k<q; k++)
  {
    MatrixXd G_Y_Phi_train_ini=as<MatrixXd>(G_Y_Phi_train_list(k));
    Map<MatrixXd> G_Y_Phi_train=Map<MatrixXd>(G_Y_Phi_train_ini.data(), G_Y_Phi_train_ini.rows(), G_Y_Phi_train_ini.cols());
    MatrixXd E=Y_valid-t_valid_mtx.col(0)*(t_train_mtx.col(0).transpose()*G_Y_Phi_train);
    for(int ncomp=0;ncomp<ncol;ncomp++)
    {
      E.noalias()-=t_valid_mtx.col(ncomp+1)*(t_train_mtx.col(ncomp+1).transpose()*G_Y_Phi_train);
      error(k,ncomp)=E.squaredNorm()/Y_valid.cols();
    }
  }
  return error;
}

///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// [[Rcpp::export]]

List C_cv_nonlinear_ff(Eigen::MatrixXd G_ini, Eigen::MatrixXd Y, List x_params, List y_params, Eigen::MatrixXi all_folds, int upper_comp, double thresh)
{
  
  int n_curves=x_params(0),  K_cv=all_folds.cols(), nsample=Y.rows();
  
  
  Map<MatrixXd> G=Map<MatrixXd>(G_ini.data(), G_ini.rows(), G_ini.cols());
  
  /////////////////////////
  MatrixXd Y_cent=Y.rowwise()-Y.colwise().mean(), Pi=Y_cent*Y_cent.transpose()/nsample/Y.cols();
  Pi=(Pi+Pi.transpose())/2;
  
  MatrixXd A=as<MatrixXd>(as<List>(y_params(6))(0)), B=as<MatrixXd>(as<List>(y_params(6))(1));
  VectorXd kappa_set=as<VectorXd>(y_params(1));
  MatrixXd Phi=as<MatrixXd>(y_params(3));
  List G_Y_Phi_list(kappa_set.size());
  
  for(int i=0; i<kappa_set.size(); i++)
  {
    double kappa=kappa_set(i);
    G_Y_Phi_list(i)=Y*(Phi.transpose()*((A + kappa*B).ldlt().solve(Phi)))/Y.cols();
  }
  
  List Pi_train_list(K_cv), Y_train_list(K_cv), Y_valid_list(K_cv),G_Y_Phi_train_list_all(K_cv), G_Y_Phi_valid_list_all(K_cv);
  for(int i_cv=0; i_cv<K_cv; i_cv++)
  {
    VectorXi omit=all_folds.col(i_cv);
    int nsample_train=nsample-omit.sum();
    MatrixXd Y_valid=MatrixXd::Zero(omit.sum(), Y.cols());
    MatrixXd Y_train=MatrixXd::Zero(nsample_train, Y.cols());
    
    int i_train=0, i_valid=0;
    for(int i=0; i<nsample; i++)
    {
      if(omit(i)==1)
      {
        Y_valid.row(i_valid)=Y.row(i);
        i_valid=i_valid+1;
      }
      else
      {
        Y_train.row(i_train)=Y.row(i);
        i_train=i_train+1;
      }
    }
    MatrixXd Y_cent_train=Y_train.rowwise()-Y_train.colwise().mean();
    MatrixXd tmp_1=Y_cent_train*Y_cent_train.transpose()/nsample_train/Y.cols();
    Pi_train_list(i_cv)=(tmp_1+tmp_1.transpose())/2;
    Y_train_list(i_cv)=Y_train;
    Y_valid_list(i_cv)=Y_valid;
    
    List tmp_train_list(kappa_set.size()), tmp_valid_list(kappa_set.size());
    for(int ii=0; ii<kappa_set.size(); ii++)
    {
      MatrixXd G_Y_Phi=as<MatrixXd>(G_Y_Phi_list(ii));
      MatrixXd tmp_train=MatrixXd::Zero(nsample_train, G_Y_Phi.cols());
      MatrixXd tmp_valid=MatrixXd::Zero(omit.sum(), G_Y_Phi.cols());
      int i_train=0, i_valid=0;
      for(int i=0; i<nsample; i++)
      {
        if(omit(i)==1)
        {
          tmp_valid.row(i_valid)=G_Y_Phi.row(i);
          i_valid=i_valid+1;
        }
        else
        {
          tmp_train.row(i_train)=G_Y_Phi.row(i);
          i_train=i_train+1;
        }
      }
      tmp_train_list(ii)=tmp_train;
      tmp_valid_list(ii)=tmp_valid;
    }
    G_Y_Phi_train_list_all(i_cv)=tmp_train_list;
    G_Y_Phi_valid_list_all(i_cv)=tmp_valid_list;
  }
  
  VectorXd tau_set=as<VectorXd>(x_params(5));
  VectorXd lambda_set=as<VectorXd>(x_params(4));
  int  size_lambda=lambda_set.size();
  int opt_K,   opt_lambda_index;
  double opt_kappa, opt_lambda, min_error = 1e20, opt_tau;
  MatrixXd opt_z, opt_T;
  VectorXd opt_max_value;
  List opt_errors;
  MatrixXd J0=as<MatrixXd>(x_params(1)), J2=as<MatrixXd>(x_params(2));
  
  for(int i_tau=0; i_tau<tau_set.size(); i_tau++)
  {
    double tau=tau_set(i_tau);
    MatrixXd tmp=J0+tau*J2;
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    MatrixXd R_ini=lltOf.matrixU();
    
    Map<MatrixXd> R_inv=Map<MatrixXd>(R_ini.data(), R_ini.rows(), R_ini.cols());
    R_inv=R_inv.inverse();
    MatrixXd R_inv_tran_G=cal_R_trans_inv(G, R_inv, n_curves);
    
    VectorXi max_comp(size_lambda);
    List max_value_list(size_lambda), z_list(size_lambda), T_list(size_lambda);
    
    double lambda, tmp_value;
    
    ///determine the max components
    for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
    {
      
      lambda=lambda_set(i_lambda);
      int p=G.rows(), r=G.cols();
      MatrixXd M_ini(r+p, r);
      Map<MatrixXd> M=Map<MatrixXd>(M_ini.data(), M_ini.rows(), M_ini.cols());
      M.topRows(r)=sqrt(lambda)*MatrixXd::Identity(r, r);
      M.bottomRows(p)=-R_inv_tran_G;
      Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
      List fit=cal_comp_without_max(Pi_eig, M, upper_comp, thresh);
      
      MatrixXd beta=as<MatrixXd>(fit["beta"]);
      MatrixXd w=beta.topRows(nsample);
      MatrixXd v_ini=beta.bottomRows(beta.rows()-nsample);
      Map<MatrixXd> v=Map<MatrixXd>(v_ini.data(), v_ini.rows(), v_ini.cols());
      MatrixXd z=pow(lambda, -0.5)*cal_R_inv(v, R_inv, n_curves);
      max_comp(i_lambda)=beta.cols();
      max_value_list(i_lambda)=as<VectorXd>(fit["max_value"]);
      z_list(i_lambda)=z;
      T_list(i_lambda)=w;
    }
    
    ///conduct cross-validation
    
    VectorXd kappa_set=as<VectorXd>(y_params(1));
    List  errors(size_lambda);
    for(int i_cv=0; i_cv<K_cv; i_cv++)
    {
      
      MatrixXd Y_train=as<MatrixXd>(Y_train_list(i_cv));
      MatrixXd Y_valid=as<MatrixXd>(Y_valid_list(i_cv));
      MatrixXd Pi_train=as<MatrixXd>(Pi_train_list(i_cv));
      
      VectorXi omit=all_folds.col(i_cv);
      int nsample_train=nsample-omit.sum();
      MatrixXd R_inv_tran_G_valid(R_inv_tran_G.rows(), omit.sum());
      MatrixXd R_inv_tran_G_train(R_inv_tran_G.rows(), nsample_train);
      
      int i_train=0, i_valid=0;
      for(int i=0; i<nsample; i++)
      {
        if(omit(i)==1)
        {
          R_inv_tran_G_valid.col(i_valid)=R_inv_tran_G.col(i);
          i_valid=i_valid+1;
        }
        else
        {
          R_inv_tran_G_train.col(i_train)=R_inv_tran_G.col(i);
          i_train=i_train+1;
        }
      }
      
      R_inv_tran_G_valid=(R_inv_tran_G_valid.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
      R_inv_tran_G_train=(R_inv_tran_G_train.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
      
      List G_Y_Phi_train_list=G_Y_Phi_train_list_all(i_cv);
      for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
      {
        lambda=lambda_set(i_lambda);
        
        MatrixXd M(nsample_train+R_inv_tran_G_train.rows(), nsample_train);
        Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), M.cols());
        M_basis.topRows(nsample_train)=sqrt(lambda)*(MatrixXd::Identity(nsample_train, nsample_train));
        M_basis.bottomRows(R_inv_tran_G_train.rows())=-R_inv_tran_G_train;
        
        Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi_train.data(), Pi_train.rows(), Pi_train.cols());
        List fit=cal_comp_with_max(Pi_eig, M_basis, max_comp(i_lambda));
        MatrixXd beta=as<MatrixXd>(fit["beta"]);
        MatrixXd T_train=beta.topRows(nsample_train);
        
        MatrixXd T_valid=pow(lambda, -0.5)*R_inv_tran_G_valid.transpose()*beta.bottomRows(beta.rows()-nsample_train);
        
        for(int i=0; i<T_train.cols(); i++)
        {
          tmp_value=(T_train.col(i)).norm();
          T_train.col(i)=T_train.col(i)/tmp_value;
          T_valid.col(i)=T_valid.col(i)/tmp_value;
        }
        
        if(i_cv==0)
        {
          errors(i_lambda)=get_cv_error_nonlinear(T_train, T_valid, Y_valid, G_Y_Phi_train_list);
        }else
        {
          errors(i_lambda)=as<MatrixXd>(errors(i_lambda))+get_cv_error_nonlinear(T_train, T_valid, Y_valid, G_Y_Phi_train_list);
        }
      }
    }
    
    
    int tmp_opt_K=0, tmp_opt_lambda_index=0;
    double tmp_opt_kappa=0, tmp_opt_lambda=0, tmp_min_error=1e20;
    for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
    {
      MatrixXd::Index minRow, minCol;
      double tmp_error=(as<MatrixXd>(errors(i_lambda))).minCoeff(&minRow, &minCol);
      if(tmp_error<tmp_min_error)
      {
        tmp_min_error=tmp_error;
        tmp_opt_lambda_index=i_lambda;
         
        tmp_opt_kappa=kappa_set(minRow);
        tmp_opt_lambda=lambda_set(i_lambda);
        tmp_opt_K=minCol+1;
      }
      // Rcout << "errors=" <<" " << as<MatrixXd>(errors(i_lambda)) << std::endl;
    }
    if(min_error>tmp_min_error)
    {
      min_error=tmp_min_error;
      opt_lambda_index=tmp_opt_lambda_index;
      opt_tau=tau_set(i_tau);
           
      opt_K= tmp_opt_K;
      opt_kappa=tmp_opt_kappa;
      opt_lambda=tmp_opt_lambda;
      opt_errors=errors;
      MatrixXd opt_z_1=as<MatrixXd>(z_list(opt_lambda_index));
      opt_z=opt_z_1.leftCols(opt_K);
      MatrixXd opt_T_1=as<MatrixXd>(T_list(opt_lambda_index));
      opt_T=opt_T_1.leftCols(opt_K);
      VectorXd opt_max_value_1=as<VectorXd>(max_value_list(opt_lambda_index));
      opt_max_value=opt_max_value_1.head(opt_K);
    }
    
    //Rcout << "########################################" <<std::endl;
  }
   
  return List::create(_["min_error"]=min_error,  _["opt_K"]=opt_K, _["opt_lambda"]=opt_lambda, _["opt_tau"]=opt_tau,
                      _["opt_kappa"]=opt_kappa,  _["opt_T"]=opt_T, _["opt_z"]=opt_z );
}
