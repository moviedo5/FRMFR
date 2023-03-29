#include <RcppEigen.h>
#include <Rcpp.h>
#include "common_function.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;
////////////////////////////////////

///////////////////////////////////////////////////////////////////////
////////////////////////////////////////

MatrixXd get_cv_error_ffsig(MatrixXd T_train, MatrixXd T_valid, MatrixXd Y_valid, List G_Y_Phi_train_list)
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


/////////////////////////////////////////////////////
// [[Rcpp::export]]
List C_cv_ff_sig(List t_x, List X, Eigen::MatrixXd Y, List x_params, List y_params, List all_folds, int upper_comp, double thresh)
{
  int n_curves=t_x.size(), nsample=Y.rows(), K_cv=all_folds.size();
  int n_x_basis=as<MatrixXd>(as<List>(x_params(3))(0)).rows();
  int total_x_basis=n_curves*n_x_basis;
  List B_x_vals=x_params(3); // L*S
  
  MatrixXd G_ini(total_x_basis, nsample);//pL*n
  Map<MatrixXd> G=Map<MatrixXd>(G_ini.data(), total_x_basis, nsample);
  List G_mean_curves_list(n_curves), G_x_list(n_curves);
  for(int i=0; i<n_curves; i++)
  {
    MatrixXd tmp_mat=as<MatrixXd>(X(i)).transpose();  //S*n
    tmp_mat=as<MatrixXd>(B_x_vals(i))*tmp_mat/tmp_mat.rows(); //L*n  B*X'/S
    G_mean_curves_list(i)=tmp_mat.rowwise().mean();
    G.block(i*n_x_basis, 0, n_x_basis, nsample)=(tmp_mat.colwise()-tmp_mat.rowwise().mean())/sqrt(nsample); // cenceted BX'/S\sqrt{n}
  }
  
  ////////Pi=\int(Y-Ybar)(Y-Ybar)\trans dt /n
  MatrixXd Y_cent=Y.rowwise()-Y.colwise().mean(), Pi=Y_cent*Y_cent.transpose()/nsample/Y.cols();
  Pi=(Pi+Pi.transpose())/2;
  
  MatrixXd A=as<MatrixXd>(as<List>(y_params(6))(0)), B=as<MatrixXd>(as<List>(y_params(6))(1));
  VectorXd kappa_set=as<VectorXd>(y_params(1));
  MatrixXd Phi=as<MatrixXd>(y_params(3));
  List G_Y_Phi_list(kappa_set.size());
  for(int i=0; i<kappa_set.size(); i++)
  {
    double kappa=kappa_set(i);
    G_Y_Phi_list(i)=Y*(Phi.transpose()*((A + kappa*B).ldlt().solve(Phi)))/Y.cols();  //n*T, scalar predictors * GYPhi(i)=wi(t)
  }
  
  List Pi_train_list(K_cv), Y_train_list(K_cv), Y_valid_list(K_cv);
  List G_Y_Phi_train_list_all(K_cv), G_Y_Phi_valid_list_all(K_cv);
  for(int i_cv=0; i_cv<K_cv; i_cv++)
  {
    VectorXi omit;
    omit=as<VectorXi>(all_folds(i_cv)).array()-1;
    int nsample_train=nsample-omit.size();
    MatrixXd Y_valid=MatrixXd::Zero(omit.size(), Y.cols());
    MatrixXd Y_train=MatrixXd::Zero(nsample_train, Y.cols());
    VectorXi ind=VectorXi::Zero(nsample);
    for(int i=0; i<omit.size(); i++)
    {
      ind(omit(i))=1;
      Y_valid.row(i)=Y.row(omit(i));
    }
    int i_train=0;
    for(int i=0; i<nsample; i++)
    {
      if(ind(i)==0)
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
      VectorXi ind=VectorXi::Zero(nsample);
      for(int i=0; i<omit.size(); i++)
      {
        ind(omit(i))=1;
        tmp_valid.row(i)=G_Y_Phi.row(omit(i));
      }
      int i_train=0;
      for(int i=0; i<nsample; i++)
      {
        if(ind(i)==0)
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
  int opt_K,   opt_lambda_index;
  double opt_kappa, opt_lambda, min_error = 1e20, opt_tau;
  MatrixXd opt_z, opt_T;
  VectorXd opt_max_value;
  List opt_errors;
  MatrixXd J0=as<List>(x_params(4))(0), J2=as<List>(x_params(4))(1);
  
  for(int i_tau=0; i_tau<tau_set.size(); i_tau++)
  {
    double tau=tau_set(i_tau);
    MatrixXd tmp=J0+tau*J2;
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    MatrixXd R_ini=lltOf.matrixU();
    Map<MatrixXd> R_inv=Map<MatrixXd>(R_ini.data(), R_ini.rows(), R_ini.cols());
    R_inv=R_inv.inverse(); //Ls*Ls
    MatrixXd R_inv_tran_G(total_x_basis,nsample);  //G:pLs*n
    //Map<MatrixXd> R_inv_tran_G=Map<MatrixXd>(R_inv_tran_G_ini.data(), total_x_basis, nsample);
    for(int i=0; i<n_curves; i++){
      R_inv_tran_G.block(i*n_x_basis,0,n_x_basis,nsample)=R_inv.transpose()*G.block(i*n_x_basis,0,n_x_basis,nsample);
    }
    
    
    //MatrixXd R_inv_tran_G=cal_R_trans_inv(G,  x_params, n_curves, tau);
    VectorXd lambda_set=as<VectorXd>(x_params(1));
    int  size_lambda=lambda_set.size();
    VectorXi max_comp(size_lambda);
    List max_value_list(size_lambda), z_list(size_lambda), T_list(size_lambda);
    
    double lambda, tmp_value;
    for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
    {
      lambda=lambda_set(i_lambda);
      MatrixXd M(nsample+total_x_basis, nsample);
      //Rcout << "nsample=" << nsample << "  total_x_basis=" << total_x_basis << std::endl;
      Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), nsample);
      M_basis.topRows(nsample)=sqrt(lambda)*(MatrixXd::Identity(nsample, nsample));
      M_basis.bottomRows(total_x_basis)=-R_inv_tran_G;
      Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
      List fit=cal_comp_without_max(Pi_eig, M_basis, upper_comp, thresh);
      
      MatrixXd beta=as<MatrixXd>(fit["beta"]);
      MatrixXd w=beta.topRows(nsample);
      int tmp_col=beta.cols();
      MatrixXd v_ini=beta.bottomRows(beta.rows()-nsample), z_ini(v_ini.rows(),tmp_col);
      Map<MatrixXd> v=Map<MatrixXd>(v_ini.data(), v_ini.rows(), tmp_col);
      Map<MatrixXd> z=Map<MatrixXd>(z_ini.data(), z_ini.rows(), tmp_col);
      for(int i=0; i<n_curves; i++){
        z.block(i*n_x_basis,0,n_x_basis, tmp_col)=R_inv * v.block(i*n_x_basis,0,n_x_basis,tmp_col);
      }
      z=pow(lambda, -0.5)*z;
      max_comp(i_lambda)=tmp_col;
      max_value_list(i_lambda)=as<VectorXd>(fit["max_value"]);
      z_list(i_lambda)=z;
      T_list(i_lambda)=w;
    }
    
    ///conduct cross-validation
    int K_cv=all_folds.size();
    VectorXd kappa_set=as<VectorXd>(y_params(1));
    List  errors(size_lambda);
    for(int fold_ind=0; fold_ind<K_cv; fold_ind++)
    {
      MatrixXd Y_train=as<MatrixXd>(Y_train_list(fold_ind));
      MatrixXd Y_valid=as<MatrixXd>(Y_valid_list(fold_ind));
      MatrixXd Pi_train=as<MatrixXd>(Pi_train_list(fold_ind));
      VectorXi omit=as<VectorXi>(all_folds(fold_ind)).array()-1;
      int nsample_train=nsample-omit.size();
      MatrixXd R_inv_tran_G_valid(R_inv_tran_G.rows(), omit.size());
      MatrixXd R_inv_tran_G_train(R_inv_tran_G.rows(), nsample_train);
      VectorXi ind=VectorXi::Zero(nsample);
      for(int i=0; i<omit.size(); i++)
      {
        ind(omit(i))=1;
        R_inv_tran_G_valid.col(i)=R_inv_tran_G.col(omit(i));
      }
      int i_train=0;
      for(int i=0; i<nsample; i++)
      {
        if(ind(i)==0)
        {
          R_inv_tran_G_train.col(i_train)=R_inv_tran_G.col(i);
          i_train=i_train+1;
        }
      }
      R_inv_tran_G_valid=(R_inv_tran_G_valid.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
      R_inv_tran_G_train=(R_inv_tran_G_train.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
      List G_Y_Phi_train_list=G_Y_Phi_train_list_all(fold_ind);
      
      
      for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
      {
        lambda=lambda_set(i_lambda);
        MatrixXd M(nsample_train+total_x_basis, nsample_train);
        Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), M.cols());
        M_basis.topRows(nsample_train)=sqrt(lambda)*(MatrixXd::Identity(nsample_train, nsample_train));
        M_basis.bottomRows(R_inv_tran_G_train.rows())=-R_inv_tran_G_train;
        Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi_train.data(), Pi_train.rows(), Pi_train.cols());
        
        List fit=cal_comp_with_max(Pi_eig, M_basis, max_comp(i_lambda));
        MatrixXd beta=as<MatrixXd>(fit["beta"]);
        MatrixXd T_train=beta.topRows(nsample_train);
        MatrixXd T_valid=pow(lambda, -0.5)*R_inv_tran_G_valid.transpose()*beta.bottomRows(beta.rows()-nsample_train);//T=(BX\tran)\trans z, z=lambda^{-0.5}R^{-1}beta.bottom
        for(int i=0; i<T_train.cols(); i++)
        {
          tmp_value=(T_train.col(i)).norm();
          T_train.col(i)=T_train.col(i)/tmp_value;
          T_valid.col(i)=T_valid.col(i)/tmp_value;
        }
        
        if(fold_ind==0)
        {
          errors(i_lambda)=get_cv_error_ffsig(T_train, T_valid, Y_valid, G_Y_Phi_train_list);
        }else
        {
          errors(i_lambda)=as<MatrixXd>(errors(i_lambda))+get_cv_error_ffsig(T_train, T_valid, Y_valid, G_Y_Phi_train_list);
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
  }
  
  for(int i=0; i<opt_T.cols(); i++)
  {
    double tmp_value=opt_T.col(i).norm();
    opt_T.col(i)=opt_T.col(i)/tmp_value;
    opt_z.col(i)=opt_z.col(i)/tmp_value;
  }
  MatrixXd tmp_mat=Y*(Phi.transpose()*((A + opt_kappa*B).inverse()))/Y.cols();
  VectorXd mu_expan_coef=tmp_mat.colwise().mean();
  List beta_expan_coef(n_curves);
  for(int i=0; i<n_curves; i++)
  {
    MatrixXd tmp_z=opt_z.block(i*n_x_basis, 0, n_x_basis, opt_z.cols())/sqrt(nsample);
    VectorXd tmp_mean=as<VectorXd>(G_mean_curves_list(i));
    MatrixXd tmp_beta=tmp_z*(opt_T.transpose()*tmp_mat);
    mu_expan_coef=mu_expan_coef-tmp_beta.transpose()*tmp_mean;
    beta_expan_coef(i)=tmp_beta;
  }
   
  return List::create(_["mu.expan.coef"]=mu_expan_coef, _["beta.expan.coef"]=beta_expan_coef, _["min_error"]=min_error,  _["opt_K"]=opt_K,
                      _["opt_lambda"]=opt_lambda, _["opt_tau"]=opt_tau,   _["opt_kappa"]=opt_kappa);
}

////////////////////////////////////////////////////////////


// [[Rcpp::export]]
List C_cv_ff_sig_with_scalar(List t_x, List X, Eigen::MatrixXd Y, Eigen::MatrixXd Z, List x_params, List y_params, List all_folds, int upper_comp, double thresh)
{
  int n_curves=t_x.size(), nsample=Y.rows(), K_cv=all_folds.size(), n_scalar=Z.cols();
  int n_x_basis=as<MatrixXd>(as<List>(x_params(3))(0)).rows();
  int total_x_basis=n_curves*n_x_basis;
  List B_x_vals=x_params(3); // L*S
  int n_G_row=total_x_basis+n_scalar;
  MatrixXd G_ini(n_G_row, nsample);//pL*n
  Map<MatrixXd> G=Map<MatrixXd>(G_ini.data(), n_G_row, nsample);
  List G_mean_curves_list(n_curves), G_x_list(n_curves);
  
  for(int i=0; i<n_curves; i++)
  {
    MatrixXd tmp_mat=as<MatrixXd>(X(i)).transpose();  //S*n
    tmp_mat=as<MatrixXd>(B_x_vals(i))*tmp_mat/tmp_mat.rows(); //L*n  B*X'/S
    G_mean_curves_list(i)=tmp_mat.rowwise().mean();
    G.block(i*n_x_basis, 0, n_x_basis, nsample)=(tmp_mat.colwise()-tmp_mat.rowwise().mean())/sqrt(nsample); // cenceted BX'/S\sqrt{n}
  }
  G.bottomRows(n_scalar)=(Z.rowwise()-Z.colwise().mean()).transpose()/sqrt(nsample);
  VectorXd Z_mean=Z.colwise().mean();
  
  ////////Pi=\int(Y-Ybar)(Y-Ybar)\trans dt /n
  MatrixXd Y_cent=Y.rowwise()-Y.colwise().mean(), Pi=Y_cent*Y_cent.transpose()/nsample/Y.cols();
  Pi=(Pi+Pi.transpose())/2;
  
  MatrixXd A=as<MatrixXd>(as<List>(y_params(6))(0)), B=as<MatrixXd>(as<List>(y_params(6))(1));
  VectorXd kappa_set=as<VectorXd>(y_params(1));
  MatrixXd Phi=as<MatrixXd>(y_params(3));
  List G_Y_Phi_list(kappa_set.size());
  for(int i=0; i<kappa_set.size(); i++)
  {
    double kappa=kappa_set(i);
    G_Y_Phi_list(i)=Y*(Phi.transpose()*((A + kappa*B).ldlt().solve(Phi)))/Y.cols();  //n*T, scalar predictors * GYPhi(i)=wi(t)
  }
  
  List Pi_train_list(K_cv), Y_train_list(K_cv), Y_valid_list(K_cv);
  List G_Y_Phi_train_list_all(K_cv), G_Y_Phi_valid_list_all(K_cv);
  for(int i_cv=0; i_cv<K_cv; i_cv++)
  {
    VectorXi omit;
    omit=as<VectorXi>(all_folds(i_cv)).array()-1;
    int nsample_train=nsample-omit.size();
    MatrixXd Y_valid=MatrixXd::Zero(omit.size(), Y.cols());
    MatrixXd Y_train=MatrixXd::Zero(nsample_train, Y.cols());
    VectorXi ind=VectorXi::Zero(nsample);
    for(int i=0; i<omit.size(); i++)
    {
      ind(omit(i))=1;
      Y_valid.row(i)=Y.row(omit(i));
    }
    int i_train=0;
    for(int i=0; i<nsample; i++)
    {
      if(ind(i)==0)
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
      VectorXi ind=VectorXi::Zero(nsample);
      for(int i=0; i<omit.size(); i++)
      {
        ind(omit(i))=1;
        tmp_valid.row(i)=G_Y_Phi.row(omit(i));
      }
      int i_train=0;
      for(int i=0; i<nsample; i++)
      {
        if(ind(i)==0)
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
  VectorXd delta_set=as<VectorXd>(x_params(8));
  VectorXd lambda_set=as<VectorXd>(x_params(1));
  int  size_lambda=lambda_set.size();
  int opt_K,   opt_lambda_index;
  double opt_kappa, opt_lambda, opt_delta, min_error = 1e20, opt_tau;
  MatrixXd opt_z, opt_T;
  VectorXd opt_max_value;
  List opt_errors;
  MatrixXd J0=as<List>(x_params(4))(0), J2=as<List>(x_params(4))(1);
  
  for(int i_tau=0; i_tau<tau_set.size(); i_tau++)
  {
    double tau=tau_set(i_tau);
    MatrixXd tmp=J0+tau*J2;
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    MatrixXd R_ini=lltOf.matrixU();
    Map<MatrixXd> R_inv=Map<MatrixXd>(R_ini.data(), R_ini.rows(), R_ini.cols());
    R_inv=R_inv.inverse(); //Ls*Ls
    MatrixXd R_inv_tran_G(n_G_row, nsample);  //G:pLs*n
    
    for(int i=0; i<n_curves; i++)
    {
      R_inv_tran_G.block(i*n_x_basis,0,n_x_basis,nsample)=R_inv.transpose()*G.block(i*n_x_basis,0,n_x_basis,nsample);
    }
    
    for(int i_delta=0; i_delta<delta_set.size(); i_delta++)
    {
      double delta=delta_set(i_delta);
      
      VectorXi max_comp(size_lambda);
      List max_value_list(size_lambda), z_list(size_lambda), T_list(size_lambda);
      
      double lambda, tmp_value;
      for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
      {
        lambda=lambda_set(i_lambda);
        R_inv_tran_G.bottomRows(n_scalar)=G.bottomRows(n_scalar)/pow((delta/lambda), 0.5);
        MatrixXd M(nsample+n_G_row, nsample);
        //Rcout << "nsample=" << nsample << "  n_G_row=" << n_G_row << std::endl;
        Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), nsample);
        M_basis.topRows(nsample)=sqrt(lambda)*(MatrixXd::Identity(nsample, nsample));
        M_basis.bottomRows(n_G_row)=-R_inv_tran_G;
        Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
        List fit=cal_comp_without_max(Pi_eig, M_basis, upper_comp, thresh);
        
        MatrixXd beta=as<MatrixXd>(fit["beta"]);
        MatrixXd w=beta.topRows(nsample);
        int tmp_col=beta.cols();
        MatrixXd v_ini=beta.bottomRows(beta.rows()-nsample), z_ini(v_ini.rows(),tmp_col);
        Map<MatrixXd> v=Map<MatrixXd>(v_ini.data(), v_ini.rows(), tmp_col);
        Map<MatrixXd> z=Map<MatrixXd>(z_ini.data(), z_ini.rows(), tmp_col);
        for(int i=0; i<n_curves; i++)
        {
          z.block(i*n_x_basis,0,n_x_basis, tmp_col)=R_inv * v.block(i*n_x_basis,0,n_x_basis,tmp_col);
        }
        z.bottomRows(n_scalar)=v.bottomRows(n_scalar)/pow((delta/lambda), 0.5);
        z=pow(lambda, -0.5)*z;
        max_comp(i_lambda)=tmp_col;
        max_value_list(i_lambda)=as<VectorXd>(fit["max_value"]);
        z_list(i_lambda)=z;
        T_list(i_lambda)=w;
      }
      
      ///conduct cross-validation
      int K_cv=all_folds.size();
      VectorXd kappa_set=as<VectorXd>(y_params(1));
      List  errors(size_lambda);
      for(int fold_ind=0; fold_ind<K_cv; fold_ind++)
      {
        MatrixXd Y_train=as<MatrixXd>(Y_train_list(fold_ind));
        MatrixXd Y_valid=as<MatrixXd>(Y_valid_list(fold_ind));
        MatrixXd Pi_train=as<MatrixXd>(Pi_train_list(fold_ind));
        for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
        {
          lambda=lambda_set(i_lambda);
          R_inv_tran_G.bottomRows(n_scalar)=G.bottomRows(n_scalar)/pow((delta/lambda), 0.5);
          VectorXi omit=as<VectorXi>(all_folds(fold_ind)).array()-1;
          int nsample_train=nsample-omit.size();
          MatrixXd R_inv_tran_G_valid(R_inv_tran_G.rows(), omit.size());
          MatrixXd R_inv_tran_G_train(R_inv_tran_G.rows(), nsample_train);
          VectorXi ind=VectorXi::Zero(nsample);
          for(int i=0; i<omit.size(); i++)
          {
            ind(omit(i))=1;
            R_inv_tran_G_valid.col(i)=R_inv_tran_G.col(omit(i));
          }
          int i_train=0;
          for(int i=0; i<nsample; i++)
          {
            if(ind(i)==0)
            {
              R_inv_tran_G_train.col(i_train)=R_inv_tran_G.col(i);
              i_train=i_train+1;
            }
          }
          R_inv_tran_G_valid=(R_inv_tran_G_valid.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
          R_inv_tran_G_train=(R_inv_tran_G_train.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
          List G_Y_Phi_train_list=G_Y_Phi_train_list_all(fold_ind);
          
          
          MatrixXd M(nsample_train+n_G_row, nsample_train);
          Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), M.cols());
          M_basis.topRows(nsample_train)=sqrt(lambda)*(MatrixXd::Identity(nsample_train, nsample_train));
          M_basis.bottomRows(R_inv_tran_G_train.rows())=-R_inv_tran_G_train;
          Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi_train.data(), Pi_train.rows(), Pi_train.cols());
          
          List fit=cal_comp_with_max(Pi_eig, M_basis, max_comp(i_lambda));
          MatrixXd beta=as<MatrixXd>(fit["beta"]);
          MatrixXd T_train=beta.topRows(nsample_train);
          
          
          MatrixXd T_valid=pow(lambda, -0.5)*R_inv_tran_G_valid.transpose()*beta.bottomRows(beta.rows()-nsample_train);//T=(BX\tran)\trans z, z=lambda^{-0.5}R^{-1}beta.bottom
          for(int i=0; i<T_train.cols(); i++)
          {
            tmp_value=(T_train.col(i)).norm();
            T_train.col(i)=T_train.col(i)/tmp_value;
            T_valid.col(i)=T_valid.col(i)/tmp_value;
          }
          
          if(fold_ind==0)
          {
            errors(i_lambda)=get_cv_error_ffsig(T_train, T_valid, Y_valid, G_Y_Phi_train_list);
          }else
          {
            errors(i_lambda)=as<MatrixXd>(errors(i_lambda))+get_cv_error_ffsig(T_train, T_valid, Y_valid, G_Y_Phi_train_list);
          }
        }
      }
      
      int tmp_opt_K=0,  tmp_opt_lambda_index=0;
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
        
      }
      
      if(min_error>tmp_min_error)
      {
        min_error=tmp_min_error;
        opt_lambda_index=tmp_opt_lambda_index;
        opt_tau=tau;
        opt_delta=delta;
         
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
    }
  }
  
  for(int i=0; i<opt_T.cols(); i++)
  {
    double tmp_value=opt_T.col(i).norm();
    
    opt_T.col(i)=opt_T.col(i)/tmp_value;
    opt_z.col(i)=opt_z.col(i)/tmp_value;
  }
  MatrixXd tmp_mat=Y*(Phi.transpose()*((A + opt_kappa*B).inverse()))/Y.cols();
  VectorXd mu_expan_coef=tmp_mat.colwise().mean();
  MatrixXd W=opt_T.transpose()*tmp_mat;
  List beta_expan_coef(n_curves);
  for(int i=0; i<n_curves; i++)
  {
    MatrixXd tmp_z=opt_z.block(i*n_x_basis, 0, n_x_basis, opt_z.cols())/sqrt(nsample);
    VectorXd tmp_mean=as<VectorXd>(G_mean_curves_list(i));
    MatrixXd tmp_beta=tmp_z*W;
    mu_expan_coef=mu_expan_coef-tmp_beta.transpose()*tmp_mean;
    beta_expan_coef(i)=tmp_beta;
  }
  MatrixXd gamma_expan_coef=(opt_z.bottomRows(n_scalar)*W).transpose()/sqrt(nsample);
  mu_expan_coef=mu_expan_coef-gamma_expan_coef*Z_mean;
  
  
  return List::create(_["mu.expan.coef"]=mu_expan_coef, _["beta.expan.coef"]=beta_expan_coef, _["gamma.expan.coef"]=gamma_expan_coef,
                      _["min_error"]=min_error,  _["opt_K"]=opt_K, _["opt_lambda"]=opt_lambda,
                      _["opt_delta"]=opt_delta, _["opt_tau"]=opt_tau,   _["opt_kappa"]=opt_kappa);
}


