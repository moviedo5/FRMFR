#include <RcppEigen.h>
#include <Rcpp.h>
#include "common_function.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;
/////////////////////////////
 
/////////////////////////////////////////////////////
// [[Rcpp::export]]
List c_cv_ff_spike(List t_x_list, List X, Eigen::MatrixXd Y, List x_params, List y_params, Eigen::MatrixXi all_folds, int upper_comp, double thresh)
{
  int ncurves=X.size(), max_level_x=x_params(3), K_cv=all_folds.cols(), nsample=Y.rows();

  ////////////////////////////intiatiation for X ///////////////
  MatrixXd wave_mat=as<MatrixXd>(x_params(4));
  VectorXd father_x=wave_mat.col(0), mother_x=wave_mat.col(1);
  VectorXi nbasis_levels_x=as<VectorXi>(x_params(6));
  int upper=as<int>(x_params(5));
  VectorXd ref_t=VectorXd::LinSpaced(father_x.size(), 0, upper);
  List G_list(ncurves);
  VectorXi ind_level_x=VectorXi::Zero(nbasis_levels_x(max_level_x));

  for(int i_level=1; i_level<=max_level_x; i_level++)
  {
    ind_level_x.segment(nbasis_levels_x(i_level-1), nbasis_levels_x(i_level)-nbasis_levels_x(i_level-1))=VectorXi::Constant(nbasis_levels_x(i_level)-nbasis_levels_x(i_level-1), i_level);
  }

  for(int i=0; i<ncurves;i++)
  {
    int ind=0;
    MatrixXd G_tmp=MatrixXd::Zero(nsample, nbasis_levels_x(max_level_x));
    ArrayXd t_x=as<ArrayXd>(t_x_list(i));
    MatrixXd X_tmp=as<MatrixXd>(X(i));
    Map<MatrixXd> Xmap=Map<MatrixXd>(X_tmp.data(), X_tmp.rows(), X_tmp.cols());
    int t_size=Xmap.cols();
    for(int k=1-upper; k<1;k++)
    {
      G_tmp.col(ind).noalias()=Xmap*eval_wave_basis(t_x-k, father_x, ref_t, upper)/t_size;
      ind=ind+1;
    }

    for(int i_level=0; i_level<=max_level_x; i_level++)
    {
      int N=pow_int(2, i_level);
      for(int k=1-upper; k<N;k++)
      {

        G_tmp.col(ind).noalias()=sqrt(N)*(Xmap*eval_wave_basis(N*t_x-k, mother_x, ref_t, upper))/t_size;
        ind=ind+1;
      }
    }
    G_list(i)=(G_tmp.rowwise()-G_tmp.colwise().mean())/sqrt(nsample);
  }

  VectorXd tau_set_x=as<VectorXd>(x_params(2));
  VectorXd alpha_set_x=as<VectorXd>(x_params(1));
  VectorXd lambda_set_x=as<VectorXd>(x_params(7));
  int total_x=0;
  for(int i_level=1; i_level<=max_level_x; i_level++)
  {
    for(int i_alpha=0; i_alpha<alpha_set_x.size(); i_alpha++)
    {
      for(int i_tau=0; i_tau<tau_set_x.size(); i_tau++)
      {
        for(int i_lambda=0; i_lambda<lambda_set_x.size(); i_lambda++)
        {
          total_x++;
        }
      }
    }
  }

  MatrixXd param_index_set_x(total_x, 4);
  int indic=0;
  for(int i_level=1; i_level<=max_level_x; i_level++)
  {
    for(int i_alpha=0; i_alpha<alpha_set_x.size(); i_alpha++)
    {
      double alpha=alpha_set_x(i_alpha);
      for(int i_tau=0; i_tau<tau_set_x.size(); i_tau++)
      {
        double tau=tau_set_x(i_tau);
        for(int i_lambda=0; i_lambda<lambda_set_x.size(); i_lambda++)
        {
          double lambda=lambda_set_x(i_lambda);

          param_index_set_x(indic,0)=i_level;
          param_index_set_x(indic,1)=alpha;
          param_index_set_x(indic,2)=tau;
          param_index_set_x(indic,3)=lambda;
          indic=indic+1;
        }
      }
    }
  }

  ////////////////////////////initiatiation for Y///////////////

  int max_level_y=y_params(3), ind=0;
  int upper_y=as<int>(y_params(5));
  MatrixXd wave_mat_y=as<MatrixXd>(y_params(4));
  VectorXd father_y=wave_mat_y.col(0), mother_y=wave_mat_y.col(1);
  VectorXi nbasis_levels_y=as<VectorXi>(y_params(6));
  MatrixXd G_Y=MatrixXd::Zero(nsample, nbasis_levels_y(max_level_y));
  ArrayXd t_y=as<ArrayXd>(y_params(0));
  int t_size=Y.cols();
  ref_t=VectorXd::LinSpaced(father_y.size(), 0, upper_y);
  MatrixXd Phi_ini(nbasis_levels_y(max_level_y), t_size);
  Map<MatrixXd> Phi=Map<MatrixXd>(Phi_ini.data(), Phi_ini.rows(), Phi_ini.cols());

  for(int k=1-upper_y; k<1;k++)
  {
    Phi.row(ind)=eval_wave_basis(t_y-k, father_y, ref_t, upper_y);
    ind=ind+1;
  }
  for(int i_level=0; i_level<=max_level_y; i_level++)
  {
    int N=pow_int(2, i_level);
    for(int k=1-upper_y; k<N;k++)
    {
      Phi.row(ind)=sqrt(N)*eval_wave_basis(N*t_y-k, mother_y, ref_t, upper_y);
      ind=ind+1;
    }
  }

  G_Y=Y*Phi.transpose()/t_size;

  VectorXi ind_level_y=VectorXi::Zero(nbasis_levels_y(max_level_y));

  for(int i_level=1; i_level<=max_level_y; i_level++)
  {
    ind_level_y.segment(nbasis_levels_y(i_level-1), nbasis_levels_y(i_level)-nbasis_levels_y(i_level-1))=VectorXi::Constant(nbasis_levels_y(i_level)-nbasis_levels_y(i_level-1), i_level);
  }

  MatrixXd Y_cent=Y.rowwise()-Y.colwise().mean(), Pi=Y_cent*Y_cent.transpose()/nsample/Y.cols();
  Pi=(Pi+Pi.transpose())/2;
  List Pi_train_list(K_cv), Y_train_list(K_cv),  Y_valid_list(K_cv),   G_Y_valid_list(K_cv);
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

    Y_cent=Y_train.rowwise()-Y_train.colwise().mean();
    MatrixXd tmp_1=Y_cent*Y_cent.transpose()/nsample_train/Y.cols();
    Pi_train_list(i_cv)=(tmp_1+tmp_1.transpose())/2;
    Y_train_list(i_cv)=Y_train;
    Y_valid_list(i_cv)=Y_valid;
  }

  VectorXd tau_set_y=as<VectorXd>(y_params(2));
  VectorXd alpha_set_y=as<VectorXd>(y_params(1));
  VectorXd lambda_set_y=as<VectorXd>(y_params(7));
  int total_y=0;
  for(int i_level=1; i_level<=max_level_y; i_level++)
  {
    for(int i_alpha=0; i_alpha<alpha_set_y.size(); i_alpha++)
    {
      for(int i_tau=0; i_tau<tau_set_y.size(); i_tau++)
      {
        for(int i_lambda=0; i_lambda<lambda_set_y.size(); i_lambda++)
        {
          total_y++;
         }
      }
    }
  }


  MatrixXd param_index_set_y(total_y, 4);
  indic=0;
  for(int i_level=1; i_level<=max_level_y; i_level++)
  {
    for(int i_alpha=0; i_alpha<alpha_set_y.size(); i_alpha++)
    {
      double alpha=alpha_set_y(i_alpha);
      for(int i_tau=0; i_tau<tau_set_y.size(); i_tau++)
      {
        double tau=tau_set_y(i_tau);
        for(int i_lambda=0; i_lambda<lambda_set_y.size(); i_lambda++)
        {
          double lambda=lambda_set_y(i_lambda);
          param_index_set_y(indic,0)=i_level;
          param_index_set_y(indic,1)=alpha;
          param_index_set_y(indic,2)=tau;
          param_index_set_y(indic,3)=lambda;
          indic=indic+1;
        }
      }
    }
  }


  /////////////////////////////// CV procedure/////////////////////////////
  /////////////////////////////// Calculate all components/////////////////////////////

  double min_cv_error=1e20;
  MatrixXd min_cv_mat;


  MatrixXd All_train_ini=MatrixXd::Zero(4*nsample, upper_comp*total_x);

  Map<MatrixXd> All_train=Map<MatrixXd>(All_train_ini.data(), All_train_ini.rows(), All_train_ini.cols());
  MatrixXd All_valid_ini=MatrixXd::Zero(nsample, upper_comp*total_x);

  Map<MatrixXd> All_valid=Map<MatrixXd>(All_valid_ini.data(), All_valid_ini.rows(), All_valid_ini.cols());
  VectorXi cv_size(K_cv), train_start(K_cv), valid_start(K_cv), All_comp(total_x), row_start(total_x);

  for(int i_cv=0; i_cv<K_cv; i_cv++)
  {
    cv_size(i_cv)=nsample-all_folds.col(i_cv).sum();
    if(i_cv==0)
    {
      train_start(i_cv)=0;
      valid_start(i_cv)=0;
    }
    else
    {
      train_start(i_cv)=cv_size.head(i_cv).sum();
      valid_start(i_cv)= i_cv*nsample-train_start(i_cv);
    }
  }
   MatrixXd errors;
   indic=0;
  List train_list(K_cv), valid_list(K_cv);
  row_start(0)=0;
   for(int i_level=1; i_level<=max_level_x; i_level++)
  {

    int col_size=nbasis_levels_x(i_level);
    MatrixXd G_tmp(ncurves*col_size, nsample);
    Map<MatrixXd> G=Map<MatrixXd>(G_tmp.data(), G_tmp.rows(), G_tmp.cols());
    VectorXi level_vec(ncurves*col_size);
    for(int j=0;j<ncurves;j++)
    {
      G.block(j*col_size, 0, col_size, nsample)=as<MatrixXd>(G_list(j)).leftCols(col_size).transpose();
      level_vec.segment(j*col_size, col_size)=ind_level_x.head(col_size);
    }
    int L=G.rows();

    List   G_valid_list(K_cv), G_train_list(K_cv);
     for(int i_cv=0; i_cv<K_cv; i_cv++)
    {
      VectorXi omit=all_folds.col(i_cv);

      int nsample_train=nsample-omit.sum();
      MatrixXd G_train(L, nsample_train), G_valid(L, omit.sum());
      int i_train=0, i_valid=0;

      for(int i=0; i<nsample; i++)
      {
        if(omit(i)==1)
        {
          G_valid.col(i_valid)=G.col(i);
          i_valid=i_valid+1;
        }
        else
        {
          G_train.col(i_train)=G.col(i);
          i_train=i_train+1;
        }
      }

      G_valid_list(i_cv)=(G_valid.colwise()-G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
      G_train_list(i_cv)=(G_train.colwise()-G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
    }

     for(int i_alpha=0; i_alpha<alpha_set_x.size(); i_alpha++)
    {
      double alpha=alpha_set_x(i_alpha);
      for(int i_tau=0; i_tau<tau_set_x.size(); i_tau++)
      {
        double tau=tau_set_x(i_tau);
        for(int i_lambda=0; i_lambda<lambda_set_x.size(); i_lambda++)
        {
          double lambda=lambda_set_x(i_lambda);

          VectorXd K_inv_half(L);
          for(int i=0; i<L; i++)
          {
            K_inv_half(i)= pow(lambda*exp(tau*(level_vec(i)-i_level))*pow(2, 2*alpha*(level_vec(i))), -0.5);
          }
           int p=G.rows(), r=G.cols();
          MatrixXd M_ini(r+p, r);
          Map<MatrixXd> M=Map<MatrixXd>(M_ini.data(), M_ini.rows(), M_ini.cols());
          M.topRows(r)=MatrixXd::Identity(r, r);
          for(int i=0; i<p; i++)
          {
            M.row(r+i)=-G.row(i)*K_inv_half(i);
          }
          Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
          List fit=cal_comp_without_max(Pi_eig, M, upper_comp, thresh);

          MatrixXd beta=as<MatrixXd>(fit["beta"]);
          MatrixXd w=beta.topRows(nsample);
          int tmp_col=beta.cols();
          MatrixXd v_ini=beta.bottomRows(p), z_ini(p, r);
          Map<MatrixXd> v=Map<MatrixXd>(v_ini.data(), v_ini.rows(), tmp_col);
          Map<MatrixXd> z=Map<MatrixXd>(z_ini.data(), z_ini.rows(), tmp_col);
          for(int j=0;j<v.cols();j++)
          {
            z.col(j)=v.col(j).array()*K_inv_half.array();
          }

          int  max_comp=tmp_col;

          for(int i_cv=0; i_cv<K_cv; i_cv++)
          {
            MatrixXd G_train=as<MatrixXd>(G_train_list(i_cv));
            MatrixXd G_valid=as<MatrixXd>(G_valid_list(i_cv));
            MatrixXd Y_train=as<MatrixXd>(Y_train_list(i_cv));
            MatrixXd Y_valid=as<MatrixXd>(Y_valid_list(i_cv));
            MatrixXd Pi_train=as<MatrixXd>(Pi_train_list(i_cv));

            p=G_train.rows(), r=G_train.cols();
            M_ini.resize(r+p, r);
            new (&M) Map<MatrixXd>(M_ini.data(), M_ini.rows(), M_ini.cols());
            M.topRows(r)=MatrixXd::Identity(r, r);
            for(int i=0; i<p; i++)
            {
              M.row(r+i)=-G_train.row(i)*K_inv_half(i);;
            }
            Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi_train.data(), Pi_train.rows(), Pi_train.cols());

            List fit=cal_comp_with_max(Pi_eig, M, max_comp);

            MatrixXd beta=as<MatrixXd>(fit["beta"]);
            v=beta.bottomRows(beta.rows()-r);
            MatrixXd z_train(v.rows(), v.cols());
            for(int j=0;j<v.cols();j++)
            {
              z_train.col(j)=v.col(j).array()*K_inv_half.array();
            }
            MatrixXd T_train=beta.topRows(r);
            MatrixXd T_valid=G_valid.transpose()*z_train;
            for(int i=0; i<T_train.cols(); i++)
            {
              double tmp_value=(T_train.col(i)).norm();
              T_train.col(i).noalias()=T_train.col(i)/tmp_value;
              T_valid.col(i).noalias()=T_valid.col(i)/tmp_value;
            }

            int ncomp=T_train.cols();
            All_train.block(train_start(i_cv), row_start(indic), T_train.rows(), ncomp)=T_train;
            All_valid.block(valid_start(i_cv), row_start(indic), T_valid.rows(), ncomp)=T_valid;

            All_comp(indic)=ncomp;
            if(indic<total_x-1)
            {
              row_start(indic+1)=row_start(indic)+ncomp;
            }

          }
          indic=indic+1;

        }
      }
    }
  }

  /////////////////////////////// calculate all cv errors /////////////////////////////


  int opt_i_x, opt_i_y, opt_comp;
  indic=0;
  MatrixXd cv_errors=MatrixXd::Constant(upper_comp, total_x, 1e30);
  int q=Y.cols();
  for(int i_level=1; i_level<=max_level_y; i_level++)
  {
    int L=nbasis_levels_y(i_level);
    MatrixXd D_ini=G_Y.leftCols(L);
    Map<MatrixXd> D=Map<MatrixXd>(D_ini.data(), D_ini.rows(), D_ini.cols());
    MatrixXd B_ini(D.rows(), L);
    Map<MatrixXd> B=Map<MatrixXd>(B_ini.data(), B_ini.rows(), B_ini.cols());
    VectorXi level_vec_y=ind_level_y.head(L);

    for(int i_alpha=0; i_alpha<alpha_set_y.size(); i_alpha++)
    {
      double alpha=alpha_set_y(i_alpha);
      for(int i_tau=0; i_tau<tau_set_y.size(); i_tau++)
      {
        double tau=tau_set_y(i_tau);
        for(int i_lambda=0; i_lambda<lambda_set_y.size(); i_lambda++)
        {
          double lambda=lambda_set_y(i_lambda);
          for(int i=0; i<L; i++)
          {
            double tmp= 1+lambda*exp(tau*(level_vec_y(i)-i_level))*pow(2, 2*alpha*(level_vec_y(i)));

            B.col(i).noalias()=D.col(i)/tmp;
          }
          MatrixXd tmp_ini=B*Phi.topRows(L);

          Map<MatrixXd> tmp=Map<MatrixXd>(tmp_ini.data(), tmp_ini.rows(), tmp_ini.cols());
          for(int i_cv=0; i_cv<K_cv; i_cv++)
          {
            VectorXi omit=all_folds.col(i_cv);

            MatrixXd tmp_train_ini(cv_size(i_cv), tmp.cols());
            Map<MatrixXd>  tmp_train=Map<MatrixXd>(tmp_train_ini.data(), tmp_train_ini.rows(), tmp_train_ini.cols());
            int i_train=0;
            for(int i=0; i<nsample; i++)
            {
              if(omit(i)==0)
              {
                tmp_train.row(i_train).noalias()=tmp.row(i);
                i_train=i_train+1;
              }
            }

            MatrixXd Y_valid=as<MatrixXd>(Y_valid_list(i_cv));
            MatrixXd E0=Y_valid.rowwise()-tmp_train.colwise().mean();
            MatrixXd E_ini(E0.rows(), E0.cols());
            Map<MatrixXd> E=Map<MatrixXd>(E_ini.data(), E_ini.rows(), E_ini.cols());
            MatrixXd All_coef_ini=tmp_train.transpose()*All_train.block(train_start(i_cv), 0, cv_size(i_cv), All_comp.sum());
            Map<MatrixXd> All_coef=Map<MatrixXd>(All_coef_ini.data(), All_coef_ini.rows(), All_coef_ini.cols());
            for(int jj=0; jj<total_x; jj++)
            {
              E=E0;
              for(int ii=0; ii<All_comp(jj); ii++)
              {
                E.noalias()-=All_valid.block(valid_start(i_cv), row_start(jj)+ii, cv_size(i_cv), 1)*All_coef.col(row_start(jj)+ii).transpose();

                if(i_cv==0)
                {
                  cv_errors(ii, jj)=E.squaredNorm()/q;
                }
                else
                {
                  cv_errors(ii, jj)+=E.squaredNorm()/q;
                }
              }
            }
          }
          MatrixXd::Index i_x_params, i_comp;
          double tmp_min_error=cv_errors.minCoeff(&i_comp, &i_x_params);
          if(tmp_min_error<min_cv_error)
          {
            min_cv_error=tmp_min_error;
            opt_i_x=i_x_params;
            opt_i_y=indic;
            opt_comp=i_comp+1;
            min_cv_mat=cv_errors;
          }
          indic=indic+1;
          cv_errors=MatrixXd::Constant(upper_comp, total_x, 1e30);

        }
      }
    }
  }

  ///////////////////////////////////////// identify the optimal parameters and calculate the coefficient functions/////////////


  int opt_i_level_y=param_index_set_y(opt_i_y,0);
  double opt_alpha_y=param_index_set_y(opt_i_y,1);
  double opt_tau_y=param_index_set_y(opt_i_y,2);
  double opt_lambda_y=param_index_set_y(opt_i_y,3);

  int L_y=nbasis_levels_y(opt_i_level_y);
  MatrixXd D_ini=G_Y.leftCols(L_y);
  Map<MatrixXd> D=Map<MatrixXd>(D_ini.data(), D_ini.rows(), D_ini.cols());
  MatrixXd B_ini(D.rows(), L_y);
  Map<MatrixXd> B=Map<MatrixXd>(B_ini.data(), B_ini.rows(), B_ini.cols());
  VectorXi level_vec_y=ind_level_y.head(L_y);
  for(int i=0; i<L_y; i++)
  {
    double tmp= 1+opt_lambda_y*exp(opt_tau_y*(level_vec_y(i)-opt_i_level_y))*pow(2, 2*opt_alpha_y*(level_vec_y(i)));
    B.col(i).noalias()=D.col(i)/tmp;
  }



  int opt_i_level_x=param_index_set_x(opt_i_x,0);
  double opt_alpha_x=param_index_set_x(opt_i_x,1);
  double opt_tau_x=param_index_set_x(opt_i_x,2);
  double opt_lambda_x=param_index_set_x(opt_i_x,3);

  int col_size=nbasis_levels_x(opt_i_level_x);
  MatrixXd G_tmp(ncurves*col_size, nsample);
  Map<MatrixXd> G=Map<MatrixXd>(G_tmp.data(), G_tmp.rows(), G_tmp.cols());
  VectorXi level_vec(ncurves*col_size);
  for(int j=0;j<ncurves;j++)
  {
    G.block(j*col_size, 0, col_size, nsample)=as<MatrixXd>(G_list(j)).leftCols(col_size).transpose();
    level_vec.segment(j*col_size, col_size)=ind_level_x.head(col_size);
  }


  int L=G.rows();


  VectorXd K_inv_half(L);
  for(int i=0; i<L; i++)
  {

      K_inv_half(i)=pow(opt_lambda_x*exp(opt_tau_x*(level_vec(i)-opt_i_level_x))*pow(2, 2*opt_alpha_x*(level_vec(i))), -0.5);

  }
  int p=G.rows(), r=G.cols();

  MatrixXd M_ini(r+p, r);
  Map<MatrixXd> M=Map<MatrixXd>(M_ini.data(), M_ini.rows(), M_ini.cols());
  M.topRows(r)=MatrixXd::Identity(r, r);
  for(int i=0; i<p; i++)
  {
    M.row(r+i)=-G.row(i)*K_inv_half(i);
  }
  Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
  List fit=cal_comp_with_max(Pi_eig, M, opt_comp);

  MatrixXd beta=as<MatrixXd>(fit["beta"]);
  MatrixXd T_train=beta.topRows(r);
  MatrixXd v=beta.bottomRows(beta.rows()-r);
  MatrixXd Z(v.rows(), v.cols());
  for(int j=0;j<v.cols();j++)
  {
    Z.col(j)=v.col(j).array()*K_inv_half.array();
  }
  for(int i=0; i<T_train.cols(); i++)
  {
    double tmp_value=(T_train.col(i)).norm();
    T_train.col(i).noalias()=T_train.col(i)/tmp_value;
    Z.col(i).noalias()=Z.col(i)/tmp_value;
  }
  MatrixXd W=T_train.transpose()*(B*Phi.topRows(L));
  List Beta_list(ncurves);
  VectorXd mu=(B*Phi.topRows(L)).colwise().mean();

  for(int i=0; i<ncurves;i++)
  {
    int ind=0;
    ArrayXd t_x=as<ArrayXd>(t_x_list(i));
    MatrixXd X_tmp=as<MatrixXd>(X(i));
    int t_size=X_tmp.cols();
    MatrixXd Phi_x_ini(col_size, t_size);
    Map<MatrixXd> Phi_x=Map<MatrixXd>(Phi_x_ini.data(), Phi_x_ini.rows(), Phi_x_ini.cols());
    int upper_x=as<int>(x_params(5));
    for(int k=1-upper_x; k<1;k++)
    {
      Phi_x.row(ind)=eval_wave_basis(t_x-k, father_x, ref_t, upper_x);
      ind=ind+1;
    }
    for(int i_level=0; i_level<=opt_i_level_x; i_level++)
    {
      int N=pow_int(2, i_level);
      for(int k=1-upper_x; k<N;k++)
      {
        Phi_x.row(ind)=sqrt(N)*eval_wave_basis(N*t_x-k, mother_x, ref_t, upper_x);
        ind=ind+1;
      }
    }

    beta=(Phi_x.transpose()*Z.block(i*col_size, 0, col_size, Z.cols()))*W/sqrt(nsample);
    VectorXd tmp_vec=X_tmp.colwise().mean();
    mu=mu-beta.transpose()*tmp_vec/t_size;
    Beta_list(i)=beta;
  }

  return List::create(_["min_cv_mat"]=min_cv_mat, _["mu"]=mu, _["Beta_list"]=Beta_list, _["opt_comp"]=opt_comp, _["opt_i_level_x"]=opt_i_level_x,
                      _["opt_alpha_x"]=opt_alpha_x, _["opt_tau_x"]=opt_tau_x,  _["opt_lambda_x"]=opt_lambda_x, _["opt_i_level_y"]=opt_i_level_y,
                      _["opt_alpha_y"]=opt_alpha_y, _["opt_tau_y"]=opt_tau_y,  _["opt_lambda_y"]=opt_lambda_y, _["param_index_set_x"]=param_index_set_x,
                      _["param_index_set_y"]=param_index_set_y, _["opt_i_x"]=opt_i_x,    _["opt_i_y"]=opt_i_y);
}


