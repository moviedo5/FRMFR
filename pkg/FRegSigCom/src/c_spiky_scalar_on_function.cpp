#include <RcppEigen.h>
#include <Rcpp.h>
#include "common_function.h"


// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;


/////////////////////////////////////////
List QR_full(Map<MatrixXd> B)
{
  int L=B.rows(), n=B.cols();
  double tol=1e-12, tol2=1e-12;
  MatrixXd C_ini=B;
  Map<MatrixXd> C=Map<MatrixXd>(C_ini.data(), L, n);
  MatrixXd D_ini=MatrixXd::Zero(L,n);
  Map<MatrixXd> D=Map<MatrixXd>(D_ini.data(), L, n);
  double max_norm_all=C.colwise().norm().maxCoeff();
  int rank=0, i_switch=-1;
  MatrixXi switch_mat(n, 2);
  
  for(int j=0;j<n;j++)
  {
    MatrixXd::Index maxindex;
    double max_norm=C.bottomRightCorner(L-j,n-j).colwise().norm().maxCoeff(&maxindex);
    if(max_norm>tol2*max_norm_all)
    {
      if(maxindex!=0)
      {
        VectorXd tmp=C.col(j);
        C.col(j)=C.col(j+maxindex);
        C.col(j+maxindex)=tmp;
        i_switch++;
        switch_mat(i_switch,0)=j;
        switch_mat(i_switch,1)=j+maxindex;
      }
      VectorXd v=C.block(j,j,L-j,1);
      if(v.tail(L-j-1).norm()>tol)
      {
        if(v(0)>0)
        {
          v(0)=v(0)+max_norm;
        }
        else
        {
          v(0)=v(0)-max_norm;
        }
        v.noalias()=v/v.norm();
        D.block(j,j,L-j,1)=v;
        C.block(j,j,L-j,n-j).noalias()-=2*v*(v.transpose()*C.block(j,j,L-j,n-j));
      }
      rank++;
    }
    else
    {
      break;
    }
  }
  
  
  for(int j=i_switch;j>=0;j--)
  {
    int j1=switch_mat(j,0), j2=switch_mat(j,1);
    VectorXd tmp=C.col(j1);
    C.col(j1)=C.col(j2);
    C.col(j2)=tmp;
  }
  MatrixXd Q_ini=MatrixXd::Identity(L, n);
  Map<MatrixXd> Q=Map<MatrixXd>(Q_ini.data(), L, n);
  for(int i=0; i<n;i++)
  {
    Q.noalias()-=2*D.col(n-i-1)*(D.col(n-i-1).transpose()*Q);
  }
  return  List::create(_["R"]=C_ini.topRows(rank), _["Q"]=Q_ini.leftCols(rank));
}




/////////////////////////////////////////
MatrixXd QR_house(Map<MatrixXd> B)
{
  int L=B.rows(), n=B.cols();
  double tol=1e-12, tol2=1e-12;
  MatrixXd C_ini=B;
  Map<MatrixXd> C=Map<MatrixXd>(C_ini.data(), L, n);
  double max_norm_all=C.colwise().norm().maxCoeff();
  int rank=0, i_switch=-1;
  MatrixXi switch_mat(n, 2);
  
  for(int j=0;j<n;j++)
  {
    MatrixXd::Index maxindex;
    double max_norm=C.bottomRightCorner(L-j,n-j).colwise().norm().maxCoeff(&maxindex);
    if(max_norm>tol2*max_norm_all)
    {
      if(maxindex!=0)
      {
        VectorXd tmp=C.col(j);
        C.col(j)=C.col(j+maxindex);
        C.col(j+maxindex)=tmp;
        i_switch++;
        switch_mat(i_switch,0)=j;
        switch_mat(i_switch,1)=j+maxindex;
      }
      VectorXd v=C.block(j,j,L-j,1);
      if(v.tail(L-j-1).norm()>tol)
      {
        if(v(0)>0)
        {
          v(0)=v(0)+max_norm;
        }
        else
        {
          v(0)=v(0)-max_norm;
        }
        v.noalias()=v/v.norm();
        C.block(j,j,L-j,n-j).noalias()-=2*v*(v.transpose()*C.block(j,j,L-j,n-j));
      }
      rank++;
    }
    else
    {
      break;
    }
  }
  for(int j=i_switch;j>=0;j--)
  {
    int j1=switch_mat(j,0), j2=switch_mat(j,1);
    VectorXd tmp=C.col(j1);
    C.col(j1)=C.col(j2);
    C.col(j2)=tmp;
  }
  
  return C_ini.topRows(rank);
}






/////////////////////////////////////////
List cv_fix_level(Map<MatrixXd> G, VectorXd Y, VectorXd alpha_set, VectorXd lambda_set, VectorXi level_vec, List all_folds, double shift, double lower_bound)
{
  int L=G.cols(), n=G.rows(), K_cv=all_folds.size();
  MatrixXd cv_error=MatrixXd::Zero(alpha_set.size(), lambda_set.size());
  MatrixXd B_ini(L, n);
  Map<MatrixXd> B=Map<MatrixXd>(B_ini.data(), L,n);
  
  for(int i_alpha=0; i_alpha<alpha_set.size(); i_alpha++)
  {
    double alpha=alpha_set(i_alpha);
    for(int i=0; i<L; i++)
    {
      double tmp;
      tmp=pow(pow(2, 2*alpha*(level_vec(i)-exp(-(level_vec(i)-shift)/alpha)))+lower_bound, -0.5);
      B.row(i).noalias()=G.col(i)*tmp;
    }
    MatrixXd R_ini=QR_house(B);
    Map<MatrixXd> R=Map<MatrixXd>(R_ini.data(), R_ini.rows(), R_ini.cols());
    for(int i_cv=0; i_cv<K_cv; i_cv++)
    {
      VectorXi omit;
      omit=as<VectorXi>(all_folds(i_cv)).array()-1;
      int n_train=n-omit.size();
      MatrixXd R2_ini(R.rows(), omit.size());
      Map<MatrixXd> R2=Map<MatrixXd>(R2_ini.data(), R2_ini.rows(), R2_ini.cols());
      MatrixXd R1_ini(R.rows(), n_train);
      Map<MatrixXd> R1=Map<MatrixXd>(R1_ini.data(), R1_ini.rows(), R1_ini.cols());
      VectorXi ind=VectorXi::Zero(n);
      VectorXd Y_train(n_train), Y_valid(omit.size());
      
      for(int i=0; i<omit.size(); i++)
      {
        ind(omit(i))=1;
        R2.col(i)=R.col(omit(i));
        Y_valid(i)=Y(omit(i));
      }
      int i_train=0;
      for(int i=0; i<n; i++)
      {
        if(ind(i)==0)
        {
          R1.col(i_train)=R.col(i);
          Y_train(i_train)=Y(i);
          i_train=i_train+1;
        }
      }
      VectorXd R1_mean=R1.rowwise().mean();
      R1=R1.colwise()-R1_mean;
      R2=R2.colwise()-R1_mean;
      
      Y_valid=Y_valid.array()-Y_train.mean();
      Y_train=Y_train.array()-Y_train.mean();
      List fit_QR=QR_full(R1);
      
      MatrixXd S1_ini=fit_QR(0), P1_ini=fit_QR(1);
      
      Map<MatrixXd> S1=Map<MatrixXd>(S1_ini.data(), S1_ini.rows(), S1_ini.cols());
      Map<MatrixXd> P1=Map<MatrixXd>(P1_ini.data(), P1_ini.rows(), P1_ini.cols());
      MatrixXd D_ini=S1*S1.transpose();
      Map<MatrixXd> D=Map<MatrixXd>(D_ini.data(), D_ini.rows(), D_ini.cols());
      
      for(int i_lambda=0; i_lambda<lambda_set.size(); i_lambda++)
      {
        double lambda=lambda_set(i_lambda);
        D.diagonal()=D.diagonal().array()+lambda;
        VectorXd Y_pred=R2.transpose()*(P1*D.ldlt().solve(S1*Y_train));
        cv_error(i_alpha, i_lambda)=cv_error(i_alpha, i_lambda)+(Y_pred-Y_valid).squaredNorm();
      }
    }
  }
  return  List::create(_["cv_error"]=cv_error);
}




/////////////////////////////////////////
// [[Rcpp::export]]
List c_cv_spiky_scalar_on_function(List t_x_list, List X, Eigen::VectorXd Y, List x_params, List all_folds)
{
  int ncurves=X.size(), max_level=x_params(3);
  int nsample=Y.size();
  MatrixXd wave_mat=as<MatrixXd>(x_params(4));
  VectorXd lower_bound_set=as<VectorXd>(x_params(7));
  VectorXd father=wave_mat.col(0), mother=wave_mat.col(1);
  VectorXi nbasis_levels=as<VectorXi>(x_params(6));
  VectorXd alpha_set=as<VectorXd>(x_params(1)),   lambda_set=as<VectorXd>(x_params(2));
  int upper=as<int>(x_params(5));
  VectorXd ref_t=VectorXd::LinSpaced(father.size(), 0, upper);
  List G_list(ncurves);
  VectorXi ind_level=VectorXi::Zero(nbasis_levels(max_level));
  double min_cv_error=1e20;
  int opt_level;
  double opt_shift, opt_lower_bound;
  VectorXi opt_level_vec;
  MatrixXd min_cv_mat;
  for(int i_level=1; i_level<=max_level; i_level++)
  {
    ind_level.segment(nbasis_levels(i_level-1), nbasis_levels(i_level)-nbasis_levels(i_level-1))=VectorXi::Constant(nbasis_levels(i_level)-nbasis_levels(i_level-1), i_level);
  }
  for(int i=0; i<ncurves;i++)
  {
    int ind=0;
    MatrixXd G_tmp=MatrixXd::Zero(nsample, nbasis_levels(max_level));
    ArrayXd t_x=as<ArrayXd>(t_x_list(i));
    MatrixXd X_tmp=as<MatrixXd>(X(i));
    Map<MatrixXd> Xmap=Map<MatrixXd>(X_tmp.data(), X_tmp.rows(), X_tmp.cols());
    int t_size=Xmap.cols();
    for(int k=1-upper; k<1;k++)
    {
      G_tmp.col(ind)=Xmap*eval_wave_basis(t_x-k, father, ref_t, upper)/t_size;
      ind=ind+1;
    }
    for(int i_level=0; i_level<=max_level; i_level++)
    {
      int N=pow_int(2, i_level);
      int m0, m1;
      for(int k=1-upper; k<N;k++)
      {
        m0=floor(k*(t_size-1)/N)-1;
        if(m0<0)
        {
          m0=0;
        }
        m1=ceil((upper+k)*(t_size-1)/N)+1;
        
        if(m1>t_size-1)
        {
          m1=t_size-1;
        }
        G_tmp.col(ind)=sqrt(N)*(Xmap.block(0,m0, nsample, m1-m0+1)*eval_wave_basis(N*t_x.segment(m0, m1-m0+1)-k, mother, ref_t, upper))/t_size;
        ind=ind+1;
      }
    }
    
    G_list(i)=G_tmp;
  }
  
  for(int i_lower=0; i_lower<lower_bound_set.size(); i_lower++)
  {
    MatrixXd cv_error;
    double lower_bound=lower_bound_set(i_lower);
    for(int i_level=1; i_level<=max_level; i_level++)
    {
      int col_size=nbasis_levels(i_level);
      MatrixXd G_tmp(nsample, ncurves*col_size);
      Map<MatrixXd> G=Map<MatrixXd>(G_tmp.data(), G_tmp.rows(), G_tmp.cols());
      VectorXi level_vec(ncurves*col_size);
      for(int j=0;j<ncurves;j++)
      {
        G.block(0, j*col_size, nsample, col_size)=as<MatrixXd>(G_list(j)).leftCols(col_size);
        level_vec.segment(j*col_size, col_size)=ind_level.head(col_size);
      }
      VectorXd shift_set=VectorXd::LinSpaced(2*i_level-1, 0.5, i_level-0.5);
      for(int i_shift=0; i_shift<shift_set.size(); i_shift++)
      {
        double shift=shift_set(i_shift);
        List fit_fix_level=cv_fix_level(G, Y, alpha_set, lambda_set, level_vec, all_folds, shift, lower_bound);
        cv_error=as<MatrixXd>(fit_fix_level(0));
        double tmp=cv_error.minCoeff();
        if((i_level>opt_level)&(opt_lower_bound==lower_bound))
        {
          if(min_cv_error*0.95>tmp)
          {
            min_cv_error=tmp;
            min_cv_mat=cv_error;
            opt_level=i_level;
            opt_shift=shift;
            opt_lower_bound=lower_bound;
            opt_level_vec=level_vec;
          }
        }
        else
        {
          if(opt_lower_bound>lower_bound)
          {
            if(min_cv_error*0.95>tmp)
            {
              min_cv_error=tmp;
              min_cv_mat=cv_error;
              opt_level=i_level;
              opt_shift=shift;
              opt_lower_bound=lower_bound;
              opt_level_vec=level_vec;
            }
          }
          else
          {
            if(min_cv_error>tmp)
            {
              min_cv_error=tmp;
              min_cv_mat=cv_error;
              opt_level=i_level;
              opt_shift=shift;
              opt_lower_bound=lower_bound;
              opt_level_vec=level_vec;
            }
          }
        }
      }
    }
  }
  
  int col_size=nbasis_levels(opt_level);
  MatrixXd G_tmp(nsample, ncurves*col_size);
  Map<MatrixXd> G=Map<MatrixXd>(G_tmp.data(), G_tmp.rows(), G_tmp.cols());
  
  VectorXi level_vec(ncurves*col_size);
  
  for(int j=0;j<ncurves;j++)
  {
    G.block(0, j*col_size, nsample, col_size)=as<MatrixXd>(G_list(j)).leftCols(col_size);
    level_vec.segment(j*col_size, col_size)=ind_level.head(col_size);
  }
  G=G.rowwise()-G.colwise().mean();
  MatrixXd::Index i_alpha, i_lambda;
  min_cv_error=min_cv_mat.minCoeff(&i_alpha, &i_lambda);
  double alpha=alpha_set(i_alpha);
  double lambda=lambda_set(i_lambda);
  
  
  int L=G.cols(), n=G.rows();
  MatrixXd B_ini(L, n);
  Map<MatrixXd> B=Map<MatrixXd>(B_ini.data(), L,n);
  VectorXd K_inv_half(L);
  for(int i=0; i<L; i++){
    K_inv_half(i)= pow(pow(2, 2*alpha*(level_vec(i)-exp(-(level_vec(i)-opt_shift)/alpha)))+opt_lower_bound, -0.5);
    B.row(i).noalias()=G.col(i)*K_inv_half(i);
  }
  
  MatrixXd R_ini=QR_house(B);
  Map<MatrixXd> R=Map<MatrixXd>(R_ini.data(), R_ini.rows(), R_ini.cols());
  MatrixXd E_ini=R*R.transpose();
  Map<MatrixXd> E=Map<MatrixXd>(E_ini.data(), E_ini.rows(), E_ini.cols());
  double mu=Y.mean();
  Y=Y.array()-mu;
  VectorXd c=B*(R.transpose()*((E*E+lambda*E).ldlt().solve(R*Y)));
  c=c.array() * K_inv_half.array();
  
  List coef(ncurves);
  
  for(int i=0; i<ncurves;i++)
  {
    int ind=0;
    
    ArrayXd t_x=as<ArrayXd>(t_x_list(i));
    int t_size=t_x.size();
    VectorXd  coef_vec=VectorXd::Zero(t_size);
    for(int k=1-upper; k<1;k++)
    {
      coef_vec.noalias()+=eval_wave_basis(t_x-k, father, ref_t, upper)*c(ind+i*col_size);
      ind=ind+1;
    }
    for(int i_level=0; i_level<=opt_level; i_level++)
    {
      int N=pow_int(2, i_level);
      for(int k=1-upper; k<N;k++)
      {
        coef_vec.noalias()+=sqrt(N)*eval_wave_basis(N*t_x-k, mother, ref_t, upper)*c(ind+i*col_size);
        ind=ind+1;
      }
    }
    MatrixXd X_tmp=as<MatrixXd>(X(i));
    mu=mu-coef_vec.dot(X_tmp.colwise().mean())/t_size;
    coef(i)=coef_vec;
  }
  return  List::create(_["error_mat"]=min_cv_mat, _["min_cv_error"]=min_cv_error, _["opt_level"]=opt_level , _["opt_shift"]=opt_shift,
                       _["opt_lower_bound"]=opt_lower_bound, _["opt_level_vec"]=level_vec, _["opt_alpha"]=alpha,  _["opt_lambda"]=lambda, _["coef"]=coef, _["mu"]=mu);
}

