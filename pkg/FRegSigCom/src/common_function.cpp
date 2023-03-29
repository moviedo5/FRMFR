#include <Rcpp.h>
#include <RcppEigen.h>
#include "common_function.h"

// [[Rcpp::depends(RcppEigen)]]


using namespace Rcpp;
using namespace std;
using namespace Eigen;

/////////////////////////////////////////////////////
 
void power_eig(Eigen::Map<Eigen::MatrixXd> Pi_eig, Eigen::Map<Eigen::MatrixXd> M, Eigen::Map<Eigen::MatrixXd> MM_inv, Eigen::MatrixXd J, double & max_value, Eigen::Map<Eigen::VectorXd>& beta_current)
{
  double tol_1=1e-10, tol_2=1e-4;
  int r=Pi_eig.rows(), p=M.rows();
  VectorXd tmp_ini=VectorXd::Constant(p, 1), b_ini=tmp_ini/tmp_ini.norm(), b_old_ini(p), tmp1_ini=VectorXd::Zero(p);
  Map<VectorXd> tmp=Map<VectorXd>(tmp_ini.data(), tmp_ini.size());
  Map<VectorXd> b=Map<VectorXd>(b_ini.data(), b_ini.size());
  Map<VectorXd> b_old=Map<VectorXd>(b_old_ini.data(), b_old_ini.size());
  Map<VectorXd> tmp1=Map<VectorXd>(tmp1_ini.data(), tmp1_ini.size());
  MatrixXd J_orth=J-M*(MM_inv*(M.transpose()*J));
  MatrixXd JJ_inv=(J_orth.transpose()*J_orth).inverse();
  int count=0, ind=0;
  double old_max_value=0;
  double e1, e2;
  MatrixXd w0=J.topRows(r);
  while((count<40)&&(ind==0))
  {
    b_old=b;
    tmp=b;
    tmp.noalias()-= M*(MM_inv*(M.transpose()*b));
    tmp.noalias()-= J_orth*(JJ_inv*(J_orth.transpose()*b));
    tmp1.head(r).noalias()=Pi_eig*tmp.head(r);
    b=tmp1;
    b.noalias()-=M*(MM_inv*(M.transpose()*tmp1));
    b.noalias()-=J_orth*(JJ_inv*(J_orth.transpose()*tmp1));
    max_value=b_old.dot(b);
    b=b/b.norm();
    e1=(M.transpose()*b).cwiseAbs().maxCoeff()+(w0.transpose()*b.head(r)).cwiseAbs().maxCoeff();
    e2=(max_value-old_max_value)/max_value;
    //e1=(b-b_old).norm();
    //e2=(b+b_old).norm();
    if((e1<tol_1)&&(e2<tol_2))
      //if((e1<1e-8)||(e2<1e-8))
    {
      ind=1;
    }
    old_max_value=max_value;
    count++;
  }
  beta_current=b;
}


void  power_eig_0(Eigen::Map<Eigen::MatrixXd> Pi_eig, Eigen::Map<Eigen::MatrixXd> M, Eigen::Map<Eigen::MatrixXd> MM_inv, double & max_value, Eigen::Map<Eigen::VectorXd>& beta_current)
//void  power_eig_0(Map<MatrixXd> Pi_eig, Map<MatrixXd> M, Map<MatrixXd> MM_inv, double & max_value, Map<VectorXd> & beta_current)//function used to calculate the first eigenvalue of Pi which is orthogonal to M
{
  double tol_1=1e-10, tol_2=1e-4;
  int r=Pi_eig.rows(), p=M.rows();
  VectorXd tmp_ini=VectorXd::Constant(p, 1), b_ini=tmp_ini/tmp_ini.norm(), b_old_ini(p), tmp1_ini=VectorXd::Zero(p);
  Map<VectorXd> tmp=Map<VectorXd>(tmp_ini.data(), tmp_ini.size());
  Map<VectorXd> b=Map<VectorXd>(b_ini.data(), b_ini.size());
  Map<VectorXd> b_old=Map<VectorXd>(b_old_ini.data(), b_old_ini.size());
  Map<VectorXd> tmp1=Map<VectorXd>(tmp1_ini.data(), tmp1_ini.size());
  int count=0, ind=0;
  double old_max_value=0;
  double e1, e2;
  while((count<40)&&(ind==0))
  {
    b_old=b;
    tmp=b;
    tmp.noalias()-=M*(MM_inv*(M.transpose()*b));
    tmp1.head(r).noalias()=Pi_eig*tmp.head(r);
    b=tmp1;
    b.noalias()-=M*(MM_inv*(M.transpose()*tmp1));
    max_value=b_old.dot(b);
    b=b/b.norm();
    e1=(M.transpose()*b).cwiseAbs().maxCoeff();
    e2=(max_value-old_max_value)/max_value;
    //e1=(b-b_old).norm();
    //e2=(b+b_old).norm();
    if((e1<tol_1)&&(e2<tol_2))
      //if((e1<1e-8)||(e2<1e-8))
    {
      ind=1;
    }
    old_max_value=max_value;
    count++;
  }
  beta_current=b;
}


//////////////////////////////////////////////

/////////////////////////////////////////////////////

List cal_comp_without_max(Map<MatrixXd> Pi_eig, Map<MatrixXd> M, int upper_comp, double thresh)
{
  int r=Pi_eig.cols(), p=M.rows();
  MatrixXd MM_inv_ini=(M.transpose()*M).inverse();
  Map<MatrixXd> MM_inv=Map<MatrixXd>(MM_inv_ini.data(), MM_inv_ini.rows(), MM_inv_ini.cols());
  VectorXd  beta_current_ini=VectorXd::Zero(p);
  Map<VectorXd> beta_current=Map<VectorXd>(beta_current_ini.data(), p);
  VectorXd max_value_vec=VectorXd::Zero(upper_comp);
  MatrixXd beta=MatrixXd::Zero(p, upper_comp);
  int n_comp;
  double denominator=Pi_eig.diagonal().sum();
  for(n_comp=0; n_comp<upper_comp; n_comp++)
  {
    double max_value;
    if(n_comp==0)
    {
      power_eig_0(Pi_eig, M, MM_inv, max_value, beta_current);
    }
    else
    {
      MatrixXd J=MatrixXd::Zero(p, n_comp);
      J.topRows(r).noalias()=beta.topLeftCorner(r,n_comp-1);
      power_eig(Pi_eig, M,  MM_inv, J, max_value, beta_current);
    }
    max_value_vec(n_comp)=max_value;
    beta.col(n_comp).noalias()=beta_current;
    if((max_value_vec(n_comp)/denominator<thresh))
      //if((max_value_vec(n_comp)/max_value_vec.head(n_comp+1).sum()<thresh))
    {break;}
    if(n_comp==upper_comp-1)
    {break;}
  }
  return List::create( _["beta"]=beta.leftCols(n_comp+1),  _["max_value"]=max_value_vec);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////////

List cal_comp_with_max(Map<MatrixXd> Pi_eig, Map<MatrixXd> M, int max_comp)
{
  int r=Pi_eig.cols(), p=M.rows();
  MatrixXd MM_inv_ini=(M.transpose()*M).inverse();
  Map<MatrixXd> MM_inv=Map<MatrixXd>(MM_inv_ini.data(), MM_inv_ini.rows(), MM_inv_ini.cols());
  VectorXd  beta_current_ini=VectorXd::Zero(p);
  Map<VectorXd> beta_current=Map<VectorXd>(beta_current_ini.data(), p);
  MatrixXd beta=MatrixXd::Zero(p, max_comp);
  int n_comp;
  for(n_comp=0; n_comp<max_comp; n_comp++)
  {
    double max_value;
    if(n_comp==0)
    {
      power_eig_0(Pi_eig, M, MM_inv, max_value, beta_current);
    }
    else
    {
      MatrixXd J=MatrixXd::Zero(p, n_comp);
      J.topRows(r).noalias()=beta.topLeftCorner(r,n_comp-1);
      power_eig(Pi_eig, M,  MM_inv, J, max_value, beta_current);
    }
    beta.col(n_comp).noalias()=beta_current;
  }
  return List::create( _["beta"]=beta);
}

///////////////////////////////////////////////////////////////////////




/////////////////////////////////////////
VectorXi extract(VectorXi x, VectorXi ind)
{
  VectorXi out(ind.sum());
  int j=0;
  for(int i=0; i<ind.size(); i++)
  {
    if(ind(i)!=0)
    {
      out(j)=x(i);
      j=j+1;
    }
  }
  return(out);
}
////////////
/////////////////
// [[Rcpp::export]]
List basic_max(Eigen::VectorXd a, double lambda)
{
  int p=a.size();
  VectorXd a_abs=a.array().abs();
  VectorXd aa=a_abs;
  sort(aa.data(), aa.data()+aa.size(), greater<double>());
  int nonzero_num=0;
  double S=0;
  for(int i=0;i<p; ++i)
  {
    S+=aa(i);
    if(aa(i)<=lambda*S/(1+i*lambda))
    {
      S-=aa(i);
      nonzero_num=i;
      break;
    }
    nonzero_num=i+1;
  }
  S=lambda*S/(1+(nonzero_num-1)*lambda);
  VectorXd x=VectorXd::Zero(p);
  int nonzero=0, ix_zero=0;
  VectorXi id_nonzero(p), id_zero(p);
  for(int i=0; i<p; i++)
  {
    if(a(i)>S)
    {
      x(i)=a(i)-S;
      id_nonzero(nonzero)=i;
      nonzero++;
    }
    else
    {
      if(a(i)<-S)
      {
        x(i)=a(i)+S;
        id_nonzero(nonzero)=i;
        nonzero++;
      }
      else
      {
        id_zero(ix_zero)=i;
        ix_zero++;
      }
    }
  }
  
  double x_norm=pow((1-lambda)*x.squaredNorm()+lambda*pow(x.lpNorm<1>(),2.0),0.5);
  x=x/x_norm;
  return List::create( _["x"]=x, _["m"]=nonzero, _["norm"]=x_norm, _["id_nonzero"]=id_nonzero.head(nonzero),
                       _["id_zero"]=id_zero.head(ix_zero));
}



/////////////////////////////////////////
VectorXd eval_wave_basis(ArrayXd t, VectorXd ref_vec, VectorXd ref_t, int upper)
{
  int  N=ref_vec.size(), m=t.size();
  VectorXd out(m);
  double scale=N/upper;
  for(int i=0; i<m; i++)
  {
    double tmp=t(i);
    if((tmp>0)&&(tmp<upper))
    {
      int k=(int)(tmp*scale);
      if(k>N-2)
      {
        k=N-2;
      }
      out(i)=scale*(ref_vec(k)*(ref_t(k+1)-tmp)+ref_vec(k+1)*(tmp-ref_t(k)));
    }
    else
    {
      out(i)=0;
    }
  }
  return out;
}

////////////////////////////
int pow_int(int base, int power)
{
  int out=1;
  for(int k=1; k<=power; k++)
  {
    out*=base;
  }
  return out;
}

////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////
// [[Rcpp::export]]

Eigen::MatrixXd find_orth_basis(Eigen::MatrixXd X)
{
  MatrixXd orthConst_mtx_ini=X;
  Map<MatrixXd> orthConst_mtx=Map<MatrixXd>(orthConst_mtx_ini.data(), orthConst_mtx_ini.rows(), orthConst_mtx_ini.cols());
  int m=X.rows();
  orthConst_mtx.row(0)=orthConst_mtx.row(0)/((orthConst_mtx.row(0)).norm());
  for(int i=1;i<m;++i)
  {
    VectorXd tmp=orthConst_mtx.bottomRows(m-i)*(orthConst_mtx.row(i-1).transpose());
    orthConst_mtx.bottomRows(m-i)=orthConst_mtx.bottomRows(m-i)-tmp*orthConst_mtx.row(i-1);
    orthConst_mtx.row(i)=orthConst_mtx.row(i)/((orthConst_mtx.row(i)).norm());
  }
  return(orthConst_mtx);
}

