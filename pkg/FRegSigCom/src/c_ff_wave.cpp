#include <RcppEigen.h>
#include <Rcpp.h>
#include "common_function.h"
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;


/////////////////////////////////////////
// [[Rcpp::export]]
double max_eig(Eigen::MatrixXd A)
{
  int p=A.rows();
  NumericVector xx = runif(p); 
  VectorXd v=as<Eigen::Map<Eigen::VectorXd> >(xx);
  // Rcout << "v="<< v.transpose() << std::endl;
  int count=1;
  v=v/v.norm();
  double value=v.transpose()*(A*v);
  double old_value=value/2.0;
  while((count<100)&((value-old_value)>1e-5*old_value))
  {
    v=A*v;
    v=v/v.norm();
    old_value=value;
    value=v.transpose()*(A*v);
    count++;
  }
  return(value);
}
 
////////////////////////////
double sign(double x)
{
  double out=-1;
  if(x>0)
  {
    out=1;
  }
  return(out);
}

/////////////////////////////////////////
// [[Rcpp::export]]
List max_H(Eigen::VectorXd gamma0, Eigen::MatrixXd Psi0, Eigen::VectorXd t0, int nvar, double tau, double lambda)
{
  double tol=1e-12;
  Map<VectorXd> gamma=Map<VectorXd>(gamma0.data(), gamma0.size());
  Map<VectorXd> t=Map<VectorXd>(t0.data(), t0.size());
  Map<MatrixXd> Psi=Map<MatrixXd>(Psi0.data(), Psi0.rows(), Psi0.cols());
  int dim_old=Psi.cols(), d=Psi.rows();
  VectorXd temp=gamma+Psi*t;   
  
  MatrixXd Psi_1_ini=Psi.topRows(d-nvar);
  Map<MatrixXd> Psi_1=Map<MatrixXd>(Psi_1_ini.data(), Psi_1_ini.rows(), Psi_1_ini.cols());
  MatrixXd Psi_2_ini=Psi.bottomRows(nvar);
  Map<MatrixXd> Psi_2=Map<MatrixXd>(Psi_2_ini.data(), Psi_2_ini.rows(), Psi_2_ini.cols());
  VectorXd b=temp.head(d-nvar);
  VectorXd a=temp.tail(nvar);
  List list_1=basic_max(a,lambda);
  VectorXd x=list_1(0), x_new;
  double nu=as<double>(list_1(2))/2, nu_new;
  int m=list_1(1), m_new;
  VectorXd h=tau*Psi_1.transpose()*b+(2*nu/(1-lambda))*Psi_2.transpose()*x, h_new;
  double H=h.squaredNorm(), H_new;
  double alpha=lambda/(1+(m-1)*lambda), alpha_new;
  int  count1=0;
  MatrixXd A_ini;
  double delta=0;
  VectorXd t_inc, t_new;
  while((H>tol) & (count1<100))
  {
    count1++;
    VectorXi I=as<VectorXi>(list_1(3));
    MatrixXd NNPsi_ini(I.size(), Psi_2.cols());
    Map<MatrixXd> NNPsi=Map<MatrixXd>(NNPsi_ini.data(), NNPsi_ini.rows(), NNPsi_ini.cols());
    VectorXd signx(I.size());
    for(int j=0; j<I.size(); j++)
    {
      signx(j)=sign(x(I(j)));
      NNPsi.row(j)=Psi_2.row(I(j));
    }
    if((I.size()<nvar/2.0) | (I.size()==nvar))
    { 
      MatrixXd K_ini= NNPsi - alpha*(signx*(signx.transpose()*NNPsi));
      Map<MatrixXd> K=Map<MatrixXd>(K_ini.data(), K_ini.rows(), K_ini.cols());
      A_ini=tau*Psi_1.transpose()*Psi_1+K.transpose()*NNPsi/(1-lambda);
    }
    else
    { 
      VectorXi I0=as<VectorXi>(list_1(4));
      MatrixXd NNPsi0_ini(I0.size(), Psi_2.cols());
      Map<MatrixXd> NNPsi0=Map<MatrixXd>(NNPsi0_ini.data(), NNPsi0_ini.rows(), NNPsi0_ini.cols());
      for(int j=0; j<I0.size(); j++)
      {
        NNPsi0.row(j)=Psi_2.row(I0(j));
      }
      temp=NNPsi.transpose()*signx;
      A_ini=(tau-1/(1-lambda))*Psi_1.transpose()*Psi_1+(MatrixXd::Identity(dim_old, dim_old)-NNPsi0.transpose()*NNPsi0-alpha*temp*temp.transpose())/(1-lambda);
    }
    MatrixXd tmp_A=A_ini.transpose();
    A_ini=(A_ini+tmp_A)/2;
    Map<MatrixXd> A=Map<MatrixXd>(A_ini.data(), A_ini.rows(), A_ini.cols());
    VectorXd deriv_1=2*A*h;
    
    if(count1==1)
    {
      delta=0.01*max_eig(A);
    }
    if(H>1e-6)
    {
      t_inc=-(A+delta*MatrixXd::Identity(A.rows(), A.rows())).llt().solve(h);
    }
    else
    {
      t_inc=-A.llt().solve(h);
    }
    t_new=t+t_inc;
    
    temp=gamma+Psi*t_new;
    b=temp.head(d-nvar);
    a=temp.tail(nvar);
    list_1=basic_max(a,lambda);
    x_new=list_1(0);
    nu_new=as<double>(list_1(2))/2;
    m_new=list_1(1);
    
    alpha_new=lambda/(1+(m_new-1)*lambda);
    h_new = tau*Psi_1.transpose()*b+(2*nu_new/(1-lambda))*Psi_2.transpose()*x_new;  
    H_new=h_new.squaredNorm();
    int count2=0;
    while((H_new>tol)&(H_new>H)&(count2<5))
    {
      count2=count2+1;
      t_inc=0.3*t_inc;
      t_new=t+t_inc;
      
      if(count2 > 4)
      {
        NumericVector xx = runif(t.size()); 
        temp=as<Eigen::Map<Eigen::VectorXd> >(xx);
         
        t_new=t+0.1*temp/temp.norm();
        temp=gamma+Psi*t_new;
        b=temp.head(d-nvar);
        a=temp.tail(nvar);
        list_1=basic_max(a,lambda);
        x_new=list_1(0);
        nu_new=as<double>(list_1(2))/2;
        m_new=list_1(1);
        alpha_new=lambda/(1+(m_new-1)*lambda);
        h_new = tau*Psi_1.transpose()*b+(2*nu_new/(1-lambda))*Psi_2.transpose()*x_new;  
        H_new=h_new.squaredNorm();
        break;
      }
      
      
      temp=gamma+Psi*t_new;
      b=temp.head(d-nvar);
      a=temp.tail(nvar);
      list_1=basic_max(a,lambda);
      x_new=list_1(0);
      nu_new=as<double>(list_1(2))/2;
      m_new=list_1(1);
      alpha_new=lambda/(1+(m_new-1)*lambda);
      h_new = tau*Psi_1.transpose()*b+(2*nu_new/(1-lambda))*Psi_2.transpose()*x_new;  
      H_new=h_new.squaredNorm();
    }
    
    t=t_new;
    h=h_new;
    nu=nu_new;
    H = H_new;
    x = x_new;
    m=m_new;
    alpha = alpha_new;
  } 
  temp=gamma+Psi*t;
  b=temp.head(d-nvar);
  a=temp.tail(nvar);
  double t1=2*nu/(1-lambda);
  double  t2=tau*b.squaredNorm();
  VectorXd v(b.size()+x.size());
  double tmp=1.0/pow(pow(t1, 2.0)+t2, 0.5);
  v.head(b.size())=pow(tau, 0.5)*tmp*b;
  v.tail(x.size())=tmp*t1*pow(tau, -0.5)*x;
  return List::create(_["x"]=v, _["t"]=t);
}

 


