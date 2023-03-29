#include <RcppEigen.h>
#include <Rcpp.h>
#include "common_function.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;

///////////////////////////////////////////////////////////////


List maxLinearH(VectorXd gamma, int nobs, int nbasis, int p, double tau, double lambda)
{
  VectorXd b=gamma.head(nobs);
  VectorXd a(p);
  VectorXd x=VectorXd::Zero(nobs+p*nbasis);

  for(int i=0;i<p;++i)
  {
    a(i)=gamma.segment(i*nbasis+nobs, nbasis).norm();
  }
  List fit_1=basic_max(a, lambda);
  VectorXd u=fit_1(0);
  VectorXi id_nonzero=fit_1(3);
  double tmp=u.dot(a);
  double normc=pow(tau*b.squaredNorm()+pow(tmp, 2.0),0.5);

  VectorXi id=fit_1(3);
  int k=id.size();
  for(int i=0;i< k;++i)
  {
    x.segment(id(i)*nbasis+nobs, nbasis)=gamma.segment(id(i)*nbasis+nobs, nbasis)*u(id(i))*tmp/a(id(i))/pow(tau, 0.5)/normc;
  }

  x.segment(0, nobs)=pow(tau,0.5)*b/normc;

  return List::create( _["x"]=x, _["a"]=a, _["u"]=u, _["normc"]=normc, _["id_nonzero"]=id);
}


///////////////////////////////////////////////////////////////



List constOptOrth_H(VectorXd gamma_ini, MatrixXd B1_ini, MatrixXd orthConst_mtx_ini, VectorXd t0, List Lambda_all, int nobs, int nbasis, int nvarX, double tau, double lambda, double tol)
{

  Map<MatrixXd> orthConst_mtx=Map<MatrixXd>(orthConst_mtx_ini.data(), orthConst_mtx_ini.rows(), orthConst_mtx_ini.cols());
  Map<MatrixXd> B1=Map<MatrixXd>(B1_ini.data(), B1_ini.rows(), B1_ini.cols());
  int  kcomp=orthConst_mtx.rows();
  Map<VectorXd> gamma=Map<VectorXd>(gamma_ini.data(), gamma_ini.size());

  Map<VectorXd> t_new=Map<VectorXd>(t0.data(), t0.size());
  VectorXd t_ini(t0.size());
  Map<VectorXd> t=Map<VectorXd>(t_ini.data(), t_ini.size());


  VectorXd tmpXd=-orthConst_mtx.transpose()*(orthConst_mtx*gamma);
  VectorXd gamma_orth=gamma+tmpXd;
  VectorXd gammat=gamma_orth, gammat_new;
  double normc, normc_new, tau_sqrt=pow(tau, 0.5);
  VectorXd x, a, u, x_new, a_new, u_new;
  List fit_2=maxLinearH(gammat, nobs, nbasis, nvarX, tau, lambda);
  x=as<VectorXd>(fit_2(0));
  a=as<VectorXd>(fit_2(1));
  u=as<VectorXd>(fit_2(2));
  normc=as<double>(fit_2(3));
  VectorXi I=as<VectorXi>(fit_2(4)), I_new;

  tmpXd=orthConst_mtx*x;
  VectorXd h=normc*tau_sqrt*tmpXd, h_new;
  double H=h.squaredNorm(), H_new;
  if(t0.lpNorm<1>()>0)
  {
  //  Rcout << "start 2.6 in constOptOrth_H"   << std::endl;
    gammat_new=gamma_orth+orthConst_mtx.transpose()*t_new;
  //  Rcout << "start 2.7 in constOptOrth_H"   << std::endl;
    fit_2=maxLinearH(gammat_new , nobs, nbasis, nvarX, tau, lambda);
    x_new=as<VectorXd>(fit_2(0));
    a_new=as<VectorXd>(fit_2(1));
    u_new=as<VectorXd>(fit_2(2));
    normc_new=as<double>(fit_2(3));
    I_new=as<VectorXi>(fit_2(4));
    h_new=normc_new*tau_sqrt*(orthConst_mtx*x_new);
    H_new=h_new.squaredNorm();
  //  Rcout << "start 2.8 in constOptOrth_H"   << std::endl;
    if(H_new<H)
    {
      gammat=gammat_new;
      t=t_new;
      x=x_new;
      a=a_new;
      u=u_new;
      normc=normc_new;
      I=I_new;
      h=h_new;
      H=H_new;
    }
  }

 // Rcout << "start 3 in constOptOrth_H"   << std::endl;
  int count1=0;
  int stop_id;
  while((H>tol) & (count1<50))
  { //Rcout << "count1 " <<" "<< count1  << std::endl;
   // Rcout << "H= " <<" "<< H  << std::endl;
    count1=count1+1;


    double alpha=lambda/(1+(I.size()-1)*lambda);
    MatrixXd Xi_ini(kcomp, I.size());
    Map<MatrixXd> Xi=Map<MatrixXd>(Xi_ini.data(), Xi_ini.rows(), Xi_ini.cols());
    VectorXd tmpvec;
    for(int i=0;i<I.size();++i)
    {
      tmpvec=x.segment(I(i)*nbasis+nobs, nbasis);
      Xi.col(i)=orthConst_mtx.block(0, I(i)*nbasis+nobs, kcomp, nbasis)*tmpvec/tmpvec.norm();
    }

    double nu2=(1-lambda)*(a.dot(u));

    MatrixXd tmp_Xi(Xi.rows(), Xi.cols());
    MatrixXd B4=MatrixXd::Zero(kcomp, kcomp);
    for(int i=0;i<I.size();++i)
    {
      double tmp_val=nu2*u(I(i))/a(I(i));
      tmp_Xi.col(i)=(1-tmp_val)*Xi.col(i);
      B4+=as<MatrixXd>(Lambda_all(I(i)))*tmp_val;
    }
    tmpvec=Xi*VectorXd::Constant(I.size(), 1.0);
    MatrixXd A_ini;
    A_ini=tau*B1+1/(1-lambda)*(tmp_Xi*Xi.transpose()-alpha*tmpvec*tmpvec.transpose()+B4);
    MatrixXd tmp_A=A_ini.transpose();
    A_ini=(A_ini+tmp_A)/2;
    Map<MatrixXd> A=Map<MatrixXd>(A_ini.data(), A_ini.rows(), A_ini.cols());

     VectorXd t_inc=-A.llt().solve(h);
    if(H>1e6)
    {
      if(H>1e8)
      {
        if(H>1e10)
        {
          t_inc=t_inc/6;
        }else
        {
          t_inc=t_inc/4;
        }
      }else
      {
        t_inc=t_inc/2;
      }
    }else
    {
      t_inc=t_inc;
    }


    MatrixXd delta_gamma= orthConst_mtx.transpose()*t_inc;

    gammat_new=gammat+delta_gamma;

    fit_2=maxLinearH(gammat_new, nobs, nbasis, nvarX, tau, lambda);
    x_new=as<VectorXd>(fit_2(0));
    a_new=as<VectorXd>(fit_2(1));
    u_new=as<VectorXd>(fit_2(2));
    normc_new=as<double>(fit_2(3));
    I_new=as<VectorXi>(fit_2(4));

    h_new=normc_new*tau_sqrt*(orthConst_mtx*x_new);
    H_new=h_new.squaredNorm();


    VectorXd gammat_tmp;
    //Rcout << "H_new= "<<" "<< H_new  << std::endl;
    //Rcout << "H= "<<" "<< H  << std::endl;
    if(H_new<H)
    {//Rcout << "start 4 in constOptOrth_H"   << std::endl;
      t=t+t_inc;
      gammat=gammat_new;
      stop_id=1;
      x=x_new;
      a=a_new;
      u=u_new;
      normc=normc_new;
      I=I_new;
      h=h_new;
      H=H_new;
    }else
    {
     //  Rcout << "start 7 in constOptOrth_H"   << std::endl;
      double factor;
      if(H_new>2*H)
      {
        if(H_new>4*H)
        {
          factor=0.125;
        }else
        {
          factor=0.25;
        }
      }else
      {
        factor=0.5;
      }
      double u1=factor;
      stop_id=0;
      int count3=0;
      while((stop_id==0)&(count3<3))
      {
        ++count3;
        gammat_tmp=gammat+u1*delta_gamma;
        fit_2=maxLinearH(gammat_tmp, nobs, nbasis, nvarX, tau, lambda);
        x_new=as<VectorXd>(fit_2(0));
        a_new=as<VectorXd>(fit_2(1));
        u_new=as<VectorXd>(fit_2(2));
        normc_new=as<double>(fit_2(3));
        I_new=as<VectorXi>(fit_2(4));

        h_new=normc_new*tau_sqrt*(orthConst_mtx*x_new);
        H_new=h_new.squaredNorm();

        if(H_new>H)
        {
          if(count3==3)
          {
            t=t+u1*t_inc;
            gammat=gammat_tmp;
            x=x_new;
            a=a_new;
            u=u_new;
            normc=normc_new;
            I=I_new;
            h=h_new;
            H=H_new;
          }else
          {
            u1=factor*u1;
          }
        }else
        {
          t=t+u1*t_inc;
          gammat=gammat_tmp;
          x=x_new;
          a=a_new;
          u=u_new;
          normc=normc_new;
          I=I_new;
          h=h_new;
          H=H_new;
          stop_id=1;

        }
      }
      // Rcout << "H_new= "<<" "<< H_new  << std::endl;
      // Rcout << "u1= "<<" "<< u1  << std::endl;
    }
    //    Rcout << "count1= "<<" "<< count1  << std::endl;
    //       Rcout << "H= "<<" "<< H  << std::endl;
  }
  return List::create(_["x"]=x, _["t"]=t, _["count1"]=count1, _["H"]=H);
}



////////////////////////////////////////////////////////////////////////


void find_orth_basis(Map<MatrixXd> &orthConst_mtx)
{
  int m=orthConst_mtx.rows();
  orthConst_mtx.row(0)=orthConst_mtx.row(0)/((orthConst_mtx.row(0)).norm());
  for(int i=1;i<m;++i)
  {
    VectorXd tmp=orthConst_mtx.bottomRows(m-i)*(orthConst_mtx.row(i-1).transpose());
    orthConst_mtx.bottomRows(m-i)=orthConst_mtx.bottomRows(m-i)-tmp*orthConst_mtx.row(i-1);
    orthConst_mtx.row(i)=orthConst_mtx.row(i)/((orthConst_mtx.row(i)).norm());
  }
}


////////////////////////////////////////////////////////////////////////

List get_comp(Map<MatrixXd> XbTransInv, MatrixXd Y_cent_ini, int K_comp, int max_comp,   int nbasis, int nvarX, double tau, double lambda, double thresh, double tol)
{

  MatrixXd Beta, Beta_return, B1;
  Map<MatrixXd> Y_cent=Map<MatrixXd>(Y_cent_ini.data(), Y_cent_ini.rows(), Y_cent_ini.cols());
  int totalbasis=nbasis*nvarX, upp_com;
  int nobs=XbTransInv.rows();
  MatrixXd new_orthConst_mtx(nobs, nobs+nbasis*nvarX), orthConst_mtx_final;
  new_orthConst_mtx.leftCols(nobs)=MatrixXd::Identity(nobs, nobs);
  new_orthConst_mtx.rightCols(totalbasis)=-XbTransInv;
  MatrixXd XY_prod_ini=Y_cent*XbTransInv;
  Map<MatrixXd> XY_prod=Map<MatrixXd>(XY_prod_ini.data(), XY_prod_ini.rows(), XY_prod_ini.cols());


  VectorXd tmpVXd, obj_val, obj_val_return, gamma(nobs+totalbasis), gamma_old, tmpvec,tempvec;
  if(K_comp==0)
  {
    Beta=MatrixXd::Zero(nbasis*nvarX, max_comp);
    obj_val=VectorXd::Zero(max_comp);
  }else
  {
    Beta=MatrixXd::Zero(nbasis*nvarX, K_comp);
    obj_val=VectorXd::Zero(K_comp);
  }
  int count=0;

  int ncomp=0;
  List Lambda_all(nvarX);
  MatrixXd tmp;
  while(ncomp<max_comp)
  {

    Map<MatrixXd> orthConst_mtx=Map<MatrixXd>(new_orthConst_mtx.data(), new_orthConst_mtx.rows(), new_orthConst_mtx.cols());
    if(ncomp==0)
    {
      find_orth_basis(orthConst_mtx);
      orthConst_mtx_final=orthConst_mtx;
      for(int i=0;i<nvarX; ++i)
      {
        tmp=orthConst_mtx.block(0, nobs+i*nbasis, orthConst_mtx.rows(), nbasis);
        Lambda_all(i)=tmp*tmp.transpose();
      }
    }
    tmp=orthConst_mtx.leftCols(nobs);
    B1=tmp*tmp.transpose();
    NumericVector xx = runif(gamma.size()); 
    gamma=as<Eigen::Map<Eigen::VectorXd> >(xx);

    gamma=gamma/gamma.norm();

    gamma_old=VectorXd::Zero(nobs+totalbasis);
    tmpvec=gamma.tail(gamma.size()-nobs);
    tempvec=VectorXd::Zero(nobs+totalbasis);
    tempvec.tail(totalbasis)=XY_prod.transpose()*(XY_prod*tmpvec);

    double err_1=(gamma.tail(totalbasis)-gamma_old.tail(totalbasis)).norm();
    double err_2=(gamma.tail(totalbasis)+gamma_old.tail(totalbasis)).norm();
    //  Rcout << "start 31 "  << std::endl;
    double err=min(err_1, err_2);
    VectorXd t0=MatrixXd::Zero(orthConst_mtx.rows(),1);
    while((count<50) & (err>gamma.tail(totalbasis).norm()*1e-5))
    {
      ++count;

      List ret=constOptOrth_H(tempvec, B1,  orthConst_mtx, t0,  Lambda_all,  nobs,  nbasis,  nvarX,  tau,  lambda, tol);
      gamma_old=gamma;
      t0=as<VectorXd>(ret(1));
      gamma=as<VectorXd>(ret(0));

      tmpvec=gamma.tail(totalbasis);
      tempvec=VectorXd::Zero(nobs+totalbasis);
      tempvec.tail(totalbasis)=XY_prod.transpose()*(XY_prod*tmpvec);
      err_1=(gamma.tail(totalbasis)-gamma_old.tail(totalbasis)).norm();
      err_2=(gamma.tail(totalbasis)+gamma_old.tail(totalbasis)).norm();

      err=min(err_1, err_2);

    }
    obj_val(ncomp)=tempvec.dot(gamma);
    // Rcout << "obj_val[ncomp]= "<<" "<< obj_val[ncomp]  << std::endl;
    // Rcout << "count "<<" "<< count  << std::endl;
    Beta.col(ncomp)=gamma.tail(totalbasis);
    if(K_comp==0)
    {
      if((ncomp==max_comp-1)||((obj_val(ncomp)/obj_val.head(ncomp+1).sum())<thresh))
      {
        Beta_return=Beta.leftCols(ncomp+1);
        upp_com=ncomp+1;
        obj_val_return=obj_val.head(upp_com);
        break;
      }
    } else
    {
      if(K_comp==ncomp+1)
      {
        Beta_return=Beta;
        upp_com=K_comp;
        obj_val_return=obj_val;
        break;
      }
    }
    ++ncomp;
    tmpVXd=VectorXd::Zero(nobs+totalbasis);
    tmpVXd.tail(totalbasis)=XbTransInv.transpose()*(XbTransInv*Beta.col(ncomp-1));
    tmpVXd=tmpVXd-orthConst_mtx.transpose()*(orthConst_mtx*tmpVXd);
    tmpVXd=tmpVXd/tmpVXd.norm();
    MatrixXd tmp(nobs+ncomp, nobs+ncomp);
    for(int i=0;i<nvarX; ++i)
    {
      // Rcout << "start= "  << std::endl;
      tmp.bottomRightCorner(nobs+ncomp-1, nobs+ncomp-1)=as<MatrixXd>(Lambda_all(i));
      tmp(0,0)=(tmpVXd.segment(nobs+i*nbasis, nbasis)).dot(tmpVXd.segment(nobs+i*nbasis, nbasis));
      tmp.bottomLeftCorner(nobs+ncomp-1,1)=orthConst_mtx.block(0, nobs+i*nbasis, orthConst_mtx.rows(), nbasis)*tmpVXd.segment(nobs+i*nbasis, nbasis);
      tmp.topRightCorner(1, nobs+ncomp-1)=(tmp.bottomLeftCorner(nobs+ncomp-1,1)).transpose();
      Lambda_all(i)=(tmp+tmp.transpose())/2;
      //    Rcout << "end= "  << std::endl;
    }
    MatrixXd old_orthConst_mtx=new_orthConst_mtx;
    new_orthConst_mtx=MatrixXd::Constant(nobs+ncomp, nobs+totalbasis, 0);
    new_orthConst_mtx.bottomRows(nobs+ncomp-1)=old_orthConst_mtx;
    new_orthConst_mtx.row(0)=tmpVXd.transpose();
    count=0;
  }


  return List::create(_["Beta"]=Beta_return, _["upp_comp"]=upp_com,  _["obj_val"]=obj_val_return);

}

////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////
VectorXd get_cv_error_hd_mf(MatrixXd T_train, MatrixXd T_valid, MatrixXd Y_train, MatrixXd Y_valid)
{
  int  ncol=T_train.cols();
  MatrixXd t_train_mtx(T_train.rows(), ncol+1);
  t_train_mtx.col(0)=MatrixXd::Constant(T_train.rows(),1, 1/sqrt(T_train.rows()));
  t_train_mtx.rightCols(ncol)=T_train;
  MatrixXd t_valid_mtx(T_valid.rows(), ncol+1);
  t_valid_mtx.col(0)=MatrixXd::Constant(T_valid.rows(),1, 1/sqrt(T_train.rows()));
  t_valid_mtx.rightCols(ncol)=T_valid;
  // Rcout << "t_train_mtx.transpose()*t_train_mtx=" << t_train_mtx.transpose()*t_train_mtx << std::endl;
  VectorXd error(ncol);
  MatrixXd E=Y_valid-t_valid_mtx.col(0)*(t_train_mtx.col(0).transpose()*Y_train);
  for(int ncomp=0;ncomp<ncol;ncomp++)
  {
    E.noalias()-=t_valid_mtx.col(ncomp+1)*(t_train_mtx.col(ncomp+1).transpose()*Y_train);
    error(ncomp)=E.squaredNorm()/Y_valid.cols();
  }
  return error;
}





// [[Rcpp::export]]
List cv_hd_msof(List X, Eigen::MatrixXd Y, List wb, List x_params,   List all_folds, int max_comp, double thresh)
{
  int nvarX=X.size(), nbasis=as<double>(x_params(4)), nobs=Y.rows();
  int nsample=Y.rows();
  double tol=1e-11, X_scale=0, X_scale_final;
  MatrixXd XbTransInv_final;
  int K_cv=all_folds.size();

  MatrixXd J=as<MatrixXd>(x_params(2)), J2=as<MatrixXd>(x_params(3)), Beta_final, normTransInv, normTransInv_final;
  VectorXd mu_y=Y.colwise().mean();
  MatrixXd tmpmat=Y.rowwise()-Y.colwise().mean();
  MatrixXd Y_cent=tmpmat.transpose();
  List mu_x(nvarX), Xb_lst(nvarX);
  X_scale=0;
  for(int i=0; i<nvarX;++i)
  {
    MatrixXd tmp_X=as<MatrixXd>(X(i));
    MatrixXd tmpmat=tmp_X.rowwise()-tmp_X.colwise().mean();
    mu_x(i)=tmp_X.colwise().mean();
    tmpmat=tmpmat*(as<MatrixXd>(wb(i)));
    for(int j=0; j<tmpmat.rows(); ++j)
    {
      X_scale=max(X_scale, tmpmat.row(j).norm());
    }
    Xb_lst(i)=tmpmat;
  }

  MatrixXd tau_lambda_set=as<MatrixXd>(x_params(6));
  VectorXd eta_set=as<VectorXd>(x_params(1));
  List  error_list(eta_set.size()*tau_lambda_set.rows());
  int opt_K;
  double opt_lambda, opt_tau, opt_eta, min_error=1e20;
  int tuning_id=0;
  for(int eta_id=0; eta_id<eta_set.size(); ++eta_id)
  {
    double eta=eta_set[eta_id];
    MatrixXd tmp=J+eta*J2;
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    MatrixXd R_ini=lltOf.matrixU();
    Map<MatrixXd> R_inv=Map<MatrixXd>(R_ini.data(), R_ini.rows(), R_ini.cols());
    normTransInv=R_inv.inverse();

    MatrixXd XbTransInv_ini=MatrixXd::Zero(nobs, nbasis*nvarX);
    for(int i=0; i<nvarX; i++)
    {
      XbTransInv_ini.block(0, i*nbasis, nobs, nbasis)=(as<MatrixXd>(Xb_lst(i)))*normTransInv;
    }
    XbTransInv_ini=XbTransInv_ini/X_scale;
    Map<MatrixXd> XbTransInv=Map<MatrixXd>(XbTransInv_ini.data(), XbTransInv_ini.rows(), XbTransInv_ini.cols());
    for(int tau_lambda_ind=0; tau_lambda_ind<tau_lambda_set.rows(); ++tau_lambda_ind)
    {
      double tau=tau_lambda_set(tau_lambda_ind, 0);
      double lambda=tau_lambda_set(tau_lambda_ind, 1);
      Rcout << "CV for eta, tau, lambda= "<<" "<<eta<<" "<<tau<<" "<<lambda   << std::endl;
      //Rcout << "start the diagnostic "   << std::endl;
      List fit_full=get_comp(XbTransInv, Y_cent, 0, max_comp,  nbasis, nvarX, tau,lambda, thresh, tol);
      VectorXi omit;


      int upp_comp=fit_full["upp_comp"];
      VectorXd errors=VectorXd::Zero(upp_comp);
      for(int fold_ind=0;fold_ind<K_cv;++fold_ind)
      {
//        Rcout << "CV fold"<<" "<< fold_ind+1 << std::endl;
        omit=as<VectorXi>(all_folds(fold_ind)).array()-1;
        int nsample_train=nsample-omit.size();
        MatrixXd Y_valid(omit.size(), Y.cols());
        MatrixXd Y_train(nsample_train, Y.cols());
        MatrixXd XbTransInv_train_ini(nsample_train, XbTransInv.cols()), XbTransInv_valid_ini(omit.size(), XbTransInv.cols());
        VectorXi ind=VectorXi::Zero(nsample);

        for(int i=0; i<omit.size(); i++)
        {
          ind(omit(i))=1;
          Y_valid.row(i)=Y.row(omit(i));
          XbTransInv_valid_ini.row(i)=XbTransInv.row(omit(i));
        }

        int i_train=0;
        for(int i=0; i<nsample; i++)
        {
          if(ind(i)==0)
          {
            Y_train.row(i_train)=Y.row(i);
            XbTransInv_train_ini.row(i_train)=XbTransInv.row(i);
            i_train=i_train+1;
          }
        }
     //   Rcout << "check 3"   << std::endl;
        Map<MatrixXd> XbTransInv_train=Map<MatrixXd>(XbTransInv_train_ini.data(), XbTransInv_train_ini.rows(), XbTransInv_train_ini.cols());
        Map<MatrixXd> XbTransInv_valid=Map<MatrixXd>(XbTransInv_valid_ini.data(), XbTransInv_valid_ini.rows(), XbTransInv_valid_ini.cols());
      //  Rcout << "check 4 "   << std::endl;
        tmpmat=Y_train.rowwise()-Y_train.colwise().mean();
        MatrixXd Y_train_cent=tmpmat.transpose();
        List fit_cv=get_comp(XbTransInv_train, Y_train_cent, upp_comp, max_comp, nbasis, nvarX, tau,lambda, thresh, tol);
        MatrixXd Beta_ini=as<MatrixXd>(fit_cv(0));
        Map<MatrixXd> Beta=Map<MatrixXd>(Beta_ini.data(), Beta_ini.rows(), Beta_ini.cols());
        MatrixXd T_train=XbTransInv_train*Beta;
        MatrixXd T_valid=XbTransInv_valid*Beta;
        for(int i=0; i<T_train.cols(); i++)
        {
          double tmp_value=(T_train.col(i)).norm();
          T_train.col(i)=T_train.col(i)/tmp_value;
          T_valid.col(i)=T_valid.col(i)/tmp_value;
        }
        errors=errors+get_cv_error_hd_mf(T_train, T_valid, Y_train, Y_valid);
      }
      error_list(tuning_id)=List::create(_["lambda"]=lambda, _["tau"]=tau, _["eta"]=eta, _["errors"]=errors);
      tuning_id++;
      MatrixXd::Index minIndex;
      double tmp_min=errors.minCoeff(&minIndex);
      //Rcout << "end 51 "  << std::endl;
      if(min_error>tmp_min)
      {
        min_error=tmp_min;
        opt_K=minIndex+1;
        opt_eta=eta;
        opt_tau=tau;
        opt_lambda=lambda;
        XbTransInv_final=XbTransInv;
        Beta_final=as<MatrixXd>(fit_full(0));
        normTransInv_final=normTransInv;
        X_scale_final=X_scale;
      }

    }
  }
  return List::create(_["XbTransInv"]=XbTransInv_final, _["X.scale"]=X_scale_final, _["Beta"]=Beta_final,  _["normTransInv"]=normTransInv_final,
                      _["mu.x"]=mu_x, _["errors"]=error_list, _["opt.K"]=opt_K,  _["opt.eta"]=opt_eta, _["opt.lambda"]=opt_lambda,
                        _["opt.tau"]=opt_tau, _["min.error"]=min_error);
}





// [[Rcpp::export]]
Eigen::MatrixXd get_pred_msof(List X_test, double X_scale, Eigen::MatrixXd XbTransInv, Eigen::MatrixXd Beta, Eigen::MatrixXd normTransInv, List mu_x, List wb, Eigen::MatrixXd Y,   int nobs, int nbasis, int nvarX)
{
  double tmp_value;
  MatrixXd XbTransInv_test_ini(nobs, nbasis*nvarX), tmpmat;

  List  Xb_lst(nvarX);

  for(int i=0; i<nvarX;++i)
  {
    tmpmat=as<MatrixXd>(X_test(i)).rowwise()- as<VectorXd>(mu_x(i)).transpose();
    Xb_lst(i)=tmpmat*(as<MatrixXd>(wb(i)));
    XbTransInv_test_ini.block(0, i*nbasis, nobs, nbasis)=(as<MatrixXd>(Xb_lst(i)))*normTransInv;
  }

  XbTransInv_test_ini=XbTransInv_test_ini/X_scale;
  Map<MatrixXd> XbTransInv_test=Map<MatrixXd>(XbTransInv_test_ini.data(), XbTransInv_test_ini.rows(), XbTransInv_test_ini.cols());

  MatrixXd T_test_ini=XbTransInv_test*Beta;
  MatrixXd T_ini=XbTransInv*Beta;
  MatrixXd T=Map<MatrixXd>(T_ini.data(), T_ini.rows(), T_ini.cols());
  MatrixXd T_test=Map<MatrixXd>(T_test_ini.data(), T_test_ini.rows(), T_test_ini.cols());

  for(int i=0; i<T.cols(); i++)
  {
    tmp_value=T.col(i).norm();
    T.col(i)=T.col(i)/tmp_value;
    T_test.col(i)=T_test.col(i)/tmp_value;
  }

  int ncol=T.cols();
  MatrixXd t_mtx(T.rows(), ncol+1);
  t_mtx.col(0)=MatrixXd::Constant(T.rows(),1, 1/sqrt(T.rows()));
  t_mtx.rightCols(ncol)=T;
  MatrixXd t_test_mtx(T_test.rows(), ncol+1);

  t_test_mtx.col(0)=MatrixXd::Constant(T_test.rows(),1, 1/sqrt(T.rows()));
  t_test_mtx.rightCols(ncol)=T_test;
  MatrixXd coef=t_mtx.transpose()*Y;
  MatrixXd Y_pred=t_test_mtx*coef;

  return Y_pred;
}

