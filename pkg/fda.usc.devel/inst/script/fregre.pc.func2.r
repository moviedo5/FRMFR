# Shuang Wu. and Hans-Georg M¡§uller 2011 Response-Adaptive Regression for Longitudinal Data
# Lifeng WANG, Hongzhe LI, and Jianhua Z. HUANG2008 Functional Response Variable selection
# Hang Tesis Weighted functional principal component regression
# Rong Chen Functional coe..cient autoregressive models estimation and tests of hypotheses
# Aguilera Test 2012 Penalized spline approaches for functional logit regression
####################################################################
####################################################################
# si se resume la respuesta: 
# scores resumidos no es valido para representar el beta
# ver kokoszka 2008
####################################################################
####################################################################
fregre.pc.func2=function (fdataobj, y, l =NULL,rn=0,l.y=NULL,all.argvals=TRUE,weights= rep(1,ny),...){
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    if (is.null(l))    {        l<-pc$l       }
    else if (length(l)>nrow(pc$rotation)) stop("Incorrect value for  argument l")
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    namesx<-fdataobj[["names"]]
   }
else {
 if (is.null(l)) l<- 1:3
 if (!is.fdata(fdataobj))    fdataobj = fdata(fdataobj)
  tt<-fdataobj[["argvals"]]
  x<-fdataobj[["data"]]
  namesx<-fdataobj[["names"]]
  pc<-fdata2pc(fdataobj,ncomp=max(l),...)
}
ny<-nrow(y);npy<-ncol(y)
if (class(y)=="pca.fd") {
   is.yfd=TRUE
   y3<-t(y$coefs)
   y<-y$socres
   }
else {
  if (is.null(l.y)) {
      is.yfd=FALSE
      if (!is.fdata(y)) y=fdata(y)
      y3<-y[["data"]]
      tty2<-tty<-y[["argvals"]]
      rtty<-y[["rangeval"]]
      namesy<-y[["names"]]
      l.y<-1:npy
   }
   else  {
print(l.y)   
      is.yfd=TRUE
      tty2<-tty<-y[["argvals"]]
      rtty<-y[["rangeval"]]
      namesy<-y[["names"]]
if (all.argvals) { 
     y.fd<-create.pc.basis(y,l.y)$fdataobj.pc   
     y3<-y.fd$data # suavizados
     l.y<-1:ncol(y3)
     }
else {  
      y.fd<-fdata2pc(y,ncomp=max(l.y))
      y3<-y.fd$x[,l.y] #scores resumidos no es válido para represntar el beta??
      tty2<-seq(rtty[1],rtty[2],len=length(l.y))
}     
print(dim(y3))      
}
}
print(l.y)
len.l.y<-length(l.y)
#f (rn==TRUE) rn<-0.05*(pc$newd[1]^2)
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); np <- ncol(x);lenl = length(l)
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    X <-xcen<-pc$fdataobj.cen 
    if (n != ny)   stop("ERROR IN THE DATA DIMENSIONS")
   namy1<-namy2<-namesy
   namy1$main<-"Coefficient of determination"
   namy2$main<-"Residual variance"
   namy1$ylab<-expression("R"^2)
   namy2$ylab<-expression("sr"^2)
   C <- match.call()
#    if (length(l) == 1)     vs<-pc$rotation$data
#    else                    vs<-t(pc$rotation$data)
   scores<-Z<-(pc$x[,l])
   cnames<-colnames(pc$x)[l]
   df<-lenl+1
   J<-min(np,lenl)
   yy<-fdata.cen(y)
   ymean<-yy$meanX$data
   ycen<-yy$Xcen$data
   W<-diag(weights)   ######
if (rn>=0) {
#print("rn")
    xmean<-pc$mean
    d<-pc$newd[l]
    D<-diag(d)
    diagJ<-diag(J)
     scores<-cbind(rep(1,n),pc$x[,l])
     mat<-rn*diag(J+1)
     mat[1,1]<-0          
#     W<-diag(n)     
print(dim(W))     
print(dim(scores))     
     S<-solve(t(scores)%*%W%*%scores)       #incluir pesos solve(W)
#     S<-solve(t(scores)%*%W%*%scores+mat)       #incluir pesos solve(W)
print("W y S")           
print(W[1:3,1:4])
print(S[1:3,1:4])
print(class(y))     
print(dim(y))     
     Cinv<-S%*%t(scores)%*%W                    #incluir pesos W repetri proceso hasta que no cambie la prediccion   
     coeff<-Cinv%*%y3
     yp<-drop(scores%*%coeff)
print(class(yp))     
print(dim(yp))
     rownames(yp)<-rownames(y3)
     yp<-fdata(yp,tty,rtty,namesy)
     e<-y-yp                  
     H<-scores%*%Cinv
     df<-sum(diag(H))
     rdf<-n-df
     coeff<-drop(coeff)
     beta.est<-drop(as.numeric(coeff[-1,1]))*pc$rotation
     beta.est$data<-colSums(beta.est$data)
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)     
     beta.l<-beta.est  
         for (j in 2:(len.l.y)) {
         beta.est<-drop(as.numeric(coeff[-1,j]))*pc$rotation
         beta.est$data<-colSums(beta.est$data)
         beta.est$names$main<-"beta.est"
         beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
         beta.l<-c(beta.l,beta.est)     
         } 
 error<-TRUE
 it<-1
 eps=.001
 err2=sum(norm.fdata(beta.l))      
################################################################################     
 while (error) {   
     W0<-diag(diag(sweep(inprod.fdata(e),1,norm.fdata(e)^2,"/")     ))    
          W0<-diag(diag(inprod.fdata(e)  ))    
#     W0<-sweep(inprod.fdata(e),1,norm.fdata(e),"/")    
#     W0<-diag(diag(cov(t(e$data))))                   #solve scores no run
#     W0<-diag(diag(inprod.fdata(e) ))     #solve scores no run
print("WWWW")           
print(W0[1:3,1:4])     
     W <- try(solve(W0),silent=TRUE)
     if (class(W)=="try-error") {
     sv<-svd(W0)
print(dim(sv$v))
     W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
     warning("Inverse of sigma computed by SVD")
    }
print(it)     
       S<-solve(t(scores)%*%W%*%scores+mat)       
print("W y S")           
print(W[1:3,1:4])
print(S[1:3,1:4])
       Cinv<-S%*%t(scores)%*%W                                   #o es sqrt
       coeff<-Cinv%*%y3
       yp<-drop(scores%*%coeff)
       rownames(yp)<-rownames(y3)
       yp<-fdata(yp,tty,rtty,namesy)
print("residuals")              
       e<-y-yp                  
print("residuals2")              
#print(e$data[1:2,1:3])             
       coeff<-drop(coeff)  
#print("coefs")              
#print(coeff[1:2,1:3])        
       beta.est<-drop(as.numeric(coeff[-1,1]))*pc$rotation      
       beta.est$data<-colSums(beta.est$data)    
       beta.est$names$main<-"beta.est"
       beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1) 
       beta.l2<-beta.est              
#print(dim(beta.l))
#print(dim(beta.est))
       for (j in 2:(len.l.y)) {
         beta.est<-drop(as.numeric(coeff[-1,j]))*pc$rotation
         beta.est$data<-colSums(beta.est$data)
         beta.est$names$main<-"beta.est"
         beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
         beta.l2<-c(beta.l2,beta.est)     
         }     
print("aaaa2error") 
print(norm.fdata(beta.l2)[1:8,1])
       err3<-sum(norm.fdata(beta.l2))
       err4<-abs(err3-err2)
#       print(norm.fdata(beta.l-beta.l2))
print(c(err2,err3,err4))
       if (err4<=eps) error<-FALSE
       else {
        it<-it+1
        beta.l<-beta.l2
        if (it==100) error<-FALSE
        err2<-err3
        }
     }   
################################################################################     
     object.lm = list()
     object.lm$coefficients <- coeff
#    yp<-fdata(object.lm$fitted.values,tty,rtty,namesy)
    if (all.argvals) e<-y-yp
#    else e<-fdata(object.lm$residuals,tty,rtty,namesy)
#     e<-y-yp
     object.lm$y <- y
     object.lm$rank <- df
     object.lm$df.residual <-  rdf
     colnames(Z)[1] = "(Intercept)"
     class(object.lm) <- "lm"
#       GCV <- sum(e^2)/(n - df)^2

    beta.est<-drop(as.numeric(coeff[-1,1]))*pc$rotation
    beta.est$data<-colSums(beta.est$data)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
    beta.l<-beta.est  
#print(l.y)   
    if (len.l.y>1) {
    if (all.argvals){    
# print("all argvals")
    for (j in 2:(len.l.y)) {
         beta.est<-drop(as.numeric(coeff[-1,j]))*pc$rotation
         beta.est$data<-colSums(beta.est$data)
         beta.est$names$main<-"beta.est"
         beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
         beta.l<-c(beta.l,beta.est)     
         }
         }
  else{
       #print(dim(pc$rotation$data))
       #print(dim(coeff[-1,]))
       #print(dim(y3))
       #print(dim(y.fd$rotation$data))
##print(dim(beta.l)) 
#print("no all argvals")
xxx<-t(t(pc$rotation$data)%*%coeff[-1,]%*%y.fd$rotation$data)
#print(dim(xxx    ))
beta.l<-fdata(xxx,tt,rtt)    
#print("sale beta.est22")
tty2<-tty         
#print(ymean)
   yp<-fdata(sweep(inprod.fdata(fdataobj,beta.l),2,ymean,"+"),tty)   
   e<-y-yp
 }          
 }
#print("ok5")    
#######################################
#######################################    
   tt.list<-list("argvals.x"=tt,"argvals.y"=tty2)
   rtt.list<-list("range.x"=rtt,"range.y"=rtty)
   names.list<-list("main"="beta.est(s,t)","xlab"=c(namesx$xlab,namesy$xlab),"ylab"=c(namesy$xlab))
#print("ok6")
   beta.l<-fdata(t(beta.l$data),tt.list,rtt.list,names.list)
#print("ko")
#######################################
#######################################    
    sume<-colSums(e$data^2)
    sr2=fdata(sume/(n-df),tty,rtty,namy2)
    r2<-fdata(1-sume/colSums(yy$Xcen$data^2),tty,rtty,namy1)
    out <- list(call = C, beta.est = beta.l,coefficients=coeff,
    fitted.values =yp,residuals = e,H=H,df = df,r2=r2,
    sr2 = sr2,l = l,rn=rn,fdata.comp=pc,lm=object.lm,
    coefs=coefficients,fdataobj = fdataobj,y = y)
}
else {
    object.lm=lm(y3~Z,x=TRUE,y=TRUE,weights=weights,...)
    object.lm$call<-object.lm$call
    df<-object.lm$rank
    rdf<-n-df
    coeff<-object.lm$coefficients
    yp<-fdata(object.lm$fitted.values,tty,rtty,namesy)
    if (all.argvals) e<-y-yp
    else e<-fdata(object.lm$residuals,tty,rtty,namesy)
#       GCV <- sum(e^2)/(n - df)^2
    beta.est<-drop(as.numeric(coeff[-1,1]))*pc$rotation
    beta.est$data<-apply(beta.est$data,2,sum)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
    beta.l<-beta.est
#    persp(ttx,tt,beta_st)
# rehacer no esta bien este calculo    
    if (len.l.y>1) {
    if (all.argvals){
    for (j in 2:(len.l.y)) {
         beta.est<-drop(as.numeric(coeff[-1,j]))*pc$rotation
         beta.est$data<-colSums(beta.est$data)
         beta.est$names$main<-"beta.est"
         beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
         beta.l<-c(beta.l,beta.est)             
         }
         }
  else{
       #print(dim(pc$rotation$data))
       #print(dim(coeff[-1,]))
       #print(dim(y3))
       #print(dim(y.fd$rotation$data))
##print(dim(beta.l)) 
#print(dim(ycen))
xxx<-t(t(pc$rotation$data)%*%coeff[-1,]%*%y.fd$rotation$data)
#print(dim(xxx    ))
 beta.l<-fdata(xxx,tt,rtt)    
#print("sale beta.est22")
tty2<-tty         
#print(ymean)
   yp<-fdata(sweep(inprod.fdata(fdataobj,beta.l),2,ymean,"+"),tty)   
   e<-y-yp
 }      
 }
#######################################
#######################################    
   tt.list<-list("argvals.y"=tt,"argvals.x"=tty2)
   rtt.list<-list("range.y"=rtt,"range.x"=rtty)
   names.list<-list("main"="beta.est(s,t)","xlab"=c(namesx$xlab,namesy$xlab),"ylab"=c(namesy$xlab))
#print("peaAt")   
#print(length(tty))
#print(length(tt))
#print(length(tty2))
#print(dim(beta.l$data))
   beta.l<-fdata(t(beta.l$data),tt.list,rtt.list,names.list)
#######################################
#######################################          
    Z=cbind(rep(1,len=n),Z)
    S=solve(t(Z)%*%W%*%Z)
    H<-Z%*%S%*%t(Z)
    sume<-matrix(colSums(e$data^2),nrow=1)     
    sr2<-r2<-NULL    
  if (length(tty2)==length(tty) & all.argvals) {
     sr2=fdata(sume/(n-df),tty,rtty,namy2)
     r2<-fdata(1-sume/colSums(yy$Xcen$data^2),tty,rtty,namy1)
  }
  out <- list(call = C, beta.est = beta.l,coefficients=object.lm$coefficients,
  fitted.values =yp,residuals = e,H=H,
  df = df,r2=r2,sr2 = sr2,l = l,l.y=l.y,rn=rn,fdata.comp=pc,lm=object.lm,
  fdataobj = fdataobj,y = y)
  }
print(it)  
 class(out) = "fregre.func"
 return(out)
}
#res3<-fregre.pc.func2(ldata$y1,ldata$y0,l.y=1:4,l=1:4,all.argvals=TRUE)

#res.pc2<-fregre.pc.func2(temp,prec,l=l,l.y=l.y,all.argvals=FALSE) #no va
