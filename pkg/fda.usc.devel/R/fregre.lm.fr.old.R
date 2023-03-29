fregre.lm.fr.old <- function(formula,data,basis.y=NULL
                         ,basis.x=NULL,basis.b=NULL
                      #   ,lambda=0,P=c(1,0,0)
                         ,...){
  # res0 <-  fregre.lm.fr(as.formula(y~x),dat)
 # formula<-as.formula(y~x)
  # data<-mdat
  lambda=0
  P=c(1,0,0)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  datanames<-names(data)
  if (any(datanames=="df")) vnf=intersect(terms,names(data$df))
  else vnf <- NULL
  vfunc=setdiff(terms,vnf)
  off<-attr(tf,"offset")
  name.coef=nam=par.fregre=beta.l=list()
  ind.name<-NULL
  kterms=1
  isfdata<-FALSE
  XX<-list()
  lambdap<-XX2<-NULL
  RR<-NULL
  penalization<-FALSE    
  if (any(datanames==response)) {
# print("respuesta funcional")
    isfdata<-TRUE#is.fdata(data[[response]])
    yfdata <- data[[response]]
    ydat <-yfdata$data       # si es de la clase fd??
    #       print(paste("Functional response:",response))
    #XX[[response]]<- yfdata$data  
    tty<-yfdata[["argvals"]]
    rtty<-yfdata[["rangeval"]]
    namy<-yfdata[["names"]]
    n<-nrow(ydat)
    yy<-fdata.cen(data[[response]])
    #tty<-yy$argvals
    #ydat <- yy$data
    if (is.null(basis.y)){
      basis.y<-create.fdata.basis(yfdata,l=1:5)
      #     print(class(basis.y))
    } 
    ycoef <- ydat
    basis.y.class <- class(basis.y)
    if (basis.y.class=="basisfd"){
        #basis.y <- basis.x[[1]]
        #names(basis.x)
        fdnames=list("time"=tty,"reps"=rownames(ydat),"values"="values")
        # print("peta basisfd1")
        yfd= Data2fd(argvals = tty, y = t(ydat),
                     basisobj = basis.y,fdnames=fdnames,lambda=0)
        #  print("peta basisfd2")
        ycoef <- t(yfd$coefs)
      }
      # basis.y <- create.pc.basis(yfdata,1:4)
      # class(basis.pc)
    if (basis.y.class=="fdata.comp"){
      #  print("peta pc1")
        ycoef <- basis.y$x[, basis.y$l]
        #  print("peta pc2")
    }
    if (basis.y.class=="character"){
      #   print("peta raw1")
        ycoef <- ydat
        #    print("peta raw2")
    }
    #XX[[response]]<- yfdata$data  
    XX[[response]]<- ycoef
  }
  else { 
    if (any(names(data[["df"]])==response)) {
      stop(paste0("No functional response,",response," object must be fdata class"))
      yy<-as.matrix(data[["df"]][response])
      XX[[response]]<-yy    
    }
    else stop("Response not found")      
  }
  
  nvnf<-length(vnf)
  nn<-1
  if (nvnf>0) {
    XX2<-XX[[vnf]]<-as.matrix(data$df[,vnf])
    nn<-nvnf+1
    lambdap<-rep(0,nn)
  }
  else lambdap<-0
  lvfunc<-length(vfunc)
  #####if (length(lambda)==1) lambda<-rep(lambda,lvfunc)
  lambda<-rep(0,lvfunc)
  #print(paste("Functional covariate/s:",vfunc))
  if (lvfunc>0) { 
    mean.list=vs.list=JJ=list()
    bsp1<-bsp2<-TRUE
    for (i in 1:length(vfunc)) {
      if (inherits(data[[vfunc[i]]],"fdata")){
        
        tt<-data[[vfunc[i]]][["argvals"]]
        rtt<-data[[vfunc[i]]][["rangeval"]]
        fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
        if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
        else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
        if (is.null(basis.b[[vfunc[i]]])& bsp1)  basis.b[[vfunc[i]]]<-create.fdata.basis(fdat)
        else    if (inherits(basis.x[[vfunc[i]]],"fdata") | basis.x[[vfunc[i]]]$type=="pls") bsp2=FALSE
        if (bsp1 & bsp2) {
          if (is.null(rownames(dat)))    rownames(fdat$data)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rownames(fdat[["data"]]),"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]]=xcc[[2]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
            int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
            basis.x[[vfunc[i]]]$nbasis<-length(int)
            basis.x[[vfunc[i]]]$dropind<-NULL
            basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
          }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
            int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
            basis.b[[vfunc[i]]]$nbasis<-length(int)
            basis.b[[vfunc[i]]]$dropind<-NULL
            basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
          }
          #1 suavizado en los datos
          if (is.null(basis.x[[vfunc[i]]]$lambda)) lambda[i]<-3e-08/diff(rtt)
          else   lambda[i]<-basis.x[[vfunc[i]]]$lambda
          #print(lambda)
          x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames,
                         lambda=lambda[i])
          #2 penalizazion en la estimaxion X'X+lambdaRn
          r=x.fd[[2]][[3]]
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          colnames(J)=colnames(Z) = basis.b[[vfunc[i]]]$names
          name.coef[[vfunc[i]]]<-paste(vfunc[i],basis.b[[vfunc[i]]]$names,sep ="")
          
          ind.name<-c(ind.name,colnames(Z))
          XX[[vfunc[i]]] <- Z
          JJ[[vfunc[i]]] <- J
          #### new lambda pen
          #       if (lambda[i]==0)    basis.b[[vfunc[i]]]$lambda<-lambda[i]
          #print("aaaa                   bbbbbbbbbbbbbbbbb")
          if (is.null(basis.b[[vfunc[i]]]$lambda)) {
            lambdap<-c(lambdap,rep(0,len=basis.b[[vfunc[i]]]$nbasis))
          }
          else {
            penalization<-TRUE
            lambdap<-c(lambdap,rep(basis.b[[vfunc[i]]]$lambda,len=basis.b[[vfunc[i]]]$nbasis))
            if (!is.null(basis.b[[vfunc[i]]]$Lfdobj))        Lfdobj1<-basis.b[[vfunc[i]]]$Lfdobj[[vfunc[i]]]
            else Lfdobj1<-vec2Lfd(c(0,0),rtt)
            RR[[vfunc[i]]]<-eval.penalty(basis.b[[vfunc[i]]],Lfdobj1,rtt)
            XX2 = cbind(XX2, Z)
          }
          nn<-c(nn,basis.b[[vfunc[i]]]$nbasis)
        }
        else {
          pc<-basis.x[[vfunc[i]]]
          l<-basis.x[[vfunc[i]]]$l
          lenl<-length(l)
          vs <- t(pc$basis$data)
          Z<-pc$x[,l,drop=FALSE]
          cnames<-paste(vfunc[i],rownames(basis.x[[vfunc[i]]]$basis$data),sep ="")     
          colnames(Z)<-rownames(basis.x[[vfunc[i]]]$basis$data)
          name.coef[[vfunc[i]]]<-cnames
          if (lambda[i]==0)    lambdap<-rep(0,len=lenl)
          else {
            colnames(Z)<-cnames
            penalization<-TRUE  
            lambdap<-c(lambdap,rep(lambda[i],len=lenl))
          }
          # print(P)
          # print(l)
          # print("peta P.penalty1")
          if (length(l)>2) RR[[vfunc[i]]]<-P.penalty(l,P)# diag(lenl)
          else RR[[vfunc[i]]]<- diag(lenl)
          # print("peta P.penalty2")
          nn<-c(nn,lenl )
          XX[[vfunc[i]]] = Z
          vs.list[[vfunc[i]]]=pc$basis
          mean.list[[vfunc[i]]]=pc$mean
          
          # new ridge regression
          #[1] "fdataobj.pc"  "basis"        "x"            "mean"         "fdataobj.cen"
          #[6] "fdataobj"     "l"            "type"
          #    xmean<-pc$mean
          #print("cooooooolllllllllllllllllnnnnnnnnnnnnnnnaaaaaaames")   
          XX2 = cbind(XX2, Z)
          #########################
        }
      }
      else {
          if (inherits(data[[vfunc[i]]],"pca.fd"))
        {
          if (is.null(basis.x[[vfunc[i]]]))
          {
            basis.x[[vfunc[i]]]<-data[[vfunc[i]]]
            data[[vfunc[i]]]<-data[[vfunc[i]]]$harmonics
          }
        }
        if (inherits(data[[vfunc[i]]],"fd")){
          fdat<-data[[vfunc[i]]]
          if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
          else   if (inherits(basis.x[[vfunc[i]]],"pca.fd")) bsp1=FALSE
          if (is.null(basis.b[[vfunc[i]]])& bsp1)
            basis.b[[vfunc[i]]]<-create.fdata.basis(fdat,
                                                    l=1:max(5,floor(basis.x[[vfunc[i]]]$nbasis/5)),type.basis=basis.x[[vfunc[i]]]$type,
                                                    rangeval=fdat$basis$rangeval)
          else          if (inherits(basis.x[[vfunc[i]]],"pca.fd")) bsp2=FALSE
          if (bsp1 & bsp2) {
            r=fdat[[2]][[3]]
            if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
            }
            if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
            }
            J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
            mean.list[[vfunc[i]]]<-mean.fd(fdat)     ### changed
            x.fd<-center.fd(fdat)                    ### changed
            Z =t(x.fd$coefs) %*% J
            colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],basis.b[[vfunc[i]]]$names,sep="")
            XX[[vfunc[i]]] = Z
            JJ[[vfunc[i]]]<-J
            #### new lambda pen
            if (lambda[i]>0) basis.b[[vfunc[i]]]$lambda<-lambda[i]
            if (is.null(basis.b[[vfunc[i]]]$lambda)) {
              lambdap<-c(lambdap,rep(0,len=basis.b[[vfunc[i]]]$nbasis))
            }
            else {
              penalization<-TRUE
              lambdap<-c(lambdap,rep(basis.b[[vfunc[i]]]$lambda,len=basis.b[[vfunc[i]]]$nbasis))
              if (!is.null(basis.b[[vfunc[i]]]$Lfdobj))        Lfdobj1<-basis.b[[vfunc[i]]]$Lfdobj[[vfunc[i]]]
              else Lfdobj1<-vec2Lfd(P,rtt)
              RR[[vfunc[i]]]<-eval.penalty(basis.b[[vfunc[i]]],Lfdobj1,rtt)    
              XX2 = cbind(XX2, Z)
            }
            nn<-c(nn,basis.b[[vfunc[i]]]$nbasis)
          }
          else {
            # ridge regression    en pca basis
            pc<-basis.x[[vfunc[i]]]
            l<-ncol(pc$scores)
            vs <- pc$harmonics$coefs
            Z<-pc$scores
            colnames(Z)<-name.coef[[vfunc[i]]]<-paste(vfunc[i],colnames(pc$x),sep ="")
            XX[[vfunc[i]]] = Z
            vs.list[[vfunc[i]]]=vs
            mean.list[[vfunc[i]]]=pc$meanfd
            ######
            # new ridge regression
            if (lambda==0)    lambdap<-rep(0,len=l)
            else {
              penalization<-TRUE
              lambdap<-c(lambdap,rep(lambda[i],len=l))
            }
            RR[[vfunc[i]]]<-diag(l)
            nn<-c(nn,lenl)
            #print("new ridge regression")
            #print(nn)    
            XX2 = cbind(XX2, Z)
          }
        }
        else stop("Please, enter functional covariate")
      }
    }  }
  ####################################
  if (penalization){
    ####### FD BASIS ########
    # print("penalization")
    # print(nn)
    nsum<-cumsum(nn)
    nnsum<-sum(nn)
    R<-matrix(0,ncol=nnsum,nrow=nnsum)
    # print(dim(R))
    # print(lvfunc)  
    # print(nsum)    
    for (i in 1:lvfunc) {       
      R[(nsum[i]+1):nsum[i+1],(nsum[i]+1):nsum[i+1]]<- RR[[vfunc[i]]]
      #  print(RR[[vfunc[i]]])
    }
    # print(RR)           
    # print(lambdap)
    lambdap<-diag(lambdap)        
    XX2 =cbind(rep(1,len=n),XX2)###########
    colnames(XX2)[1]="(Intercept)"
    if (nvnf>0) colnames(XX2)[2:(1+nvnf)]=vnf
    Sb=t(XX2)%*%XX2+lambdap%*%R
    Lmat    <- chol(Sb)          # as fRegress
    Lmatinv <- solve(Lmat)
    Cinv<- Lmatinv %*% t(Lmatinv)
    Sb2=Cinv%*%t(XX2)
    DD<-t(XX2)%*%ydat
    S=XX2%*%Sb2
    yp=S%*%ydat
    b.est=Sb2%*%ydat
    bet<-Cinv%*%DD
    df=sum(nn)+1
    e=ydat-yp
    a.est=b.est[1,]
    z=list()
    rownames(b.est)<-colnames(XX2)
    z$coefficients<-b.est
    z$R<-R
    z$rank<-df
    z$df.residual<-n-df
    # print(vnf)       
    if (nvnf>0) colnames(XX2)[2:(1+nvnf)]=vnf
    # print(XX2[1,])
    # print(colnames(XX2)) 
    # print("oe")       
    coeff<-b.est[-1,,drop=F]
    if (isfdata) {       z$residuals<-fdata(e,tty,rtty,namy)
    z$fitted.values <- fdata(yp,tty,rtty,namy)
    }
    else {
      z$residuals <- e
      z$fitted.values <- yp
    }                      
  }
  else {
    par.fregre$data=XX
    ndatos<-nrow(XX[[response]])
    nx<-ncol(XX[[response]])
    
    # print(XX[[1]])
    print(class(XX))
    print(class(XX[[1]]))
    print(class(XX[[2]]))
     print(names(XX))
     print(formula)
    z=lm(formula=as.formula(formula),data=XX,...) ################################
    #print(z)   
    if (is.vector(z$coefficients)) {
      z$coefficients<-matrix(z$coefficients,ncol=1)
      colnames(z$coefficients)<- rownames(basis.y$basis$data)
      # print("PC coefs")
      # print(z$coefficients)
      # print(name.coef)
      rownames(z$coefficients)<- c("(Intercept)",unlist(name.coef))
      #rownames(z$coefficients)<- colnames(XX)
      # print(names(XX))
      
    }
    npy<-NCOL(z$coefficients)
    z$call<-z$call[1:2]
    df<-z$rank
    coeff<-z$coefficients

    if (basis.y.class=="basisfd"){
      z$yfit<-z$fitted.values  
      yhatfd <- fd(t(z$fitted.values), basis.y)
      z$fitted.values <- fdata(yhatfd,tty,rtty,namy)
      z$residuals<-yfdata-z$fitted.values
      #z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
    }
    if (basis.y.class=="fdata.comp"){
       z$yfit<-z$fitted.values  
       # yhatfdata <- z$fitted.values %*% (diag(length(basis.y$l))*basis.y$values[basis.y$l]) %*% basis.y$basis$data
       yhatfdata <- z$fitted.values %*% basis.y$basis$data
       # print(dim(yhatfdata))
       # print("fdata.comp")
       yhatfdata<-sweep(yhatfdata,2,matrix(basis.y$mean$data,ncol=1),"+")      
       z$fitted.values <- fdata(yhatfdata,tty,rtty,namy)
       z$residuals<-yfdata-z$fitted.values
      
      #z$residuals<-fdata(z$residuals,tty,rtty,namy)
      #z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
      
      #pc.fdata<-pc$u[,l,drop=FALSE]%*%(diag(lenl)*pc$d[l])%*%vs[l,,drop=FALSE]
      #pc.fdata<-sweep(pc.fdata,2,matrix(pc$mean$data,ncol=1),"+")
      
    }
    if (basis.y.class=="character"){
      z$yfit<-z$fitted.values  
      z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
      z$residuals<-fdata(z$residuals,tty,rtty,namy)
    }
    
    }
  
  for (i in 1:length(vfunc)) { 
    #      if (penalization)      name.coef2<-name.coef[[vfunc[i]]]
    #      else name.coef2<-paste(vfunc[i],name.coef[[vfunc[i]]],sep="")
    name.coef2<-name.coef[[vfunc[i]]]
    if (bsp1) {
      # print(bsp1)
      # print(name.coef2)
      # print(coeff)      
      if (isfdata) {  
        beta.l[[vfunc[i]]]=fd((coeff[name.coef2,]),basis.b[[vfunc[i]]])
        #rownames(coeff[name.coef2,])<-name.coef[[vfunc[i]]]
      }
      else {beta.l[[vfunc[i]]]=fd(z[["coefficients"]][name.coef2],basis.b[[vfunc[i]]])}                 
    }
    else{
      if (inherits(data[[vfunc[i]]],"fdata")){
        # print("calculo betas")
        # print(isfdata)  
        # print(name.coef2)
        # print("calculo betas2")
        # print(name.coef)    
        # print(coeff[,1:2])
        # print(isfdata)
      # ***************************************************************************
        if (isfdata) {beta.est<-drop(as.numeric(coeff[name.coef2,1]))*vs.list[[vfunc[i]]] }
        else {beta.est<-z$coefficients[name.coef2]*vs.list[[vfunc[i]]]        }    
        # print("calculo betas")
        
        # print(name.coef)
        
        beta.est$data<-apply(beta.est$data,2,sum)
        beta.est$names$main<-"beta.est"
        beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
        beta.l[[vfunc[i]]]<-beta.est
        #npy<-NCOL(beta.est)
        #print(beta.est)             
        if (npy>1) {
          #for (j in 2:npy) {
          for (j in 1:npy) {
            if (isfdata) beta.est<-drop(as.numeric(coeff[name.coef2,j]))*vs.list[[vfunc[i]]]
            beta.est$data<-apply(beta.est$data,2,sum)
            beta.est$names$main<-"beta.est"
            beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
            if  (basis.x[[vfunc[1]]]$type=="pls") {
              if (pc$norm)  {
                sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 2, var))
                beta.est$data<-  beta.est$data/sd.X
              }      
            }                 
            beta.l[[vfunc[i]]]<-c(beta.l[[vfunc[i]]],beta.est)
          }
        }
      }
      else {
        #           beta.l[[vfunc[i]]]=fd((coeff[name.coef2,]),basis.b[[vfunc[i]]])
        #            rownames(coeff[name.coef2,])<-name.coef[[vfunc[i]]]
        
        
        #            beta.est<-z$coefficients[name.coef[[vfunc[i]]],]*t(vs.list[[vfunc[i]]])
        #            beta.est<-apply(beta.est,2,sum)
        if (length(name.coef2)<basis.x[[vfunc[i]]]$harmonics$basis$nbasis)
          basis.x[[vfunc[i]]]$harmonics$basis$dropind<-(length(name.coef2)+1):basis.x[[vfunc[i]]]$harmonics$basis$nbasis
        beta.l[[vfunc[i]]]<-fd(z$coefficients[name.coef2,],basis.x[[vfunc[i]]]$harmonics$basis)
      }
    }
  }
  if (isfdata){  
    #sume<-colSums(z$residuals$data^2)}
  #else {
    #print(dim(z$residuals))
    #yhatmat = eval.fd(y$argvals, alphafd) %*% matrix(1, 1, ncurves)+   eval.fd(y$argvals,  xbetafd)
    #z$residuals = eval.fd(tty, z$residuals)
    #z$residuals<-fdata(t(z$residuals),tty,rtty,namy)
    #z$residuals <- fd(z$coefficients[name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
    #sume<-colSums(z$residuals$data^2)
  }
  
  #z$sr2=fdata(sume/(n-df),tty,rtty,namy)
  #z$r2<-fdata(1-sume/colSums(yy$Xcen$data^2),tty,rtty,namy)
  #}
  # else { #caso resp escalar
  #  z$sr2<-sum(z$residuals^2)/z$df.residual
  #  ###### z$Vp=z$sr2*S
  #}
  z$beta.l=beta.l
  z$formula=pf
  z$mean=mean.list
  z$formula.ini=formula
  z$basis.y=basis.y
  z$basis.x=basis.x
  z$basis.b=basis.b
  z$y<-data[[response]]
  z$JJ<-JJ
  z$data=z$data
  z$XX=XX
  z$vs.list=vs.list
  eY=z$y-mean(z$y)
  nullrss <- sum(norm.fdata(eY))
  rss <- sum(norm.fdata(z$residuals))
  ndf <- NROW(z$coefficients)
  z$sr2 <- rss/z$df.residual
  z$r2 = 1 - rss/nullrss
  z$response <- response
  # class(z)<-c(class(z),"fregre.fr")
  class(z)<-c("fregre.lm.fr",class(z))
  z
}

   

