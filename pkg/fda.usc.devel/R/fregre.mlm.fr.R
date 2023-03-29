#' @aliases fregre.mlm.fr
#' @rdname fregre.lm.fr
#' @export
fregre.mlm.fr <- function(formula,data,basis.y=NULL
                         ,basis.x=NULL
                         ,lambda=NULL,P=NULL
                         ,...){
  # res0 <-  fregre.mlm.fr(as.formula(y~x),dat)
  # formula<-as.formula(y~x)
  # data<-mdat
  ##### lambda=0
  ##### P=c(1,0,0)
  # res2 <-  fregre.mlm.fr(as.formula(y~x),mdat,basis.y=basis.pc)
  # formula <- as.formula(y~x)
  # data <- mdat
  # basis.y=basis.pc
  # basis.x=NULL
  # lambda=NULL ;P=NULL
  
  beta2d <- list() # se guardan las beta funcionales en 2D (superficies)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  datanames<-names(data)
  if (any(datanames=="df")) 
    vnf=intersect(terms,names(data$df))  else vnf <- NULL
  vfunc=setdiff(terms,vnf)
  off<-attr(tf,"offset") 
  name.coef=nam=par.fregre=beta.l=list()
  ind.name<-NULL
  kterms <- 1
  isfdata <- FALSE
  XX <- list()
  lambdap <- XX2 <- NULL
  RR <- NULL
  penalization <- FALSE    
  # opciones para la respuesta (Bsp, pc o raw)
  if (any(datanames==response)) {
    # print("respuesta funcional")
    isfdata<-TRUE#is.fdata(data[[response]])
    yfdata <- data[[response]]
    ydat <- yfdata$data       # si es de la clase fd??
    #       print(paste("Functional response:",response))
    #XX[[response]]<- yfdata$data  
    tty <- yfdata[["argvals"]]
    rtty <- yfdata[["rangeval"]]
    namy <- yfdata[["names"]]
    n <- nrow(ydat)
    yy <- fdata.cen(data[[response]])
    #tty<-yy$argvals
    #ydat <- yy$data
    if (is.null(basis.y)){
      basis.y<-create.fdata.basis(yfdata,l=1:7)
      #     print(class(basis.y))
    } 
    ycoef <- ydat
    basis.y.class <- class(basis.y)
    
    if (basis.y.class == "basisfd"){
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
      ycoef <- basis.y$coefs[, basis.y$l]
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
  ####################  
  nvnf<-length(vnf)
  nn<-1
  if (nvnf>0) {
    XX2<-XX[[vnf]]<-as.matrix(data$df[,vnf])
    nn<-nvnf+1
    lambdap<-rep(0,nn)
  }
  else lambdap<-0
  lvfunc<-length(vfunc)
  # penalty <- TRUE
  # if (is.null(P) | is.null(lambda)){
  #   P <- as.list(numeric(lvfunc))
  #   names(P)<- vfunc
  #   lambda <- P
  #   penalty <- FALSE
  # }
  
  #####if (length(lambda)==1) lambda<-rep(lambda,lvfunc)
  #lambda<-rep(0,lvfunc)
  #print(paste("Functional covariate/s:",vfunc))
  if (lvfunc>0) { 
    mean.list<-vs.list<-list()
    
    out <- fdata2model.fr(vfunc = vfunc, vnf=vnf, response=response, XX=XX
                          , data=data, basis.x = basis.x,
                          pf=pf, tf=tf, lambda=lambda, P=P)
    
    XX <- out$XX
    bsp1 <- out$bsp1
    name.coef <- out$name.coef
    name.coef2 <- unlist(name.coef)
    basis.x <- out$basis.x
    vs.list <- out$vs.list
    mean.list <- out$mean.list
  }
  ####################################
  if (penalization){
    
    stop("not implemented yet")
    # z=list()
    # rownames(b.est)<-colnames(XX2)
    # z$coefficients<-b.est
    # #z$R<-R
    # z$rank<-df
    # z$df.residual<-n-df
    # if (nvnf>0) colnames(XX2)[2:(1+nvnf)]=vnf
    # coeff<-b.est[-1,,drop=F]
    # if (isfdata) {       z$residuals<-fdata(e,tty,rtty,namy)
    # z$fitted.values<-fdata(yp,tty,rtty,namy)
    # }    else {
    #   z$residuals<-e
    #   z$fitted.values<-yp
    # }                      
    }
  else { # MLM 
    #print("MLM")
    # par.fregre$data=XX
    # ndatos<-nrow(XX[[response]])
    # nx<-ncol(XX[[response]])
    z=lm(formula=as.formula(formula),data=XX,...) 
    if (is.vector(z$coefficients)) {
      z$coefficients<-matrix(z$coefficients,ncol=1)
      colnames(z$coefficients)<- paste(response,".",rownames(basis.y$basis$data),sep="")
      # print("PC coefs")
      # print(z$coefficients)
      # print(name.coef)
      rownames(z$coefficients)<- c("(Intercept)",unlist(name.coef))
      #rownames(z$coefficients)<- colnames(XX)
      # print(names(XX))
      
    }
    npy <- NCOL(z$coefficients)
    z$call <- z$call[1:2]
    df <- z$rank
    colnames(z$coefficients) <- paste(response,".",colnames(z$coefficients),sep="")
    rownames(z$coefficients)[-1] <- name.coef2
    coeff <- z$coefficients
    
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
    name.coef2 <- name.coef[[vfunc[i]]]
    if (bsp1) {
      if (isfdata) {  
        beta.l[[vfunc[i]]]=fd((coeff[name.coef[[vfunc[i]]],]),basis.x[[vfunc[i]]])
        }
      else {beta.l[[vfunc[i]]] <- fd(z[["coefficients"]][name.coef2],basis.x[[vfunc[i]]])}                 
    }
    else{
      
      #if (class(data[[vfunc[i]]])[1] == "fdata") {
        if (inherits(data[[vfunc[i]]],"fdata")) {            
        if (isfdata) {
          beta.est <- drop(as.numeric(coeff[name.coef2, 
                                            1])) * vs.list[[vfunc[i]]]
        }
        else {
          beta.est <- z$coefficients[name.coef2] * vs.list[[vfunc[i]]]
        }
        beta.est$data <- apply(beta.est$data, 2, sum)
        beta.est$names$main <- "beta.est"
        beta.est$data <- matrix(as.numeric(beta.est$data), 
                                nrow = 1)
        beta.l[[vfunc[i]]] <- beta.est
        if (npy > 1) {
          for (j in 2:npy) {
            if (isfdata) 
              beta.est <- drop(as.numeric(coeff[name.coef2, 
                                                j])) * vs.list[[vfunc[i]]]
            beta.est$data <- apply(beta.est$data, 2, 
                                   sum)
            beta.est$names$main <- "beta.est"
            beta.est$data <- matrix(as.numeric(beta.est$data), 
                                    nrow = 1)
            if (basis.x[[vfunc[1]]]$type == "pls") {
              if (basis.x[[vfunc[1]]]$norm) {
                sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 
                                   2, var))
                beta.est$data <- beta.est$data/sd.X
              }
            }
            beta.l[[vfunc[i]]] <- c(beta.l[[vfunc[i]]], 
                                    beta.est)
          }
        }
      }
      else {
        if (length(name.coef2)<basis.x[[vfunc[i]]]$harmonics$basis$nbasis)
          basis.x[[vfunc[i]]]$harmonics$basis$dropind<-(length(name.coef2)+1):basis.x[[vfunc[i]]]$harmonics$basis$nbasis
        beta.l[[vfunc[i]]]<-fd(z$coefficients[name.coef2,],basis.x[[vfunc[i]]]$harmonics$basis)
      }
    }
    
    #if (bsp1 & class(basis.y)=="basisfd") {  
    if (bsp1 & inherits(basis.y,"basisfd")) {            
      
      cc <- beta.l[[vfunc[i]]]$coefs
      bi <-bifd(cc,  basis.x[[vfunc[i]]], basis.y)
      dd <- eval.bifd(data[[vfunc[i]]]$argvals,data[[response]]$argvals,bi) # modificar grid.fdata para que haga lo mismo
      beta2d[[vfunc[i]]] <-  fdata(dd,list(data[[vfunc[i]]]$argvals,data[[response]]$argvals),fdata2d=TRUE)
    } else{
      beta2d <- NULL
      #beta2d[[vfunc[i]]] <-  fdata(beta.l[[vfunc[i]]]$data
      #      ,list(data[[vfunc[i]]]$argvals,data[[response]]$argvals),fdata2d=TRUE)
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
  z$beta.l <- beta.l
  z$formula <- pf
  z$mean <- mean.list
  
  z$formula.ini <- formula
  z$basis.y <- basis.y
  z$basis.x <- basis.x
  z$y<-data[[response]]
  #z$JJ<-JJ
  z$data <- z$data
  z$XX <- XX
  z$vs.list <- vs.list
  eY <- z$y-mean(z$y)
  nullrss <- sum(norm.fdata(eY)^2)
  rss <- sum(norm.fdata(z$residuals)^2)
  ndf <- NROW(z$coefficients)
  z$sr2 <- rss/z$df.residual
  z$r2  <-  1 - rss/nullrss
  z$response <- response
  # class(z)<-c(class(z),"fregre.fr")
  z$beta2d <- beta2d
  z$rn <- FALSE
  z$basis.y.class <- basis.y.class
  class(z)<-c("fregre.mlm.fr",class(z))
  z
}
