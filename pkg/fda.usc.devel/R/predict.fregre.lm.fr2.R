# repite el modelo para cada columna del y (Raw o basis)
predict.fregre.lm.fr.old <- function(object, newx = NULL,
                                  type = "response",...){
  if (is.null(object)) stop("No fregre.lm.fr object entered")
  #print(1)
  # object <-res
  # newx <- ldata
  # type <- "response"
  if (is.null(newx)) {
    if (type == "effects"){
      # fake  = predict.lm(object, type = "terms", ...) 
      # yp <- effect.gam(object,fake)
    } else{
      # yp  = predict.gam(object, type = type, ...)
      yp <- object$fitted.values
    }
    return(yp)
  } else {
    #print(2)
    
    data=newx
    basis.x=object$basis.x
    basis.y=object$basis.y
    formula=object$formula.ini
    tf <- terms.formula(formula)
    terms <- attr(tf, "term.labels")
    if (length(terms)==0) return(rep(object$coefficient,length=nrow(newx[[1]])) ) 
    nt <- length(terms)
    if (attr(tf, "response") > 0) {
      response <- as.character(attr(tf, "variables")[2])
    }
    vtab<-rownames(attr(tf,"factors"))
    name.coef <- NULL
    vfunc <- object$vfunc
    vnf <- object$vnf
    nnf<-length(vnf)
    
    if (!is.null(vnf)) {
      first=FALSE
      XX=data.frame(data[["df"]][,c(vnf)])
      names(XX)=vnf
    } else {  
      first=TRUE
    }
    lenfunc<-length(vfunc)
    bsp1 <- object$bsp
    raw <- object$raw
    
    if (lenfunc>0) {
      k=1
      mean.list=vs.list=JJ=list()
      for (i in 1:lenfunc) {
        if(class(newx[[vfunc[i]]])[1]=="fdata"){
          tt<-data[[vfunc[i]]][["argvals"]]
          rtt<-data[[vfunc[i]]][["rangeval"]]
          fdataobj<-data[[vfunc[i]]]
          fdat<-data[[vfunc[i]]];      dat<-fdataobj$data
          if (nrow(dat)==1) rwn<-NULL         else rwn<-rownames(dat)
          xaux <- fdata2basis(data[[vfunc[i]]],basis.x[[vfunc[i]]])
          name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
          Z <- xaux$coefs
          if (first) {
            XX=Z
            first=FALSE
          }   else {
            XX = cbind(XX,Z)} }
      }
    }
    
    yfdata <- object$data[[response]]
    tty <- yfdata[["argvals"]]
    rtty <- yfdata[["rangeval"]]
    namy <- yfdata[["names"]]
    npy<-NCOL(yfdata)
    if (!is.data.frame(XX)) 
      XX=data.frame(XX)
    if (is.null(object$basis.y)){
      yp<-object$fitted.values
      yp$data<-matrix(NA,nrow=nrow(XX),ncol=ncol(object$fitted.values))
      for (i in 1:npy)  {
        # print(i)
        yp$data[,i]=predict.lm(object=object$result[[i]],newdata=XX,type=type,...)
      }  
      return(yp)
    }
    if (bsp1) {
      # print(bsp1)
      npy <- basis.y$nbasis
      yp <- matrix(NA,NROW(XX),npy)
      for (i in 1:npy)  {
        yp[,i]=predict.lm(object=object$result[[i]],
                           newdata=XX,type=type,...)
      }  
      yp <- fd(t(yp),basis.y)
      yp <- fdata(yp,tty,rtty,namy)
      return(yp)
    } else {
      # print("raw")
      #if (!raw) #si es raw ya se ha devuelto
      npy <- NROW(object$basis$data)
      yp <- matrix(NA,NROW(XX),npy)
      for (i in 1:npy)  {
        yp[,i]=predict.lm(object=object$result[[i]],
                           newdata=XX,type=type,...)
      }  
      # print("aa")
      yp <- gridfdata(yp,object$basis)
    }                          
    return(yp)
  }
}

