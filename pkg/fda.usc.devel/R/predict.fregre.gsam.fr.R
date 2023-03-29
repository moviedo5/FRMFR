#' @aliases predict.fregre.sam.fr
#' @rdname predict.fregre.lm.fr
#' @export 
predict.fregre.sam.fr <- function(object, newx = NULL,
                                type = "response",...){
 if (is.null(object)) stop("No fregre.gsam.fr object entered")
  #print(1)
# object <-res
# newx <- ldata
# type <- "response"
  if (is.null(newx)) {
    if (type == "effects"){
      fake  = predict.gam(object, type = "terms", ...) 
      yp <- effect.gam(object,fake)
    } else{
      yp  = predict.gam(object, type = type, ...)
    }
    return(yp)
  } else {
    #print(2)
     
 data=newx
 basis.x=object$basis.x
 basis.y=object$basis.y
 formula=object$formula.ini
 tf <- terms.formula(formula, specials = c("s", "te", "t2"))
 terms <- attr(tf, "term.labels")
 if (length(terms)==0) return(rep(object$coefficient,length=nrow(newx[[1]])) ) 
 special <- attr(tf, "specials")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
 }
  vtab<-rownames(attr(tf,"factors"))
  gp <- interpret.gam(formula)
  #print(gp$smooth.spec)
  len.smo<-length(gp$smooth.spec)
  name.coef <- NULL
  vfunc <- object$vfunc
  vnf <- object$vnf
  nnf<-length(vnf)
  
  if (!is.null(vnf)) {
   first=FALSE
   #XX=NULL
   XX=data.frame(data[["df"]][,c(vnf)])
   names(XX)=vnf
  } else {  
    first=TRUE
  }
  lenfunc<-length(vfunc)
  bsp1 <- object$bsp
  bspy <- object$bspy # 0 es raw, 1 es pc y 2 es bsp/fou
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
    #  if (basis.x[[vfunc[i]]]$type=="pc" 
    #      | basis.x[[vfunc[i]]]$type=="pls")
    #    bsp1=FALSE      else bsp1 <- TRUE
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
      #print(i)
        yp$data[,i]=predict.gam(object=object$result[[i]],newdata=XX,type=type,...)
      }  
    return(yp)
  }
 if (bspy==2) {
  npy <- basis.y$nbasis
  yp <- matrix(NA,NROW(XX),npy)
  for (i in 1:npy)  {
    yp[,i]=predict.gam(object=object$result[[i]],
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
    yp[,i]=predict.gam(object=object$result[[i]],
                       newdata=XX,type=type,...)
  }  
  # print("aaaaaaaaaaaaa")
  # print(class(object$basis))
  # print(dim(yp))
  # print("aaaaaaaaaaaaa")
  # print(object$mean)
  # print(class(object$basis.y))
  #print(inherits(object$basis.y,"fdata"))
  if ( inherits(object$basis.y,"fdata.comp"))       
  yp <- gridfdata(yp,object$basis.y$basis,object$basis.y$mean)
  else   yp <- fdata(yp,tty,rtty,namy)
 }                          
#  if (type == "effects"){
#   fake  = predict.gam(object, newdata = XX, type = "terms",...)
#   yp <- effect.gam(object,fake)
#  } else{
#   yp <- predict.gam(object = object, newdata = XX, type = type,...)
# }                                                       
return(yp)
}
}

# reslin1=fregre.sam.fr(ffs, data=ltrain, basis.y=b.y, basis.x=b.x)
# prlin1=predict.fregre.sam.fr(reslin1,ltrain)

# pred.raw <- predict.fregre.sam.fr(res.pc,ldat)
# 
# res.pc=fregre.sam.fr(x.d2~+x+s(x.d1),data=ldat,#family=gaussian(),
#                      basis.x=xpc,
#                      basis.y=ypc
# )


# traceback()
# plot(res.raw$fitted.values)
# plot(pred.raw)
# 
# plot(pred.raw-res.raw$fitted.values)
# 
# pred.bsp <- predict.fregre.sam.fr(res.bsp,ldata)
# class(pred.bsp)
# plot(pred.bsp)
# plot(pred.bsp-res.bsp$fitted.values)
# 
# pred.pc <- predict.fregre.sam.fr(res.pc,ldata)
# plot(pred.pc)
# plot(pred.pc-res.bsp$fitted.values)


# data(aemet)
# basis.y<-create.pc.basis(aemet$logprec,1:5)
# res <- fregre.gsam.fr(logprec~s(temp),data=aemet,basis.y=basis.y)
# pred<-predict.fregre.sam.fr(res,aemet)
# plot(aemet$logprec[12])
# lines(res$fitted.values[12],col=3,lwd=2)
# lines(pred[12],col=6,lwd=2)
# 
# res <- fregre.gsam.fr(logprec~s(temp),data=aemet,basis.y=NULL)
# pred<-predict.fregre.sam.fr(res,aemet)
# lines(res$fitted.values[12],col=4,lwd=2)
# lines(pred[12],col=2,lwd=2)
# 
# basis.y<-create.bspline.basis(aemet$wind.speed$rangeval,81)
# res <- fregre.gsam.fr(logprec~s(temp),data=aemet,basis.y=basis.y)
# pred<-predict.fregre.sam.fr(res,aemet)
# lines(res$fitted.values[12],col=4,lwd=2)
# lines(pred[12],col=4,lwd=2)

# data(aemet)
# basis.y<-create.pls.basis(aemet$logprec,aemet$logprec,1:5)
# res <- fregre.gsam.fr(logprec~province,data=aemet,basis.y=basis.y)
# pred<-predict.fregre.sam.fr(res,aemet)
# no funciona 
# lines(res$fitted.values[12],col=3,lwd=2)
# lines(pred[12],col=6,lwd=2)



magic.post.proc <- function(X,object,w=NULL)
{ 
  V <- tcrossprod(object$rV)
  if (!is.null(w)) 
  { if (is.matrix(w)) WX <- X <- w%*%X else 
    WX <- as.vector(w)*X # use recycling rule to form diag(w)%*%X cheaply 
  } else {WX <- X}
  M <- WX%*%V  
  XWX <- crossprod(object$R) #t(X)%*%WX
  F <- Ve <- V%*%XWX
  edf1 <- rowSums(t(Ve)*Ve) 
  Ve <- Ve%*%V*object$scale 
  B <- X*M
  rm(M)
  hat <- rowSums(B) 
  H <- diag(X%*%V%*%t(WX))
  edf <- colSums(B) 
  Vb <- V*object$scale
  rm(V)
  list(Ve=Ve,Vb=Vb,hat=hat,edf=edf,edf1=2*edf-edf1,F=F,H=B)
}

