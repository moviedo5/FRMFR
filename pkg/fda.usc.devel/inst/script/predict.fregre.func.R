predict.fregre.func<-function(object,new.fdataobj=NULL,...){
if (is.null(object)) stop("No fregre.fd object entered")
#if (is.null(new.fdataobj)) stop("No newx entered")
if (is.null(new.fdataobj)) return(object$fitted.values)
if (!is.fdata(new.fdataobj)) new.fdataobj=fdata(new.fdataobj,object$fdataobj[["argvals"]],object$fdataobj[["rangeval"]],object$fdataobj[["names"]])
y=object$y
y.mat<-y$data
gg<-1:nrow(new.fdataobj)
nas<-apply(new.fdataobj$data,1,count.na)
if (any(nas)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   new.fdataobj$data<-new.fdataobj$data[bb,]
   gg<-gg[bb]
   }
newx<-new.fdataobj[["data"]]
tt<-new.fdataobj[["argvals"]]
rtt<-new.fdataobj[["rangeval"]]
nn <- nrow(new.fdataobj)
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn


 if (object$call[[1]]=="fregre.pc.func" | object$call[[1]]=="fregre.pc.func2"| object$call[[1]]=="fregre.pc.func3") {
 a1<-object$coefficients[1,]
# b2<-new.fdataobj$data%*%t(object$beta.est$data)
 beta.est<-fdata(object$beta.est$data,object$beta.est$argvals[[1]],object$beta.est$rangeval[[1]])
 b1<-inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],beta.est)#/(ncol(newx)-1)
 yp<-sweep(b1,2,a1,"+")

# b1<-inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$fdata.comp$rotation)#/(ncol(newx)-1)
# print(b1)
# print(a1)
# yp<-predict(object$lm,b1)
 
#glm           name.coef<-paste(vfunc[i], ".",rownames(object$basis.x[[vfunc[i]]]$basis$data),sep ="")
#glm          Z<- inprod.fdata(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]],object$vs.list[[vfunc[i]]])
#glm          colnames(Z)<-name.coef 
    
#          name.coef<-paste(vfunc[i], ".",rownames(object$basis.x[[vfunc[i]]]$basis$data),sep ="")
#          Z<- inprod.fdata(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]],object$vs.list[[vfunc[i]]])
#          colnames(Z)<-name.coef
#          npy<-length(object$result)
#          n<-nrow(XX[[i]])
#          yp<-matrix(NA,n,npy)
#          for (i in 1:npy) {
#            yp[,i]=predict.glm(object=object$result[[i]],newdata=XX,type=type,...)
#          }
           
  return(fdata(yp,y$argvals,y$rtt,y$names))
 }
 else {
  if (object$call[[1]]=="fregre.pls" ) {
  a1<-object$coefficients[1]*rep(1,len=nrow(newx))
  object$beta.est$data<-matrix(object$beta.est$data,nrow=1)
#  b2<-new.fdataobj$data%*%t(object$beta.est$data)
  b1<-inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$beta.est)/(ncol(newx)-1)
  yp<- a1+b1
 }
 else {
 if (object$call[[1]]=="fregre.basis.func" || object$call[[1]]=="fregre.basis.func.cv"){
  x=newx
  basis.x=object$basis.x.opt             #
  xcen<-fdata.cen(new.fdataobj,object$mean)[[1]]
	x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x)
  C=t(x.fd$coefs)
#  if (is.vector(object$b.est)) object$b.est<-matrix(object$b.est,ncol=1,nrow=length(object$b.est))
  yp=sweep(C%*%object$J%*%object$b.est,2,object$a.est,"+")
  return(fdata(yp,y$argvals,y$rtt,y$names))
  }
 else {
 if (object$call[[1]]=="fregre.np.func" || object$call[[1]]=="fregre.np.func.cv"){
 x=object$fdataobj
 h=object$h.opt
 n = nrow(x)
 nn = nrow(newx)
 np <- ncol(x)
# if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 par.S<-object$par.S
 bs<-as<-list()
 Ker=object$Ker
 #   par.metric<-list()
   par.metric<-attr(object$mdist,"par.metric")
   par.metric[["fdata1"]]<-new.fdataobj
   par.metric[["fdata2"]]<-x
#   parm<-attr(object$mdist,"par.metric")
#   lenpar<-length(parm)
#names(par.metric[3:(2+lenpar)])<-names(attr(object$mdist,"par.metric"))
  a1<-attr(object$mdist,"call")
  a2<-attr(object$par.S,"call")
  nmdist <- do.call(a1,par.metric)
  par.S$tt<-nmdist
  par.S$cv=FALSE
  H<-do.call(a2,par.S)
#print(y.mat)
  yp=H%*%y.mat
  return(fdata(yp,y$argvals,y$rtt,y$names))

 }     }
 }     }
#rownames(yp)=rownames(newx)
yp<-drop(yp)
names(yp)<-gg
return(yp)
}


