#' @aliases predict.fregre.gsam.con
#' @rdname predict.fregre.lm.fr
#' @noRd 
#' @export 
predict.fregre.gsam.con <- function(object, newx = NULL,
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
  
  # yfdata <- object$data[[response]]
  yfdata <- object$data[[response]]
  tty <- yfdata[["argvals"]]
  rtty <- yfdata[["rangeval"]]
  namy <- yfdata[["names"]]
  npy<-NCOL(yfdata)
  XX <- data.frame(y=1:NROW(data[[response]]))
  first=TRUE
  if (!is.null(vnf)) {
   first=FALSE
   #XX=NULL
   XX=cbind(XX,data[["df"]][,c(vnf)])
   names(XX)[-1]=vnf
  }
  lenfunc<-length(vfunc)
  # bsp1 <- object$bsp
  # bspy <- object$bspy # 0 es raw, 1 es pc y 2 es bsp/fou
  raw <- object$raw
  
  if (lenfunc>0) {
  k=1
  }
  
 
  if (!is.data.frame(XX)) 
    XX=data.frame(XX)
  if (is.null(object$basis.y)){
    yp<-object$fitted.values
    yp$data<-matrix(NA,nrow=nrow(XX),ncol=ncol(object$fitted.values))
    
    #for (i in 1:npy)  {
     # yp$data[,i]=predict.gam(object=object$result[[i]],newdata=XX,type=type,...)
    #}
    
    #print(head(XX))
    h <- object$h
    # print(h)
    # print("h")
#######################    
    for (i in 1:npy){
      # i<-1
      # cat(i)
      # print(dim(XX))
      # print(names(XX))
      if (h == 0)         { 
        for (j in 1:lenfunc){
          # print(j)
          # print(dim(data[[vfunc[j]]]$data))
          nam.var <- paste(vfunc[j],".h",i,sep="",collpase="")
          # print("ccc")
          XX[nam.var] <- data[[vfunc[j]]]$data[,i]
          # print("bbb")
          }
      } else{
        # hh <-1
        # print("hhhhhhhh")
        pf0 <-pf
        hh <- (i-h):(i+h)
        hh<-hh[hh>0 & hh<=npy]
        for (j in 1:lenfunc){
          #Xmat <- cbind(Xmat,data[[vfunc[j]]]$data[,hh])
          #nam.var <- paste(vfunc[j],".",i,"h",hh,sep="",collpase="")
          nam.var <- paste(vfunc[j],".h",i,sep="",collpase="")
          XX[nam.var] <- data[[vfunc[j]]]$data[,i]
          XX["delta"] <- factor(rep(0,len=NROW(XX)),levels=object$result[[1]]$xlevels$delta)
          
        }
      }
      yp$data[,i]=predict.gam(object=object$result[[i]],
                              newdata=XX,type=type,...)
    }
    ##############################
  return(yp)
  }
  # print("raw")
  npy <- NROW(object$basis$data)
  yp <- matrix(NA,NROW(XX),npy)
  for (i in 1:npy)  {
    yp[,i]=predict.gam(object=object$result[[i]],
                       newdata=XX,type=type,...)
  }  
  yp <- fdata(yp,tty,rtty,namy)
 }                          
return(yp)
}
# ########
# ind=1:129
# x=fdata.deriv(absorp[ind,])
# x1=fdata.deriv(x,1)
# y <- fdata.deriv(x1)
# ldat<-ldata(df=tecator$y[ind,],"y"=y,"x"=x,"x1"=x1)
# 
# x=fdata.deriv(absorp[-ind,])
# x1=fdata.deriv(x,1)
# newy <- fdata.deriv(x1)
# ldat2<-ldata(df=tecator$y[-ind,],"y"=newy,"x"=x,"x1"=x1)
# 
# res1=fregre.csam.fr(y~s(x),ldat,h=1)
# pred1 <- predict.fregre.csam.fr(res1,ldat2)
# sum(norm.fdata(pred1-newy))
# 
# res2=fregre.csam.fr(y~s(x),ldat,h=5)
# pred2 <- predict.fregre.csam.fr(res2,ldat2)
# sum(norm.fdata(pred2-newy))
# 
# res3=fregre.csam.fr(y~s(x),ldat,h=3)
# pred3 <- predict.fregre.csam.fr(res3,ldat2)
# sum(norm.fdata(pred3-newy))
# 
# bspx <- list("x"=create.fdata.basis(ldat$x,l=1:11))
# res.fr1=fregre.sam.fr(y~s(x),ldat,basis.x=bspx)
# pred.fr1 <- predict(res.fr1,ldat2)
# 
# b.x <- list("x"=create.pc.basis(ldat$x,1))
# res.fr2=fregre.sam.fr(y~s(x),ldat,basis.x=b.x)
# pred.fr2 <- predict(res.fr,ldat2)
# sum(norm.fdata(pred.fr1-newy))
# sum(norm.fdata(pred.fr2-newy))
# 
# b.y <- create.pc.basis(ldat$y,1:11)
# b.x <- list("x"=create.pc.basis(ldat$x,1:8))
# res.fr2=fregre.sam.fr(y~s(x),ldat,basis.x=b.x,basis.y=b.y)
# summary(res.fr2$result[[1]])
# pred.fr2 <- predict(res.fr2,ldat2)
# sum(norm.fdata(pred.fr1-newy))
# sum(norm.fdata(pred.fr2-newy))
# 
# # cambio en el num de bases la base no parece afectar 
# # res1=fregre.csam.fr(x.d2~+s(x),ldat,h=1,family=gaussian())
# #pred1 <- predict.fregre.csam.fr(res1,ldat)
# #plot(res1$fitted.values-pred1)
# # res1=fregre.csam.fr(x.d2~+Protein+s(x),ldat,h=1,family=gaussian())
# # pred1 <- predict.fregre.csam.fr(res1,ldat)
# # plot(res1$fitted.values-pred1)
# res0 <-  fregre.csam.fr(as.formula(y~s(x2)),mdat)
# pred0 <- predict.fregre.csam.fr(res0,mdat[1:111,row=T])
# plot(res0$fitted.values[1:111,]-pred0)
# 
# res2=fregre.csam.fr(x.d2~+s(x),ldat,h=4,family=gaussian())
# pred2 <- predict.fregre.csam.fr(res2,ldat[33:44,row=T])
# plot(res2$fitted.values[33:44,]-pred2)
# 
# ## Not run: 
# library(fda.usc.devel)
# data(tecator)
# absorp=tecator$absorp.fdata
# ind=1:129
# x=fdata.deriv(absorp[ind,])
# y=fdata.deriv(x,1)
# x2 <- fdata.deriv(y)
# mdat<-mfdata("y"=y,"x"=x,"x2"=x2)
# mdat<-ldata("y"=y,"x"=x,"x2"=x2)
# res0 <-  fregre.lm.fr(as.formula(y~x+x2),mdat)
# pred0 <- predict(res0,mdat)
# plot(res0$fitted.values-pred0)
# plot(res0$fitted.values-pred0)
# 
# res0 <-  fregre.sam.fr(as.formula(y~x+x2),mdat) #casca el specials2
# res0 <-  fregre.sam.fr(as.formula(y~x+s(x2)),mdat)
# pred0 <- predict(res0,mdat)
# plot(res0$fitted.values-pred0)
# 
# 
# 
# res0 <-  fregre.csam.fr(as.formula(y~s(x2)),mdat)
# pred0 <- predict.fregre.csam.fr(res0,mdat)
# plot(res0$fitted.values-pred0)
# plot(res0$fitted.values-pred0)
