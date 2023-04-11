#' @aliases predict.fregre.gsam.cl
#' @rdname predict.fregre.lm.fr
#' @noRd  
#' @export 
predict.fregre.gsam.cl <- function(object, newx = NULL,
                                type = "response",...){
 if (is.null(object)) stop("No fregre.gsam.dl object entered")
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
  yp<-object$fitted.values
  yp$data<-matrix(NA,nrow=nrow(XX),ncol=ncol(object$fitted.values))
  h <- object$h
#######################    
  for (i in 1:npy){
      # i<-1
      # cat(i)
      # print(dim(XX))
      # print(names(XX))
      if (h == 0 | i <= h)         { 
        for (j in 1:lenfunc){
          # print(j)
          # print(dim(data[[vfunc[j]]]$data))
          nam.var <- paste(vfunc[j],".t",i,sep="",collpase="")
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
          
          nam.var <- paste(vfunc[j],".t",i,sep="",collpase="")
          nam.var1 <- paste(vfunc[j],".t",i,"d",h,sep="",collpase="")
          
          XX[nam.var] <- data[[vfunc[j]]]$data[,i]
          XX[nam.var1] <- data[[vfunc[j]]]$data[,i]-data[[vfunc[j]]]$data[,i-h]
        }
      }
    # print(head(XX))
    # print(summary(object$result[[i]]))
      yp$data[,i]=predict.gam(object=object$result[[i]],
                              newdata=XX,type=type,...)
    }
    ##############################
  return(yp)

  # # print("raw")
  # npy <- NROW(object$basis$data)
  # yp <- matrix(NA,NROW(XX),npy)
  # for (i in 1:npy)  {
  #   yp[,i]=predict.gam(object=object$result[[i]],
  #                      newdata=XX,type=type,...)
  # }  
  # yp <- fdata(yp,tty,rtty,namy)
 }                          
return(yp)
}

