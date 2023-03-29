# # versiones antiguas que podrian funcionar pero se llama npy veces al lm
# fregre.lm.fr.pruebas <- function(formula,data,basis.y=NULL
#                           ,basis.x=NULL
#                           ,lambda=NULL,P=NULL
#                           ,...){
#   # res0 <-  fregre.mlm.fr(as.formula(y~x),dat)
#   # formula<-as.formula(y~x)
#   # data<-mdat
#   ##### lambda=0
#   ##### P=c(1,0,0)
#   # res2 <-  fregre.mlm.fr(as.formula(y~x),mdat,basis.y=basis.pc)
#   # formula <- as.formula(y~x)
#   # data <- mdat
#   # basis.y=basis.pc
#   # basis.x=NULL
#   # lambda=NULL ;P=NULL
#   #if (is.null(basis.y)) raw=TRUE
#   #else raw=FALSE
#   beta2d <- list() # se guardan las beta funcionales en 2D (superficies)
#   tf <- terms.formula(formula)
#   terms <- attr(tf, "term.labels")
#   nt <- length(terms)
#   if (attr(tf, "response") > 0) {
#     response <- as.character(attr(tf, "variables")[2])
#     pf <- rf <- paste(response, "~", sep = "")
#   } else pf <- rf <- "~"
#   vtab<-rownames(attr(tf,"factors"))
#   datanames<-names(data)
#   if (any(datanames=="df")) 
#     vnf=intersect(terms,names(data$df))  else vnf <- NULL
#   vfunc=setdiff(terms,vnf)
#   off<-attr(tf,"offset") 
#   name.coef=nam=par.fregre=beta.l=list()
#   ind.name<-NULL
#   kterms <- 1
#   isfdata <- FALSE
#   XX <- list()
#   lambdap <- XX2 <- NULL
#   RR <- NULL
#   penalization <- FALSE    
#   # opciones para la respuesta (Bsp, pc o raw)
#   if (any(datanames==response)) {
#     # print("respuesta funcional")
#     isfdata<-TRUE#is.fdata(data[[response]])
#     yfdata <- data[[response]]
#     ydat <- yfdata$data       # si es de la clase fd??
#     #       print(paste("Functional response:",response))
#     #XX[[response]]<- yfdata$data  
#     tty <- yfdata[["argvals"]]
#     rtty <- yfdata[["rangeval"]]
#     namy <- yfdata[["names"]]
#     n <- nrow(ydat)
#     yy <- fdata.cen(data[[response]])
#     #tty<-yy$argvals
#     #ydat <- yy$data
#      if (is.null(basis.y)) basis.y.class <- "raw"
#      else basis.y.class <-class(basis.y)
#     # if (is.null(basis.y)){
#     #  basis.y<-create.fdata.basis(yfdata,l=1:7)
#       #     print(class(basis.y))
#     #} 
#     ycoef <- ydat
#     #basis.y.class <- class(basis.y)
#     
#     if (basis.y.class == "basisfd"){
#       #basis.y <- basis.x[[1]]
#       #names(basis.x)
#       fdnames=list("time"=tty,"reps"=rownames(ydat),"values"="values")
#       # print("peta basisfd1")
#       yfd= Data2fd(argvals = tty, y = t(ydat),
#                    basisobj = basis.y,fdnames=fdnames,lambda=0)
#       #  print("peta basisfd2")
#       ycoef <- t(yfd$coefs)
#     }
#     
#     
#     # basis.y <- create.pc.basis(yfdata,1:4)
#     # class(basis.pc)
#     if (basis.y.class=="fdata.comp"){
#       #  print("peta pc1")
#       ycoef <- basis.y$coefs[, basis.y$l]
#       #  print("peta pc2")
#     }
#     
#     if (basis.y.class=="character"){
#       #   print("peta raw1")
#       ycoef <- ydat
#       #    print("peta raw2")
#     }
#     #XX[[response]]<- yfdata$data  
#     XX[[response]]<- ycoef
#   }
#   else { 
#     if (any(names(data[["df"]])==response)) {
#       stop(paste0("No functional response,",response," object must be fdata class"))
#       yy<-as.matrix(data[["df"]][response])
#       XX[[response]]<-yy    
#     }
#     else stop("Response not found")      
#   }
#   print("aaaaaaaaaaaaaaaaa")
#   ####################  
#   nvnf<-length(vnf)
#   nn<-1
#   if (nvnf>0) {
#     XX2<-XX[[vnf]]<-as.matrix(data$df[,vnf])
#     nn<-nvnf+1
#     lambdap<-rep(0,nn)
#   }
#   else lambdap<-0
#   lvfunc<-length(vfunc)
#   # penalty <- TRUE
#   # if (is.null(P) | is.null(lambda)){
#   #   P <- as.list(numeric(lvfunc))
#   #   names(P)<- vfunc
#   #   lambda <- P
#   #   penalty <- FALSE
#   # }
#   
#   #####if (length(lambda)==1) lambda<-rep(lambda,lvfunc)
#   #lambda<-rep(0,lvfunc)
#   print(paste("Functional covariate/s:",vfunc))
#   if (lvfunc>0) { 
#     mean.list<-vs.list<-list()
#     print(names(XX))
#     print("antes 2model")
#     out <- fdata2model.fr(vfunc = vfunc, vnf=vnf, response=response, XX=XX
#                           , data=data, basis.x = basis.x, basis.y=basis.y, 
#                           pf=pf, tf=tf, lambda=lambda, P=P)
#     
#     XX <- out$XX
#     print(names(XX))
#     print("despues 2model")
#     bsp1 <- out$bsp1
#     name.coef <- out$name.coef
#     name.coef2 <- unlist(name.coef)
#     basis.x <- out$basis.x
#     vs.list <- out$vs.list
#     mean.list <- out$mean.list
#   }
#   ####################################
#   if (penalization){
#     
#     stop("not implemented yet")
#     # z=list()
#     # rownames(b.est)<-colnames(XX2)
#     # z$coefficients<-b.est
#     # #z$R<-R
#     # z$rank<-df
#     # z$df.residual<-n-df
#     # if (nvnf>0) colnames(XX2)[2:(1+nvnf)]=vnf
#     # coeff<-b.est[-1,,drop=F]
#     # if (isfdata) {       z$residuals<-fdata(e,tty,rtty,namy)
#     # z$fitted.values<-fdata(yp,tty,rtty,namy)
#     # }    else {
#     #   z$residuals<-e
#     #   z$fitted.values<-yp
#     # }                      
#   }
#   else { 
#     
#     #z=lm(formula=as.formula(formula),data=XX,...) 
#     print(3)
#     print(out$pf)
#     formula <- as.formula(out$pf)
#    
#     npy <- NCOL(ydat)
#     zfitted <- zresid <- matrix(NA,n,npy)
#     pcoef<-zcoef<-NULL
#     result<-list()
#     Xmat <- data.frame(ydat[,1,drop=F])
#     #Xmat <- NULL
#     if (length(vnf)>0) {
#       Xmat<-data.frame(Xmat,data$df[,vnf])
#       names(Xmat)<-c(response,vnf)        
#       #Xmat <- data$df[,vnf,drop=F]
#     }  else       names(Xmat) <- response                   
#     
#     print("names(Xmat)1")
#     print(names(Xmat))
#     print(class(Xmat))
#     
#     print("names(XX)2")
#     print(names(XX))
#     print(class(XX))
#     
#     if (is.list(XX)){
#       for (i in 1:length(XX)){
#         Xmat<-cbind(Xmat,XX[[i]])
#         }
#     } else   Xmat <- (cbind(Xmat,XX))
#     print("names(Xmat)3")
#     print(colnames(Xmat))
#     
# 
#     print("name.coef4")
#     print(name.coef)
#     
#     print("name.coef5")
#     print(name.coef2)
#     
#     
#     print("names(ydat)6")
#     print(colnames(ydat))
#     print(class(ydat))
#     print(dim(ydat))
#     print(response)
#     Xmat[,response] <- ydat[,1]
#     print("names(ydat)7")
#     print(formula)
#     
#     z=lm(formula=formula,data=Xmat)  
#     print("saleeeeeeeeeeeeeee")
#     zcoef <- z$coefficients
#     zfitted[,i] <-z$fitted.values
#     #  print(names(z))
#     for (i in 2:npy){      
#       #result[[i]] <- z
#       zcoef <- rbind(zcoef, qr.coef(z$qr,ydat[,i]))
#       zfitted[,i] <-qr.fitted(z$qr,ydat[,i], k = z$qr$rank)
#       # H[,i] <-  hatvalues(z)
#     }
#     #colnames(H) <- nam.y
#     print(rownames(zcoef))# <- nam.y
#     #  print(bsp1)
#     out <- list()
#     #print(raw)
#     if (!isfdata){
#     if (bsp1) {
#       print("bsp1")
#       print(bsp1)
#       print(dim(zfitted))
#       print(basis.y)
#       zfitted <- fd(t(zfitted),basis.y)
#       out$fitted.values <- fdata(zfitted,tty,rtty,namy)
#       # H <- fd(t(H),basis.y)
#       # H <- fdata(H,tty,rtty,namy)
#     } else {
#       if (inherits(basis.y,"fdata")
#         out$fitted.values <- gridfdata(zfitted,out$basis)
#       #else  out$fitted.values <- fdata(zfitted,tty,rtty,namy)
#     }}else  out$fitted.values <- fdata(zfitted,tty,rtty,namy)
#     
#     print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
#     
#     # if (is.vector(z$coefficients)) {
#     #   z$coefficients<-matrix(z$coefficients,ncol=1)
#     #   colnames(z$coefficients)<- paste(response,".",rownames(basis.y$basis$data),sep="")
#     #   # print("PC coefs")
#     #   # print(z$coefficients)
#     #   # print(name.coef)
#     #   rownames(z$coefficients)<- c("(Intercept)",unlist(name.coef))
#     #   #rownames(z$coefficients)<- colnames(XX)
#     #   # print(names(XX))
#     #   
#     # }
#     # npy <- NCOL(z$coefficients)
#     # z$call <- z$call[1:2]
#     # df <- z$rank
#     # colnames(z$coefficients) <- paste(response,".",colnames(z$coefficients),sep="")
#     # rownames(z$coefficients)[-1] <- name.coef2
#     # coeff <- z$coefficients
#     # 
#     # if (basis.y.class=="basisfd"){
#     #   z$yfit<-z$fitted.values  
#     #   yhatfd <- fd(t(z$fitted.values), basis.y)
#     #   z$fitted.values <- fdata(yhatfd,tty,rtty,namy)
#     #   z$residuals<-yfdata-z$fitted.values
#     #   #z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
#     # }
#     # if (basis.y.class=="fdata.comp"){
#     #   z$yfit<-z$fitted.values  
#     #   # yhatfdata <- z$fitted.values %*% (diag(length(basis.y$l))*basis.y$values[basis.y$l]) %*% basis.y$basis$data
#     #   yhatfdata <- z$fitted.values %*% basis.y$basis$data
#     #   # print(dim(yhatfdata))
#     #   # print("fdata.comp")
#     #   yhatfdata<-sweep(yhatfdata,2,matrix(basis.y$mean$data,ncol=1),"+")      
#     #   z$fitted.values <- fdata(yhatfdata,tty,rtty,namy)
#     #   z$residuals<-yfdata-z$fitted.values
#     #   
#     #   #z$residuals<-fdata(z$residuals,tty,rtty,namy)
#     #   #z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
#     #   
#     #   #pc.fdata<-pc$u[,l,drop=FALSE]%*%(diag(lenl)*pc$d[l])%*%vs[l,,drop=FALSE]
#     #   #pc.fdata<-sweep(pc.fdata,2,matrix(pc$mean$data,ncol=1),"+")
#     #   
#     # }
#     # if (basis.y.class=="character"){
#     #   z$yfit<-z$fitted.values  
#     #   z$fitted.values<-fdata(z$fitted.values,tty,rtty,namy)
#     #   z$residuals<-fdata(z$residuals,tty,rtty,namy)
#     # }
#     
#   }
#   print("bbbb")
#   print(isfdata)
#   
#   for (i in 1:length(vfunc)) { 
#     #      if (penalization)      name.coef2<-name.coef[[vfunc[i]]]
#     #      else name.coef2<-paste(vfunc[i],name.coef[[vfunc[i]]],sep="")
#     name.coef2 <- name.coef[[vfunc[i]]]
#     if (bsp1) {
#       if (isfdata) {  
#         print("cccccccccccccccccccccccccccccccc1")
#         print(dim(zcoef))
#         print(class(zcoef))
#         print( colnames(zcoef))
#         print( rownames(zcoef))
#         print(vfunc)
#         print(name.coef)
#         print(name.coef[[vfunc[i]]])
#         #print(zcoef[,name.coef[[vfunc[i]]]])
#         print("cccccccccccccccccccccccccccccccc2")
#         beta.l[[vfunc[i]]]=fd(t(zcoef[,name.coef[[vfunc[i]]],drop=F]),
#                               basis.x[[vfunc[i]]])
#         print("cccccccccccccccccccccccccccccccc3")
#         #plot(beta.l[[vfunc[i]]])
#         
#       }
#       else {beta.l[[vfunc[i]]] <- fd(z[["coefficients"]][name.coef2],
#                                      basis.x[[vfunc[i]]])}                 
#     }
#     else{
#       
#       if (inherits(data[[vfunc[i]]],"fdata"){
#         print("bbb2")
#         if (isfdata) {
#           print(dim(zcoef))
#           print(name.coef)
#           beta.est <- drop(as.numeric(zcoef[,name.coef2]))* vs.list[[vfunc[i]]]
#         }
#         else {
#           beta.est <- z$coefficients[name.coef2] * vs.list[[vfunc[i]]]
#         }
#         beta.est$data <- apply(beta.est$data, 2, sum)
#         beta.est$names$main <- "beta.est"
#         beta.est$data <- matrix(as.numeric(beta.est$data), 
#                                 nrow = 1)
#         beta.l[[vfunc[i]]] <- beta.est
#         if (npy > 1) {
#           for (j in 1:npy) {
#             if (isfdata) 
#               beta.est <- drop(as.numeric(zcoef[,name.coef2])) * vs.list[[vfunc[i]]]
#             beta.est$data <- apply(beta.est$data, 2, 
#                                    sum)
#             beta.est$names$main <- "beta.est"
#             beta.est$data <- matrix(as.numeric(beta.est$data), 
#                                     nrow = 1)
#             if (basis.x[[vfunc[1]]]$type == "pls") {
#               if (basis.x[[vfunc[1]]]$norm) {
#                 sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 
#                                    2, var))
#                 beta.est$data <- beta.est$data/sd.X
#               }
#             }
#             beta.l[[vfunc[i]]] <- c(beta.l[[vfunc[i]]], 
#                                     beta.est)
#           }
#         }
#       }
#       else {
#         if (length(name.coef2)<basis.x[[vfunc[i]]]$harmonics$basis$nbasis)
#           basis.x[[vfunc[i]]]$harmonics$basis$dropind<-(length(name.coef2)+1):basis.x[[vfunc[i]]]$harmonics$basis$nbasis
#         beta.l[[vfunc[i]]]<-fd(z$coefficients[name.coef2,],basis.x[[vfunc[i]]]$harmonics$basis)
#       }
#     }
#     print(class(basis.y))
#     if (bsp1 & class(basis.y)=="basisfd") {  
#       
#       print("biiiiifffffffff")
#       cc <- beta.l[[vfunc[i]]]$coefs
#       print("biiiiifffffffff1")
#       print(dim(cc))
#       print((basis.x[[vfunc[i]]]$nbasis))
#       print((basis.y$nbasis))
#       print(dim(zcoef[,name.coef[[vfunc[i]]],drop=F]))
#       # bi <- bifd(t(cc[,]),  basis.x[[vfunc[i]]], basis.y)
#       print("biiiiifffffffff2")
#      ## dd <- eval.bifd(data[[vfunc[i]]]$argvals,data[[response]]$argvals,bi) # modificar grid.fdata para que haga lo mismo
#       print("biiiiifffffffff3")
#       #beta2d[[vfunc[i]]] <-  fdata(dd,list(data[[vfunc[i]]]$argvals,data[[response]]$argvals),fdata2d=TRUE)
#       beta2d[[vfunc[i]]] <-NULL
#     } else{
#       beta2d <- NULL
#       print(data[[vfunc[i]]]$argvals)
#       print(data[[response]]$argvals)
#       print(dim(zcoef[,name.coef[[vfunc[i]]]]))
#       print(names(beta.l[[vfunc[i]]]))
#       #beta2d[[vfunc[i]]] <-  fdata(#zcoef[,name.coef[[vfunc[i]]]]
#        # (beta.l[[vfunc[i]]]$coefs)
#         #    ,list(data[[vfunc[i]]]$argvals,data[[response]]$argvals),fdata2d=TRUE)
#     }
#   }
#   if (isfdata){  
#     #sume<-colSums(z$residuals$data^2)}
#     #else {
#     #print(dim(z$residuals))
#     #yhatmat = eval.fd(y$argvals, alphafd) %*% matrix(1, 1, ncurves)+   eval.fd(y$argvals,  xbetafd)
#     #z$residuals = eval.fd(tty, z$residuals)
#     #z$residuals<-fdata(t(z$residuals),tty,rtty,namy)
#     #z$residuals <- fd(z$coefficients[name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
#     #sume<-colSums(z$residuals$data^2)
#   }
#   print("eeeeeeeeeeeeeeeeeeeeeeee")
#   #z$sr2=fdata(sume/(n-df),tty,rtty,namy)
#   #z$r2<-fdata(1-sume/colSums(yy$Xcen$data^2),tty,rtty,namy)
#   #}
#   # else { #caso resp escalar
#   #  z$sr2<-sum(z$residuals^2)/z$df.residual
#   #  ###### z$Vp=z$sr2*S
#   #}
#   z$H <-hatvalues(z)
#   z$beta.l <- beta.l
#   z$formula <- pf
#   z$mean <- mean.list
#   
#   z$formula.ini <- formula
#   z$basis.y <- basis.y
#   z$basis.x <- basis.x
#   z$y<-data[[response]]
#   #z$JJ<-JJ
#   z$data <- z$data
#   z$XX <- XX
#   z$vs.list <- vs.list
#   print("fffffffffffff")
#   eY <- z$y-mean(z$y)
#   nullrss <- sum(norm.fdata(eY)^2)
#   print(class(z$residuals))
#   z$fitted.values <- out$fitted.values 
#   z$residuals <- z$y- out$fitted.values
#   rss <- sum(norm.fdata(z$residuals)^2)
#   ndf <- NROW(z$coefficients)
#   z$sr2 <- rss/z$df.residual
#   z$r2  <-  1 - rss/nullrss
#   z$response <- response
#   # class(z)<-c(class(z),"fregre.fr")
#   z$beta2d <- beta2d
#   z$rn <- FALSE
#   
#   z$basis.y.class <- basis.y.class
#   class(z)<-c("fregre.lm.fr",class(z))
#   z
# }
# 
# 
# 
# predict.fregre.lm.fr.recu <- function(object, newx = NULL,
#                                       type = "response",...){
#   if (is.null(object)) stop("No fregre.lm.fr object entered")
#   #print(1)
#   # object <-res
#   # newx <- ldata
#   # type <- "response"
#   if (is.null(newx)) {
#     if (type == "effects"){
#       # fake  = predict.lm(object, type = "terms", ...) 
#       # yp <- effect.gam(object,fake)
#     } else{
#       # yp  = predict.gam(object, type = type, ...)
#       yp <- object$fitted.values
#     }
#     return(yp)
#   } else {
#     #print(2)
#     
#     data=newx
#     basis.x=object$basis.x
#     basis.y=object$basis.y
#     formula=object$formula.ini
#     tf <- terms.formula(formula)
#     terms <- attr(tf, "term.labels")
#     if (length(terms)==0) return(rep(object$coefficient,length=nrow(newx[[1]])) ) 
#     nt <- length(terms)
#     if (attr(tf, "response") > 0) {
#       response <- as.character(attr(tf, "variables")[2])
#     }
#     vtab<-rownames(attr(tf,"factors"))
#     name.coef <- NULL
#     vfunc <- object$vfunc
#     vnf <- object$vnf
#     nnf<-length(vnf)
#     
#     if (!is.null(vnf)) {
#       first=FALSE
#       XX=data.frame(data[["df"]][,c(vnf)])
#       names(XX)=vnf
#     } else {  
#       first=TRUE
#       XX <- NULL
#     }
#     lenfunc<-length(vfunc)
#     raw <- object$raw
#     
#     if (lenfunc>0) {
#       k=1
#       mean.list=vs.list=JJ=list()
#       for (i in 1:lenfunc) {
#         if(class(newx[[vfunc[i]]])[1]=="fdata"){
#           tt<-data[[vfunc[i]]][["argvals"]]
#           rtt<-data[[vfunc[i]]][["rangeval"]]
#           fdataobj<-data[[vfunc[i]]]
#           fdat<-data[[vfunc[i]]];      dat<-fdataobj$data
#           if (nrow(dat)==1) rwn<-NULL         else rwn<-rownames(dat)
#           xaux <- fdata2basis(data[[vfunc[i]]],basis.x[[vfunc[i]]])
#           name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
#           Z <- xaux$coefs
#           if (first) {
#             XX=Z
#             first=FALSE
#           }   else {
#             XX = cbind(XX,Z)} }
#       }
#     }
#     
#     yfdata <- object$data[[response]]
#     tty <- yfdata[["argvals"]]
#     rtty <- yfdata[["rangeval"]]
#     namy <- yfdata[["names"]]
#     npy<-NCOL(yfdata)
#     if (!is.data.frame(XX)) 
#       XX=data.frame(XX)
#     if (is.null(object$basis.y)){
#       print("raw")
#       yp<-object$fitted.values
#       yp$data<-matrix(NA,nrow=nrow(XX),ncol=ncol(object$fitted.values))
#       for (i in 1:npy)  {
#         # print(i)
#         object$result$coefficients <- object$coefficients[i,]
#         yp$data[,i]=predict.lm(object=object$result,newdata=XX,type=type,...)
#         # yp = object$a.est * rep(1, len = nn) + Z %*%  object$b.est
#       }  
#       return(yp)
#     }
#     if (bsp) {  #Esta variable no esta definida!!
#       print(bsp)
#       npy <- basis.y$nbasis
#       yp <- matrix(NA,NROW(XX),npy)
#       for (i in i:npy)  {
#         object$result$coefficients <- object$coefficients[,i]
#         yp[,i]=predict.lm(object=object$result,newdata=XX,type=type,...)
#       }  
#       yp <- fd(t(yp),basis.y)
#       yp <- fdata(yp,tty,rtty,namy)
#       return(yp)
#     } else {
#       print("raw2")
#       #if (!raw) #si es raw ya se ha devuelto
#       npy <- NROW(object$basis$data)
#       yp <- matrix(NA,NROW(XX),npy)
#       for (i in 1:npy)  {
#         yp[,i]=predict.lm(object=object$result[[i]],
#                           newdata=XX,type=type,...)
#       }  
#       # print("aa")
#       yp <- gridfdata(yp,object$basis)
#     }                          
#     return(yp)
#   }
# }
