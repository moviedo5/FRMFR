
###################################################################
# 20210329
# No fd (sólo fdata class)
# no basis.b
# y puede ser:
# (basis.y.class=="basisfd") y se coge los coeficientes de la base
# (basis.y.class=="fdata.comp") y se cogen las puntuaciones
# (basis.y.class=="character") y se usan todos los punos
###################################################################
fdata2model.fr <- function(vfunc, vnf, response, data, XX,
                           # basis.y=NULL, # no se usa
                           basis.x = NULL, pf,tf
                           , lambda=NULL, P=NULL){
  #print("entra fdata2model.fr")
  #######################
  
  kterms <- 1
  vs.list = mean.list=name.coef=nam=beta.l=list()
  bsp1 <- TRUE
  #  XX[[response]] <- data[[response]] coeficientes
  if (length(vnf) > 0) {
    for ( i in 1:length(vnf)){
      print(paste("Non functional covariate:",vnf[i]))
      if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
      else pf <- paste(pf, vnf[i], sep = "")
      XX[[vnf[i]]] <- data[["df"]][,vnf[i]]
      kterms <- kterms + 1
    }
    if   (attr(tf,"intercept")==0) {
      pf<- paste(pf,-1,sep="")
    }
  }
  lpenalty <- ipenalty <- list()
  lenfunc <- length(vfunc)>0
  if (lenfunc) {
    for (i in 1:length(vfunc)) {
      xaux <- data[[vfunc[i]]]
      if (class(xaux)[1]=="fdata"){
        if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]] <- create.fdata.basis(xaux,l=1:5)    else  
          if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
          
          if (is.null(P[[vfunc[i]]])) P[[vfunc[i]]] <- 0
          
          if (is.null(lambda[[vfunc[i]]])) lambda[[vfunc[i]]] <- 0
          aux <- vfunc2aux(xaux, nam.func=vfunc[i], kterms=kterms, pf=pf,
                           basis.x=basis.x[[vfunc[i]]], #basis.b = basis.x[[vfunc[i]]],
                           lambda=lambda[[vfunc[i]]],P=P[[vfunc[i]]],bsp1=bsp1)
          kterms <- aux$kterms
          pf <- aux$pf
          vs.list[[vfunc[i]]] = aux$vs.list
          mean.list[[vfunc[i]]] = aux$mean.list
          XX[[vfunc[i]]] <- aux$Z # si es MLM
          #                  XX<-cbind(XX, Z)  # si es LM
          ##################basis.x=basis.x se puede poner arriba  y quitar de vfunc2aux
          # penalty=lambda0
          #bsp1=bsp1
          name.coef[[vfunc[i]]]=aux$name.coef
          lpenalty[[vfunc[i]]]=aux$lpenalty
          ipenalty[[vfunc[i]]]=aux$ipenalty 
      }      else {       stop("Option not implemented")      }
    }  }
  else pf <- tf
  pf <- as.formula(pf)
  #if (!is.data.frame(XX)) XX=data.frame(XX)
  #colnames(XX)[1] <-  response
  # print("sale fdata2model.fr")
  return(list(pf=pf, vs.list=vs.list, mean.list=mean.list, XX=XX,
              basis.x=basis.x
              , name.coef=name.coef, bsp1=bsp1
              ,lpenalty=lpenalty, ipenalty=ipenalty))
}


###################################
vfunc2aux <- function(x, nam.func, kterms=1, pf,
                      basis.x=NULL, 
                      lambda=NULL, P=NULL, bsp1){
  # print("vfunc2aux")
  tt<- x$argvals
  rtt<-x$rangeval
  fdat<-x
  dat<-x$data
  
  
  if (bsp1) {
    if (is.null(rownames(dat)))    rownames(fdat$data) <- 1:nrow(dat)
    fdnames <- list("time"=tt,"reps"=rownames(dat),"values"="values")
    xcc <- fdata.cen(x)
    mean.list <- xcc[[2]]
    # print("bsp11111111111111111111111")
    # print(bsp1);     print(mean.list)
    if (!is.null(basis.x$dropind)) {
      int <- setdiff(1:basis.x$nbasis,basis.x$dropind)
      basis.x$nbasis <- length(int)
      basis.x$dropind <- NULL
      basis.x$names <- basis.x$names[int]
    }
    x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x, fdnames=fdnames)
    r=x.fd[[2]][[3]]
    J = inprod(basis.x, basis.x) 
    Z = t(x.fd$coefs) %*% J
    name.coef <- paste(nam.func,".",basis.x$names,sep="")
    colnames(J) <- colnames(Z) <- name.coef
    #XX = cbind(XX,Z)
    for ( j in 1:length(colnames(Z))){
      if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
      else pf <- paste(pf, colnames(Z)[j], sep = "")
      kterms <- kterms + 1
      lenl<- NCOL(Z)
    }
    vs.list <- J
  }        else {    #PC o PLS
    #          print(PC)
    l <- basis.x$l
    lenl<-length(l)
    vs <- t(basis.x$basis$data)
    Z <- basis.x$coefs[,l,drop=FALSE]     
    colnames(Z) = name.coef <- paste(nam.func, ".",colnames(Z),sep ="")      
    #XX <- cbind(XX,Z)
    vs.list <- basis.x$basis
    mean.list <- basis.x$mean
    for ( j in 1:lenl){
      if (kterms >= 1)  pf <- paste(pf, "+", name.coef[j], sep = "")
      else pf <- paste(pf, name.coef[j], sep = "")
      kterms <- kterms + 1
    }       
  }  
  if (is.null(P) ) P <- 0
  if (is.null(lambda)) lambda <- 0
  lpenalty <- lambda * P.penalty(1:lenl,P) 
  ipenalty <- (kterms-lenl+1):kterms
  # print("sale vfunc2aux")  
  return(list(kterms=kterms, pf=pf,
              vs.list=vs.list, mean.list=mean.list, Z=Z,#XX=XX,
              name.coef=name.coef,
              bsp1=bsp1  ,lpenalty=lpenalty, ipenalty=ipenalty
              ,rn=FALSE))
}

