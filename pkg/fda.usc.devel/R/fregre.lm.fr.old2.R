# repite el modelo para cada columna del y (Raw o basis)
fregre.lm.fr.old2 <- function (formula
                             #, family = gaussian()
                             , data = list(), weights = NULL, 
                             basis.y = NULL, basis.x = NULL, ...) 
  {

    nam.data <- names(data)
    nam.df <- names(data$df)
    nam.func <- setdiff(nam.data,"df")
    # print(nam.data)  ;print(nam.df);print(nam.func)
    tf <- terms.formula(formula)
    terms <- attr(tf, "term.labels")
    nt <- length(terms)
    if (attr(tf, "response") > 0) {
      response <- as.character(attr(tf, "variables")[2])
      pf <- rf <- paste(response, "~", sep = "")
    }   else pf <-rf <- "~"
    vtab <- rownames(attr(tf, "fac"))
    nterms <- length(terms)
    fnf1 <- fnf2 <- fnf <- bs.dim1 <- bs.dim2 <- vfunc2 <- vfunc <- vnf <- NULL
    func <- nf <-  rep(0, nterms)
    names(func) <- names(nf) <- terms
    ndata <- length(data) - 1
    covar <- setdiff(nam.data,response)
    nam.df
    nam.func
    for (i in 1:length(covar)){
      if (covar[i] %in% nam.df) vnf<-c(vnf,covar[i])
      if (covar[i] %in% nam.func) vfunc<-c(vfunc,covar[i])
    }
    nnf<-length(vnf)
    nfunc<-length(vfunc)
    bsp1 <-  TRUE
    raw <- FALSE
    
    
    ######## manejando la respuesta  
    if (any(nam.data==response)) {
      isfdata<-is.fdata(data[[response]])
      if (isfdata) {
        #    print(paste("Functional response:",response))
        yfdata <-  data[[response]]
        
        if (is.null(basis.y)) {
          y <-yfdata$data
          bsp1 <- FALSE
          raw <- TRUE
          aux<-list("coefs"=NULL,"basis"=NULL)
          nam.y <- colnames(y)
          if (is.null(nam.y)) nam.y <-  yfdata$argvals
          
        }      else {
          aux <- fdata2basis(yfdata,basis.y) # si PCA que vaya centrada! meanY.list()
          y <- aux$coefs
          if (basis.y$type == "pc" | basis.y$type == "pls") 
            bsp1 <- FALSE
          nam.y <- colnames(y)
          #if (is.null(nam.y)) nam.y <-  yfdata$argvals
        }
        ndatos <- NROW(y)
        npy<-NCOL(y)
        
        #  XX[[response]] <- data[[response]][["data"]]
        tty <- yfdata[["argvals"]]
        rtty <- yfdata[["rangeval"]]
        namy <- yfdata[["names"]]
      }     else stop("Response must be of fdata class")
    }  else stop("Response not found in data object")     
    ########################################
    # fnf2<-rep(0,nterms)
    # fnf1 <- fnf <- nterms
    
    name.coef = nam = par.fregre = beta.l = list()
    kterms = 1
    if (nnf > 0) {
      XX = data[["df"]]
      if (attr(tf, "intercept") == 0) {
        pf <- paste(pf, -1, sep = "")
      }
      for (i in 1:nnf) {
        pf <- paste(pf, "+", vnf[i], sep = "")
        kterms <- kterms + 1
      }
    }   else {
      XX <- NULL
      #   XX = data.frame(data[["df"]][, response])
      #   names(XX) = response
    }
    lenfunc <- length(vfunc) 
    ifunc <- lenfunc > 0
    mean.list = basis.list = list()
    if (ifunc) {
      k = 1
      for (i in 1:lenfunc ) {
        if (is(data[[vfunc[i]]], "fdata")) {
          tt <- data[[vfunc[i]]][["argvals"]]
          rtt <- data[[vfunc[i]]][["rangeval"]]
          fdat <- data[[vfunc[i]]]
          nms <- data[[vfunc[i]]]$names
          #dat <- data[[vfunc[i]]]$data
          if (is.null(basis.x[[vfunc[i]]])) 
            basis.x[[vfunc[i]]] <- 
            create.fdata.basis(fdat, l = 1:7) else
              if (basis.x[[vfunc[i]]]$type == "pc" | basis.x[[vfunc[i]]]$type == "pls") 
                bsp1 = FALSE
          #         if (bsp1) {
          xaux <- fdata2basis(data[[vfunc[i]]],basis.x[[vfunc[i]]])
          name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
          Z <- xaux$coefs
          lencoef <- length(colnames(Z))
          XX = cbind(XX, Z)
          
          for (j in 1:lencoef) {
            pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
            kterms <- kterms + 1
          }       
          basis.list[[vfunc[i]]] <- xaux$basis
          # J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          #   vs.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$basis
          if (!bsp1 & !raw) 
            mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$mean
          else {
            xcc <- fdata.cen(data[[vfunc[i]]])
            mean.list[[vfunc[i]]] = xcc[[2]]
          }
        }      else {
          stop("Please, enter functional covariate of fdata class object")
        }
      }
    }
    #par.fregre$formula=as.formula(pf)
    #par.fregre$data=XX
    formula <- as.formula(pf)
    nx <- ncol(XX)
    Ymat<-y
    zfitted <- zresid <- matrix(NA,ndatos,npy)
    pcoef<-zcoef<-NULL
    result<-list()
    Xmat <- data.frame(y[,1,drop=F])
    if (length(vnf)>0) {
      Xmat<-data.frame(Xmat,data$df[,vnf])
      names(Xmat)<-c(response,vnf)         
    }  else       names(Xmat) <- response                   
    # for (i in 2:length(XX)){
    #    XX1<-data.frame(XX1,XX[[i]])
    #  }
    Xmat <- (cbind(Xmat,XX))
    # print(class(zfitted))
    
    for (i in 1:npy){
      Xmat[,response] <- Ymat[,i]
      z=lm(formula=formula,data=Xmat)  
      result[[i]] <- z
      
      # ss <- summary.sam(z)
      zcoef <- rbind(zcoef, z$coefficients)
      zfitted[,i] <-z$fitted.values
      # H[,i] <-  hatvalues(z)
    }
    #colnames(H) <- nam.y
    rownames(zcoef) <- nam.y
    #  print(bsp1)
    out <- list()
    if (bsp1) {
      zfitted <- fd(t(zfitted),basis.y)
      out$fitted.values <- fdata(zfitted,tty,rtty,namy)
      # H <- fd(t(H),basis.y)
      # H <- fdata(H,tty,rtty,namy)
    } else {
      if (inherits(basis.y,"fdata"))
        out$fitted.values <- gridfdata(zfitted,aux$basis)
      else  out$fitted.values <- fdata(zfitted,tty,rtty,namy)
    }
    out$fitted.values$names$main <- "Fitted values"
    #z$fitted.values <- zfitted #fdata(zfitted,tty,rtty,namy)
    out$residuals <- yfdata - out$fitted.values 
    out$residuals$names$main <- "Residuals"  
    out$coefficients <- zcoef
    out$edf.coef <- pcoef
    out$result<-result
    out$formula <- pf
    out$mean <- mean.list
    out$formula.ini=formula
    out$basis.x=basis.x
    out$basis.y=basis.y
    out$data=data
    out$basis <- aux$basis
    #out$basis.list <- basis.list
    out$raw <- raw
    out$bsp <- bsp1
    out$vfunc <- vfunc;   out$nnf <- nnf;   out$vnf <- vnf
    out$H <-hatvalues(z)
    class(out) <- c("fregre.lm.fr",class(out))
    # tabla con elementos significativos 
    # por fila los modelos y por columnas las mismas
    # componentes (edf si s p-valor por *)
    out
  }
  
