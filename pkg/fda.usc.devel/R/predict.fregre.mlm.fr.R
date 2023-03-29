#' @aliases predict.fregre.mlm.fr
#' @rdname predict.fregre.lm.fr
#' @export 
predict.fregre.mlm.fr <- function (object, newx = NULL, type = "response",
                                   se.fit = FALSE, scale = NULL, df = NULL, interval = "none", 
                                   level = 0.95, weights = 1, pred.var = res.var/weights, ...) 
{
  
  if (is.null(object)) 
    stop("No fregre.lm object entered")
  if (is.null(newx)) {
    if (type == "effects"){
      fake  = predict.lm(object, type = "terms", se.fit = se.fit, 
                         interval = interval, level = level, weights = weights, 
                         pred.var = pred.var, df = df, scale = scale, ...)
      yp <- effect.fake(object,fake)
    } else{
      yp = predict.lm(object, type = type, se.fit = se.fit, 
                      interval = interval, level = level, weights = weights, 
                      pred.var = pred.var, df = df, scale = scale, ...)
    }
    return(yp)
  }
  else {
    
    data = newx
    yfdata <- object$y
    tty <- yfdata[["argvals"]]
    rtty <- yfdata[["rangeval"]]
    namy <- yfdata[["names"]]
    
    basis.x = object$basis.x
    basis.y = object$basis.y
    formula = object$formula.ini
    tf <- terms.formula(formula)
    terms <- attr(tf, "term.labels")
    nt <- length(terms)
    vtab <- rownames(attr(tf, "factors"))
    vnf = intersect(terms, names(data$df))
    vfunc2 = setdiff(terms, vnf)
    vint = setdiff(terms, vtab)
    vfunc = setdiff(vfunc2, vint)
    off <- attr(tf, "offset")
    beta.l = list()
    kterms = 1
    
    if (attr(tf, "response") > 0) {
      #response <- as.character(attr(tf, "variables")[2])
      response <- object$response
      pf <- rf <- paste(response, "~", sep = "")
    }     else pf <- rf <- "~"
    if (attr(tf, "intercept") == 0) {
      print("No intecept")
      pf <- paste(pf, -1, sep = "")
    }
    if (length(vnf) > 0) {
      first = FALSE
      for (i in 1:length(vnf)) {
        if (kterms > 1) 
          pf <- paste(pf, "+", vnf[i], sep = "")
        else pf <- paste(pf, vnf[i], sep = "")
        kterms <- kterms + 1
      }
      if (attr(tf, "intercept") == 0) {
        pf <- paste(pf, -1, sep = "")
      }
      mf <- as.data.frame(model.matrix(formula(pf), data$df))
      vnf2 <- names(mf)[-1]
      for (i in 1:length(vnf2)) pf <- paste(pf, "+", vnf2[i], 
                                            sep = "")
      XX <- mf
    } else {
      pf2 <- paste(pf, "1", sep = "")
      #XX <- data.frame(model.matrix(formula(pf2), data$df))
      #XX <- data.frame(model.matrix(formula(pf2), data[[response]]))
      XX <- NULL
      first = TRUE
    }
    
    if (length(vnf) > 0) {
      spm <- matrix(object$coefficients[names(XX)], ncol = 1)
      yp <- as.matrix(XX) %*% spm
    }
    else yp <- object$coefficients[1] * rep(1, len = nrow(newx[[vfunc[1]]]))
    XX <- list()
    if (length(vfunc) > 0) {
      for (i in 1:length(vfunc)) {
        # print(vfunc[i])
        if (class(data[[vfunc[i]]])[1] == "fdata") {
          fdataobj <- data[[vfunc[i]]]
          x.fd <- fdataobj[["data"]]
          tt <- fdataobj[["argvals"]]
          rtt <- fdataobj[["rangeval"]]
          if (!object$basis.x[[vfunc[i]]]$type == "pc" & 
              !object$basis.x[[vfunc[i]]]$type == "pls") {
            # print("Fixed basis")
            
            x.fd = Data2fd(argvals = tt, y = t(fdata.cen(fdataobj, 
                                                         object$mean[[vfunc[i]]])[[1]]$data),
                           basisobj = basis.x[[vfunc[i]]], 
                           fdnames = rownames(x.fd))
            # print("petaaaaa")
            r = x.fd[[2]][[3]]
            J <- object$vs.list[[vfunc[i]]]
            Z = t(x.fd$coefs) %*% J
            colnames(Z) = colnames(J)
          }          else {
            # print("PC o PLS")
            name.coef <- paste(vfunc[i], ".", rownames(object$basis.x[[vfunc[i]]]$basis$data), 
                               sep = "")
            newXcen <- fdata.cen(fdataobj, object$mean[[vfunc[i]]])[[1]]
            
            
            if (object$basis.x[[vfunc[i]]]$type == "pls") {
              if (object$basis.x[[vfunc[i]]]$norm) {
                sd.X <- sqrt(apply(object$basis.x[[vfunc[i]]]$fdataobj$data, 
                                   2, var))
                newXcen$data <- newXcen$data/(rep(1, 
                                                  nrow(newXcen)) %*% t(sd.X))
              }
            }
            Z <- inprod.fdata(newXcen, object$vs.list[[vfunc[i]]])
            colnames(Z) <- name.coef
          }
          if (first) {
            XX[[vfunc[i]]] = Z
            first = FALSE
          }
          else XX[[vfunc[i]]] = Z #cbind(XX, Z)
          colnames(XX[[vfunc[i]]])<-colnames(object$XX[[vfunc[i]]])
        }
      }
    }
    
    # representar en un base y y quedarnos con los coeficientes
    #predict.lm(object = object, data.frame(object$XX$x0))
    nn <- nrow(XX)
    # if (!is.data.frame(XX))       XX = data.frame(XX)
    if (!object$rn) {
      # print("no penalty")
      res.var <- object$sr2
      if (type == "effects"){
        fake  = predict.lm(object, newdata = XX, type = "terms", se.fit = se.fit, 
                           interval = interval, level = level, weights = weights, 
                           pred.var = pred.var, df = df, scale = scale, ...)
        return(effect.fake(object,fake))
      } else{
        class(object)<-class(object)[-1]
        aa<-predict(object = object, newdata = object$XX)
        yp <- predict(object = object, newdata = XX, 
                      type = type, se.fit = se.fit, interval = interval, 
                      level = level, weights = weights, pred.var = pred.var, 
                      df = df, scale = scale, , ...)
        ###pred <- fdata(pred,object$fitted.values$argvals,object$fitted.values$rangeval)
        
        basis.y.class <- object$basis.y.class
        if (basis.y.class=="basisfd"){
          yhatfd <- fd(t(yp), basis.y)
          yp <- fdata(yhatfd,tty,rtty,namy)
        }
        if (basis.y.class=="fdata.comp"){
          yhatfdata <- yp %*% basis.y$basis$data
          yhatfdata<-sweep(yhatfdata,2,matrix(basis.y$mean$data,ncol=1),"+")      
          yp <- fdata(yhatfdata,tty,rtty,namy)
        }
        if (basis.y.class=="character"){
          yp <- fdata(yp,tty,rtty,namy)
        }
        
        return(yp)
      }
    }
    else {
      if (type!="response") warning("Only response type implemented for penalization")
      for (i in 1:length(vfunc)) {
        if (object$call[[1]] == "fregre.pls") 
          return(predict.lm(object = object, newdata = XX, 
                            type = type, se.fit = se.fit, ...))
        if (object$basis.x[[vfunc[i]]]$type == "pc") {
          object$beta.l[[vfunc[i]]]$data <- matrix(object$beta.l[[vfunc[i]]]$data, 
                                                   nrow = 1)
          b1 <- inprod.fdata(fdata.cen(newx[[vfunc[i]]], 
                                       object$mean.list[[vfunc[i]]])[[1]], object$beta.l[[vfunc[i]]])
          yp <- yp + b1
        }
        else {
          xcen <- fdata.cen(newx[[vfunc[i]]], object$mean.list[[vfunc[i]]])[[1]]
          x.fd = Data2fd(argvals = xcen$argvals, y = t(xcen$data), 
                         basisobj = object$basis.x[[vfunc[i]]])
          C = t(x.fd$coefs)
          cnames <- colnames(object$vs.list[[vfunc[i]]])
          b.est <- matrix(object$coefficients[cnames], 
                          ncol = 1)
          b1 <- C %*% object$vs.list[[vfunc[i]]] %*% b.est
          yp <- yp + b1
        }
      }
      XX2 <- as.matrix(cbind(rep(1, len = nn), XX))
      predictor <- drop(yp)
      
      
      
      
      # if (se.fit || interval != "none") {
      #   ip <- rowSums((XX2 %*% object$Vp * XX2))
      #   res.var <- object$sr2
      #   df <- object$df.residual
      #   if (interval != "none") {
      #     tfrac <- qt((1 - level)/2, df)
      #     hwid <- tfrac * switch(interval, confidence = sqrt(ip), 
      #                            prediction = sqrt(ip + pred.var))
      #     predictor <- cbind(predictor, predictor + hwid %o% 
      #                          c(1, -1))
      #     colnames(predictor) <- c("fit", "lwr", "upr")
      #   }
      # }
      # if (se.fit) {
      #   se <- sqrt(ip)
      #   return(list(fit = predictor, se.fit = se, df = df, 
      #               residual.scale = sqrt(res.var)))
      # }
      # else return(predictor)
    }
  }
  return(drop(yp))
}


