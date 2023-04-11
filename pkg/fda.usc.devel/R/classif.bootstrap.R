#' @title Functional Classification using Bootstrap
#' 
#' @description Computes functional classification using convex and local convex bootstrap for functional (and non functional)
#' explanatory variables 
#' 
#' @aliases classif.bootstrap classif.bootstrap2
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param data \code{list} that containing the variables in the model. The first item in the
#'  \code{data} list is called \emph{"df"} and is a data frame with the response and non 
#'  functional explanatory variables, as  \code{\link{glm}}.  Functional covariates of class 
#'  \code{fdata} or \code{fd} are introduced in' the following items in the \code{data} list.
#' @param B \code{integer}, number of boostrap resample,
#' @param weights Initial weights:   
#' \itemize{
#'  \item \code{character}, if \code{'equal'} (by default)  provides same weights for each observation; 
#'  if \code{'inverse'} provides inverse-probability of weighting.   
#' \item \code{numeric},  the user can provides a vector of length \code{n} with weigths for each observation.
#' }
#' @param boot \code{character}, type of boostrap resample,  Values are "global" and "local". See \code{\link{group.bootstrap}} function for mor details.
#' @param par.boot, list with convex boostrap parameter, see \code{\link{group.bootstrap}} function for mor details.
#' @param classif Name of classifier, by default \code{"classif.glm"} then \code{\link{classif.glm}} function is used.
#' Other  classifiers:  \code{\link{classif.gsam}}, \code{\link{classif.rpart}}, \code{\link{classif.nnet}} and \code{\link{classif.svm}}, etc.
#' @param par.classif List of parameters for the classifier.
#' @param verbose \code{logical}, if \code{TRUE} print relevant information for each iteration.
#' @param \dots Further arguments passed to or from other methods.

#' @return Return a list of length sample\code{B}, tha corresond to the \code{B} fitted models. A convex bootstrap sample is generated to fit each 
#'  model, then the bootstrap re is compued. 
#   (by default is \code{nb=100}.  
#'  The model fitted has the same  result of the classification method plus:\cr
#' \itemize{
#'  \item \code{predicted.values}, \code{vector} of size \code{n},  with the predicted response using the fitted model with  \code{nb} bootstrap sample.
#'  \item \code{alpha}, \code{numeric}, error rate value used to ensamble the \code{B} models. 
#'  \item \code{error}, \code{numeric}, the error values of missclassification  (by default is the accuracy).
##' }
#' @noRd 
#' @details 
#' The ensamble (final) classifier is  a linear combination of all of the \code{B} classifiers
#' using the \eqn{\alpha}{\alpha} coefficients.  Breiman coefficient \eqn{\alpha}=1/2 ln((1-err)/err) is computed for each classifier.
#'  
#' @author Febrero-Bande, M. and Oviedo de la Fuente, M.
#' 
#' @seealso See Also as: \code{\link{classif.glm}},  \code{\link{classif.gsam}},  \code{\link{classif.svm}} and other classifiers.\cr Alternative method:
#' \code{\link{classif.adaboost}}.
#' 
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York. 
#' 
#' Breiman, L. (1998): 'Arcing classifiers'. \emph{The Annals of Statistics}, Vol 26, 3, pp. 801-849.
#' 
#' @keywords internal
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn <- phoneme[["learn"]]
#' glearn <- phoneme[["classlearn"]]
#' mtest <- phoneme[["test"]]
#' gtest <- phoneme[["classtest"]]
#' dataf <- data.frame(glearn)
#' # Estimation
#' ij <- c(41:70, 101:200, 241:250)
#' dat = ldata("df" = dataf[ij, , drop = F], "x" = mlearn[ij, ])
#' res1 <- classif.glm(glearn ~ x, dat)
#' # Naive bootstrap
#' res2 <- classif.bootstrap(glearn ~ x, dat, classif = "classif.glm") 
#' # res2 <- classif.bootstrap(glearn ~ x, dat, classif = "classif.glm"
#' ,par.boot=list("smo"=0)) 
#' # Smooth bootstrap
#' res3 <- classif.bootstrap(glearn ~ x, dat, classif = "classif.glm"
#' ,par.boot=list("smo"=0.5)) 
#' # Global convex bootstrap
#' res4 <- classif.bootstrap(glearn ~ x, dat, classif = "classif.glm"
#' ,par.boot=list("Nhull"=4)) 
#' # Local convex bootstrap
#' res5 <- classif.bootstrap(glearn ~ x, dat, classif = "classif.glm"
#' ,par.boot=list("Nhull"=4,"Nnbh"=12)) 
#' # Prediction
#' newdat = list("df" = dataf, "x" = mtest)
#' pred1 <- predict(res1, newdat)
#' pred2 <- predict(res2, newdat)
#' pred3 <- predict(res3, newdat)
#' pred4 <- predict(res4, newdat)
#' pred5 <- predict(res5, newdat)
#' cat2meas(pred1, newdat$df$glearn)
#' cat2meas(pred2, newdat$df$glearn)
#' cat2meas(pred3, newdat$df$glearn)
#' cat2meas(pred4, newdat$df$glearn)
#' cat2meas(pred5, newdat$df$glearn)
#' diag(table(pred1, newdat$df$glearn))/50
#' diag(table(pred2, newdat$df$glearn))/50
#' diag(table(pred3, newdat$df$glearn))/50 
#' diag(table(pred4, newdat$df$glearn))/50
#' diag(table(pred5, newdat$df$glearn))/50  
#' }

# function (x, response, boot = "convex", nb = NULL, smo = 0.05, 
#               Nhull = 4, Nnbh = NULL, weights, metric = NULL) 
  
# @export  classif.bootstrap
classif.bootstrap <- function(formula, data, weights="equal", B = 50
         , par.boot = list( nb = NULL, smo = 0, Nhull  =NULL, Nnbh = NULL) 
         , classif="classif.glm", par.classif, verbose = FALSE, ...) {
  #if (is.null(par.boot$B)) par.boot$B <- 50
  #if (is.null(par.boot$nb)) par.boot$nb <- 1000
  #if (is.null(par.boot$smo)) par.boot$smo <- 0.05
  #if (is.null(par.boot$Nhull)) par.boot$Nhull <- 4
  nam.boot <- names(par.boot)
  if (!("nb" %in% nam.boot )) par.boot[["nb"]] <- NULL
  #if (!("smo" %in% nam.boot)) par.boot[["smo"]] <- 0 # si
  if (!("Nhull" %in% nam.boot)) par.boot[["Nhull"]] <- NULL
  if (!("Nnbh" %in% nam.boot)) par.boot[["Nnbh"]] <- NULL
  if (!("metric" %in% nam.boot)) par.boot[["metric"]] <- NULL
  data0 <- data
  data <- data$df
  response <- as.character(formula[[2]])
  vardep <- data[,response]
  n <- length(vardep)
  nclases <- nlevels(vardep)
  #guardarpesos <- array(0, c(n,B)) #ensamble con B
  #w     <-   pesos <- rep(1/n,n)
  newdata0 <-data0
  #resample[["weights"]]<-weights4class(vardep,type=resample$type)
  if (is.null(weights)) weights="equal"
  if (is.character(weights)) {
    weights<-weights4class(vardep,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
#print("wei");  print(weights[1:8])
  par.boot$response <- response
  par.boot[["x"]]<-data0
  data.ori<-data0
  if (missing(par.classif))
    par.classif <- list()
  if (is.null(par.classif$weights)) {
    par.classif$weights="equal"
  }
  par.classif <-c(list(formula=formula),par.classif)
  list.fit <- list() 
  # if (is.null(par.boot$weights)) {    par.boot$weights=weights  }
  alpha.vec<-numeric(B)
  for (i in 1:B){
     #cat(i)
    #print(names(par.boot))
    #print(par.boot[])
    par.classif$data <- do.call("group.bootstrap",par.boot)
    #if (boot=="global"){
      # print("convex global")
      #ik <- which(names(par.boot) =="Nnbh")
      #par.classif$data <- do.call("group.bootstrap",par.boot)
      # function (x, response, boot = "convex", nb = NULL, smo = 0.05, 
      #               Nhull = 4, Nnbh = NULL, weights, metric = NULL)       
      #  par.classif$data <- group.bootstrap(x=par.boot$x, response=par.boot$response
      #                   , nb=par.boot$nb, smo=par.boot$smo, Nhull=par.boot$Nhull
      #                , Nhull=par.boot$Nnbh,weights=weights)
      # print(dim(data0[[1]]))
       #(data,response,J=4,expand=1,B=5000,weigths)
    # } else{
    #print("local block bootstrap")
    # if (is.null(par.boot$metric)) metric <- NULL
      # par.classif$data <- group.bootstrap(x=par.boot$x
      # , response=par.boot$response, weights=weights,boot=boot
      # , nb=par.boot$nb, smo=par.boot$smo, Nhull=par.boot$Nhull
      # , Nnbh=par.boot$Nnbh, metric = metric)
     # par.classif$data <- do.call("group.bootstrap",par.boot)
      # }
  
 if (!is.null(par.classif$basis.x)){
   #print("creando base pc1")
   basis.aux <- par.classif$basis.x
   nam.x <- names(par.classif$basis.x)
   for (j in 1:length(par.classif$basis.x)){
     if (is(par.classif$basis.x[[nam.x[j]]],"fdata.comp")){
       if (par.classif$basis.x[[nam.x[j]]]$type=="pc"){
           basis.aux[[nam.x[j]]] <- create.pc.basis(par.classif$data[[nam.x[j]]]
                              ,par.classif$basis.x[[nam.x[j]]]$l)
       }}     #falta incluir para PLS basiss
   }
   basis.aux -> par.classif$basis.x
   }

    list.fit[[i]] <- do.call(classif,par.classif) 
    list.fit[[i]]$predicted.values <- predict(list.fit[[i]],data0)
#    cat(i," ",list.fit[[i]]$max.prob)
    #print(length(vardep));    print(length(list.fit[[i]]$predicted.values))
    
    aux<-cat2alpha(yobs = vardep, 
                   ypred = list.fit[[i]]$predicted.values,
                   weights = weights)
    list.fit[[i]]$alpha = aux$alpha
    list.fit[[i]]$error = aux$error
    alpha.vec[i]<-aux$alpha
#    list.fit[[i]]$predicted.prob= cat2accuracy(vardep,list.fit[[i]]$predicted.values)
    #print(list.fit[[i]]$max.prob)
      #weigthed.mean(list.fit[[i]]$predicted.values != vardep, par.S$w)
  }#fin B
  ny <- levels(vardep)
  # out2glm <- classifKgroups(vardep, prob.group, ny)
  # yest <- out2glm$yest
  #prob1 <- out2glm$prob1
  #prob.group <- out2glm$prob.group

  C <- match.call()
  output<- list("group"=vardep)
  #output$prob.classification <- 
  output$list.fit <- list.fit
  output$alpha <- alpha.vec
  output$call <- C[1:2]
  output$C <- C
  output$levels <-ny
  output$prob <- .5
  output$B <- B
  class(output) <- "classif"
  
  out2 <- predict.classif.bootstrap(output, data0, type = "probs") 
  output$group.est <- out2[[1]]
  output$prob.group <- out2[[2]]
  tab <- table( output$group.est,  output$group)
  output$prob.classification =  diag(tab)/colSums(tab)
  output$max.prob <-  mean(output$group.est ==  vardep)
  output$misclassification <- 1-output$max.prob
  return(output)
}

# subset.fdata <- fda.usc.devel:::subset.fdata
# subset.ldata <- fda.usc.devel:::subset.ldata
# cat2alpha <- fda.usc.devel:::cat2alpha
# predict.classif.bootstrap <- fda.usc.devel:::predict.classif.bootstrap
# mfdata.cen <- fda.usc.devel:::mfdata.cen

# @export  classif.bootstrap2
# @keywords internal
# @noRd
classif.bootstrap2 <- function(formula, data, weights="equal", B = 50
                               , par.boot = list( nb = NULL, smo = 0, Nhull  =NULL, Nnbh = NULL) 
                               , classif="classif.glm", par.classif, verbose = FALSE, ...) {
  #if (is.null(par.boot$B)) par.boot$B <- 50
  #if (is.null(par.boot$nb)) par.boot$nb <- 1000
  #if (is.null(par.boot$smo)) par.boot$smo <- 0.05
  #if (is.null(par.boot$Nhull)) par.boot$Nhull <- 4
  nam.boot <- names(par.boot)
  if (!("nb" %in% nam.boot )) par.boot[["nb"]] <- NULL
  #if (!("smo" %in% nam.boot)) par.boot[["smo"]] <- 0 # si
  if (!("Nhull" %in% nam.boot)) par.boot[["Nhull"]] <- NULL
  if (!("Nnbh" %in% nam.boot)) par.boot[["Nnbh"]] <- NULL
  if (!("metric" %in% nam.boot)) par.boot[["metric"]] <- NULL
  data0 <- data
  data <- data$df
  response <- as.character(formula[[2]])
  vardep <- data[,response]
  n <- length(vardep)
  nclases <- nlevels(vardep)
  #guardarpesos <- array(0, c(n,B)) #ensamble con B
  #w     <-   pesos <- rep(1/n,n)
  newdata0 <-data0
  #resample[["weights"]]<-weights4class(vardep,type=resample$type)
  if (is.null(weights)) weights="equal"
  if (is.character(weights)) {
    weights <- weights4class(vardep, type=weights)
    
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  #print("wei");  print(weights[1:8])
  par.boot$response <- response
  par.boot[["x"]]<-data0
  data.ori<-data0
  if (missing(par.classif))
    par.classif <- list()
  if (is.null(par.classif$weights)) {
    par.classif$weights="equal"
  }
  par.classif <-c(list(formula=formula),par.classif)
  list.fit <- list() 
  # if (is.null(par.boot$weights)) {    par.boot$weights=weights  }
  alpha.vec <- numeric(B)
  
  par.fda.usc <- eval(parse(text="fda.usc.devel:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  #group.bootstrap<-fda.usc.devel:::group.bootstrap
  #ncores <- ops.fda.usc()$ncores
  inumgr <- icount(B)
  list.fit <-list()
  list.fit <-foreach(i = inumgr, .packages='fda.usc.devel') %dopar% {
    #library(fda.usc)
    #cat2alpha <- fda.usc.devel:::cat2alpha
    par.classif$data <- do.call("group.bootstrap",par.boot)
    # for (i in 1:B){
    #par.classif$data <- laux[[i]]
    if (!is.null(par.classif$basis.x)){
      #print("creando base pc1")
      basis.aux <- par.classif$basis.x
      nam.x <- names(par.classif$basis.x)
      for (j in 1:length(par.classif$basis.x)){
        if (is(par.classif$basis.x[[nam.x[j]]],"fdata.comp")){
          if (par.classif$basis.x[[nam.x[j]]]$type=="pc"){
            basis.aux[[nam.x[j]]] <- create.pc.basis(par.classif$data[[nam.x[j]]]
                                                     ,par.classif$basis.x[[nam.x[j]]]$l)
          }}     #falta incluir para PLS basiss
      }
      basis.aux -> par.classif$basis.x
    }
    lfit <- do.call(classif,par.classif) 
    lfit$predicted.values <- predict(lfit,data0)
#    lfit
#  }
#  list.fit <- laux
#  for (i in 1:B){
    aux<-cat2alpha(yobs = vardep, ypred = lfit$predicted.values, weights = weights)
    lfit$alpha = aux$alpha
    lfit$error = aux$error
    lfit
  }
  for (i in 1:B){    
    alpha.vec[i] <- list.fit[[i]]$alpha
  }#fin B
  ny <- levels(vardep)
  # out2glm <- classifKgroups(vardep, prob.group, ny)
  # yest <- out2glm$yest
  #prob1 <- out2glm$prob1
  #prob.group <- out2glm$prob.group
  
  C <- match.call()
  output<- list("group"=vardep)
  #output$prob.classification <- 
  output$list.fit <- list.fit
  output$alpha <- alpha.vec
  output$call <- C[1:2]
  output$C <- C
  output$levels <-ny
  output$prob <- .5
  output$B <- B
  class(output) <- "classif"
  out2 <- predict.classif.bootstrap(output, data0, type = "probs") 
  output$group.est  <- out2[[1]]
  output$prob.group <- out2[[2]]
  tab <- table( output$group.est, output$group)
  output$prob.classification =  diag(tab)/colSums(tab)
  output$max.prob <- mean(output$group.est == vardep)
  output$misclassification <- 1-output$max.prob
  return(output)
}

#fda.usc.devel:::classif.bootstrap2<-classif.bootstrap2
# subset.fdata <- fda.usc.devel:::subset.fdata
# subset.ldata <- fda.usc.devel:::subset.ldata
# cat2alpha <- fda.usc.devel:::cat2alpha
# predict.classif.bootstrap <- fda.usc.devel:::predict.classif.bootstrap
# mfdata.cen <- fda.usc.devel:::mfdata.cen
#  group.bootstrap <- fda.usc.devel:::group.bootstrap

# pboot <- list("nb"=100,"Nnbh"=4)
# 
# ops.fda.usc()
# set.seed(1)
# t1<-proc.time()
# res2 <- classif.bootstrap(glearn ~ x, dat, B=50, classif = "classif.glm",par.boot=pboot)
# t2<-proc.time()
# 
# set.seed(1)
# t3<-proc.time()
# res3 <- classif.bootstrap2(glearn ~ x, dat, B=50, classif = "classif.glm",par.boot=pboot)
# t4<-proc.time()
# traceback()
# t2-t1
# t4-t3
# 
# summary(res2)
# summary(res3)
