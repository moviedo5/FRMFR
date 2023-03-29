#################################################################
# Warping of adabag::boosting function
#################################################################

#' @title Functional Classification using AdaBoost Algorithm
#' 
#' @description Computes classification using functional (and non functional)
#' explanatory variables by AdaBoost algorithm.  
#' The ensamble classifier is  a linear combination of all of the (weak) classifiers
#'  using the \eqn{\alpha}{\alpha} weight updating coefficients.
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param data \code{list} that containing the variables in the model. 
#' The first element in the \code{data} list is called \emph{"df"} and it is a data
#' frame with the response and non functional explanatory variables, as \code{\link{glm}}. 
#' Functional covariates of class \code{fdata}  are introduced in the following items in the \code{data} list.
#' @param classif Name of classifier, by default \code{"classif.glm"} then \code{\link{classif.glm}} function is used.
#' Other  classifiers:  \code{\link{classif.gsam}}, \code{\link{classif.rpart}}, \code{\link{classif.nnet}} and \code{\link{classif.svm}}, etc.
#' @param par.classif List of parameters for the classifier.
#' @param B \code{integer}, maximum number of iterations.
#' @param alpha.boost \code{character}.   Three types: 
#' \itemize{
#'  \item Breiman (by default), \eqn{\alpha}=1/2 ln((1-err)/err)
#'  \item Freund, \eqn{\alpha}=ln((1-err)/err)
#'  \item Zhu, \eqn{\alpha}=ln((1-err)/err) + ln(nclasses-1)
#' }
#' where \code{err} is the misclassification error.
#' @param weights Initial weights:   
#' \itemize{
#'  \item \code{character}, if \code{'equal'} (by default)  provides same weights for each observation; 
#'  if \code{'inverse'} provides inverse-probability of weighting.   
#' \item \code{numeric},  the user can provides a vector of length \code{n} with weigths for each observation.
#' }
#' @param verbose \code{logical}, if \code{TRUE} prints relevant information for each iteration.
#' @param \dots Further arguments passed to or from other methods.

#' @aliases classif.adaboost 
#' @return Return the same result of the classification method plus:
#' \itemize{
#' \item \code{group.est.boost}, \code{matrix} of \code{n X B},  the fitted response (by row) in each iteration (by column).
#' \item \code{alpha.boost}, \code{vector} of length \code{B} with the  \eqn{\alpha} coefficients.  
#' The predicted classes are a linear combination of all of the classifiers weighted by \code{alpha.boost}.
#' \item \code{votes.boost}, \code{matrix} of \code{n X g},  with   votes values (by row)  for each class group  (by column).
#' \item \code{prob.boost}, \code{matrix} of \code{n X g},   boosting probabilities of each observation (by row) for each class group  (by column).
#' \item \code{weigths.rep}, \code{matrix} of \code{n X B},  with  the weigths values (by row) in each iteration (by column).
#' \item \code{list.fit}, \code{list} with the \code{B} fitted models. 
#' }
#' 
#' @author Febrero-Bande, M. and Oviedo de la Fuente, M.
#' @seealso See Also as: \code{\link{classif.glm}},  \code{\link{classif.gsam}},  \code{\link{classif.svm}} and other classifiers.\cr Alternative method:
#' \code{\link{classif.bootstrap}}.

#' @references 
#' Breiman, L. (1998): 'Arcing classifiers'. \emph{The Annals of Statistics}, Vol 26, 3, pp. 801-849.\cr
#' Freund, Y. and Schapire, R.E. (1996): "Experiments with a new boosting algorithm". In Proceedings of the Thirteenth International Conference on Machine Learning, pp. 148-156, Morgan
#' Kaufmann.
#' 
#' Zhu, J., Zou, H., Rosset, S. and Hastie, T. (2009): 'Multi-class AdaBoost.  \emph{Statistics and Its
#' Interface}, 2, pp. 349-360.
#' 
#' Alfaro, E., Gamez, M. and Garcia, N. (2013): 'adabag: An R Package for Classification with
#' Boosting and Bagging'.  \emph{Journal of Statistical Software}, Vol 54, 2, pp. 1-35.

#' @keywords classif
#' @note Adapted version from \code{adabag::boosting} function.
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
#' dat = list("df" = dataf[ij, , drop = F], "x" = mlearn[ij, ])
#' res1 <- classif.glm(glearn ~ x, data = dat)
#' res2 <- classif.adaboost(glearn ~ x,data = dat,
#'                    classif = "classif.glm", B = 25)
#' res1$max.prob;res2$max.prob

#' # Prediction
#' newdat = list("df" = dataf, "x" = mtest)
#' pred1 <- predict(res1, newdat)
#' pred2 <- predict(res2, newdat)
#' cat2meas(pred1, newdat$df$glearn); cat2meas(pred2, newdat$df$glearn)
#' diag(table(pred1, newdat$df$glearn))/50
#' diag(table(pred2, newdat$df$glearn))/50
#' }
#' 
#' @export classif.adaboost
classif.adaboost<-function(formula, data,
         classif="classif.glm", par.classif, 
         B=100, alpha.boost="Breiman", weights="equal",
         verbose = FALSE, ...) {
#  print("entra classifidfkadfk adaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
#print(B)
  # formula <- glearn ~ x
  # data <- ldata
  # classif="classif.rpart"
  # classif="classif.glm"
  # par.classif<-list()
  # B=100
  # alpha.boost="Breiman"
  # boos=FALSE
  C <- match.call()
  #Exigimos que alpha.boost sea uno de esos dos valores
  if (!(as.character(alpha.boost) %in% c("Freund","Breiman","Zhu"))){
    stop("alpha.boost must be 'Freund', 'Breiman' or 'Zhu' ")
  }
  formula<- as.formula(formula)
  
  data0 <- data
  data <- data$df
  response <- as.character(formula[[2]])
  vardep <- data[,response]
  n <- length(vardep)
  nclases <- nlevels(vardep)
  if (is.character(weights)) {
    weights <- weights4class(vardep ,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  pesos<-w<-weights
  guardarpesos <- array(0, c(n,B)) #para ver los pesos de las observaciones
  newdata0 <-data0
  #data<-cbind(pesos, data) #Los pesos en rpart deben ser una columna del dataframe
  list.fit <- list() #Creamos una lista para guardar los modelos estimados
  pond <- rep(0,B) # Un vector donde guardaremos la ponderacion de cada arbol.
  pred<- data.frame(rep(0,n))

  if (missing(par.classif))
    par.classif <- list()
  par.classif <-c(list(formula=formula,data=data0),par.classif)
  
  maxerror<-0.5
  eac <- 0.0001 # minimum fraction of error above 0 or under maxerror 
  c.max<- log((1-eac)/eac) # se puede cambiar por 10
  
  if (alpha.boost=="Zhu"){		maxerror<-1-1/nclases		}
  # repeticiones
  pesos<-pesos/sum(pesos)
  for (m in 1:B) {
# print(m) 
  #Creamos muestras boostrap utilizando los pesos
#  	w<<- pesos
  	par.classif$weights <- pesos
   	fit <- do.call(classif,par.classif)
  	list.fit[[m]] <- fit	#Guardamos los B modelos estimados pero se podr?a parar sino se mejora o se empeora
	  flearn <- fit$group.est
    ind<-as.numeric(vardep != flearn) #Crear un vector indicador 
    err<- mean(pesos*ind)         #Calcula el error ponderado en esa iteracion
    c<- log((1-err)/err)
    if (alpha.boost=="Breiman"){	c <- .5*c	}
	  if (alpha.boost=="Zhu"){	c <- .5*c+log(nclases-1)	}
   # c[c<=0] <- 0
    c[c>=c.max] <- c.max
	  guardarpesos[,m] <- pesos
    pesos <- pesos*exp(c*ind)
    #pesos <- pesos/sum(pesos)
   
    # Si el error no es menor que la regla por defecto los pesos se deberian inicializrr, o no
    # si no 0<err<0.5 se ajustan los valores de a 3 y 0.001

    if (verbose) cat("Iter: ",m,"error",err,"maxerror",maxerror,"c",c,"sum pesos",sum(pesos),"\n")
  if (c==0 ) break
    if (err>=maxerror) {
      if (verbose) print("error>maxerror")
      pesos <- weights
      maxerror<-maxerror-eac
      c <- log((1-maxerror)/maxerror)
      if (alpha.boost=="Breiman"){      	c <- (1/2)*c      	}
      if (alpha.boost=="Zhu"){      	c <- c + log(nclases-1) 	}
    } 
    pond[m]<- c 		
    if (m==1) {
      pred <- flearn
    }    else { pred <- data.frame(pred,flearn) }
    if (err==0) {
#      if (verbose) print("error 0,    pesos 1/n")
      pesos <- weights
      c <- c.max
      break
    }
    
    
    } # fin repeticiones
#print("sale1");print(dim(pred));print(pond)
  classfinal <- array(0, c(n,nlevels(vardep)))
  lev <-levels(vardep)
  prob2<-prob1 <- nlev <- nlevels(vardep)#ngroup <-
  
  if (c<=0 )     {
    pond<- pond[1:(m-1)]
    list.fit<-list.fit[1:(m-1)]
    guardarpesos<-guardarpesos[,1:(m-1)]
  }
  if (c>=c.max)     {    pond<- pond[1:m];    list.fit<-list.fit[1:m];    guardarpesos<-guardarpesos[,1:m]  }
  
  for (i in 1:nlev){
    classfinal[,i] <- matrix(as.numeric(pred==lev[i]),nrow=n)%*%as.vector(pond)
  }
  # print("sal bucle2")   
  #2019-04-10   
  #predclass<-apply(classfinal,1,FUN=select, vardep.summary=summary(vardep))
  predclass <-lev[apply(classfinal,1,which.max)] 
  predclass <- factor(predclass,levels=lev)
  #print(length(predclass));  print(dim(classfinal))
  #Para que devuelva las probabilidades a posteriori
  votosporc <- classfinal/apply(classfinal,1,sum)
  max.prob=mean(predclass==vardep)
  tab <- table(predclass, vardep)
  prob.group <- votosporc
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:nlev) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- lev
  colnames(prob.group) <- lev
  output<-list(formula=formula,data=data,group=vardep,group.est =predclass,
               prob.classification=prob1,prob.group=prob.group, 
               max.prob=max.prob,C=C,
               group.est.boost=pred,alpha.boost=pond, votes.boost=classfinal,
              prob.boost=votosporc, weigths.rep=guardarpesos, list.fit = list.fit)
  
  tf <- terms.formula(formula)
  output$terms <- attr(tf, "term.labels")
  output$call <-list.fit[[1]]$C[1:2]
  #if (resample$type=="bootstrap") hacer prediccion?
  class(output)<-"classif"  
  return( output)
}




