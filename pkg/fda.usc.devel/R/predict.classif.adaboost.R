
# @title Predicts from a fitted object using bootstrap or adaboost classifier.
# 
# @description Classifier of multivariate and functional data by bootstrap or adaboost classifier.
# Returns the predicted classes using a previously trained model.
# 
# @param object \code{object} estimated by \code{classif} method.
# @param newdata New functional explanatory data of \code{ldata} class.
# @param type Type of prediction ("class or probability of each group
# membership").
# @param \dots Further arguments passed to or from other methods.
# 
# @return If type="class", produces a vector of predictions.
# If type="probs", a list with the following components is returned: 
# \itemize{
# \item \code{group.pred} the vector of predictions.
# \item \code{prob.group} the matrix of predicted probability by factor level. 
# }
# 
# @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
# 
# @seealso See Also as: \code{\link{classif.adaboost}} and  \code{\link{classif.bootstrap}}
# 
# @keywords classif
# 
# @rdname predict.classif.adaboost
# @export 
predict.classif.bootstrap<- function (object, newdata, 
                                      type = "class",...) 
{#, q.best=0
 #print("pred  classif bootstrap")
  B <- object$B
#print(names(object))  
  group<-object$group
  pond <- object$alpha
  object <- object$list.fit 
  #group <- object[[1]]$group
  n <- length(group)
  nlev <- nlevels(group)
  lev <- levels(group)
  numgr2 <- NROW(newdata[[1]])
  classfinal <- array(0, c(numgr2, nlev))
  pred <- data.frame(rep(0, numgr2))
  pond <- rep(1,len=B)# poner el peso de cada classificador
  #print("aaa1")  
  # revisar codigo probs
  # if (type == "probs") {
  #   for (m in 1:B) {
  #     p <- predict.classif(object[[m]], newdata, type = "probs")$prob.group[, 
  #                                                                           2]
  #     if (m > 1) 
  #       pred <- data.frame(pred, p)
  #     else pred <- p
  #   }
  #   classfinal <- pred# as.numeric(as.matrix(pred) %*% matrix(pond))
  #   predclass <- factor(ifelse(classfinal < 0.5, lev[1], 
  #                              lev[2]), levels = lev)
  # }
  # else {
  for (m in 1:B) {
  #  cat(m)
    #  if (object[[1]]$C[[1]] == "classif.ML") 
    # p <- predict.classif.ML(object[[m]], newdata, type = "class")
    #else
    p <- predict.classif(object[[m]], newdata, type = "class")
    if (m > 1) 
      pred <- data.frame(pred, p)
    else pred <- p
    #pond[m]<-object[[m]]$alpha
  }
  names(pred)<-1:B
  #qu<-quantile(pond,probs=q.best)
  #pond[pond<qu]<-0 
  pond<-pond/sum(pond)
#print("aaa2")
  for (i in 1:nlev) {
    classfinal[, i] <- matrix(as.numeric(pred == lev[i]), 
                              nrow = numgr2) %*% as.vector(pond)
  }
  #print(classfinal)
  #print(pond)
  predclass <- lev[apply(classfinal, 1, which.max)]
  predclass <- factor(predclass, levels = lev)
  #print("aaa3")  
  if (type=="class")
    return(predclass)  else
      return(list("group.pred"=predclass,"prob.group"=classfinal,
                  "votes"=pred))
}

pred2boot<- function (object, newdata, 
                                      type = "class",...) 
{#, q.best=0
  #print("pred  classif bootstrap")
  B <- object$B
  #print(names(object))  
  group<-object$group
  pond <- object$alpha
  object <- object$list.fit 
  #group <- object[[1]]$group
  n <- length(group)
  nlev <- nlevels(group)
  lev <- levels(group)
  numgr2 <- NROW(newdata[[1]])
  classfinal <- array(0, c(numgr2, nlev))
  pred <- data.frame(rep(0, numgr2))
  pond <- rep(1,len=B)# poner el peso de cada classificador
  #print("aaa1")  
  # revisar codigo probs
  # if (type == "probs") {
  #   for (m in 1:B) {
  #     p <- predict.classif(object[[m]], newdata, type = "probs")$prob.group[, 
  #                                                                           2]
  #     if (m > 1) 
  #       pred <- data.frame(pred, p)
  #     else pred <- p
  #   }
  #   classfinal <- pred# as.numeric(as.matrix(pred) %*% matrix(pond))
  #   predclass <- factor(ifelse(classfinal < 0.5, lev[1], 
  #                              lev[2]), levels = lev)
  # }
  # else {
  for (m in 1:B) {
    #  cat(m)
    #  if (object[[1]]$C[[1]] == "classif.ML") 
    # p <- predict.classif.ML(object[[m]], newdata, type = "class")
    #else
    p <- predict.classif(object[[m]], newdata, type = "class")
    if (m > 1) 
      pred <- data.frame(pred, p)
    else pred <- p
    #pond[m]<-object[[m]]$alpha
  }
  names(pred)<-1:B
  #qu<-quantile(pond,probs=q.best)
  #pond[pond<qu]<-0 
  pond<-pond/sum(pond)
  #print("aaa2")
  for (i in 1:nlev) {
    classfinal[, i] <- matrix(as.numeric(pred == lev[i]), 
                              nrow = numgr2) %*% as.vector(pond)
  }
  #print(classfinal)
  #print(pond)
  predclass <- lev[apply(classfinal, 1, which.max)]
  predclass <- factor(predclass, levels = lev)
  #print("aaa3")  
  if (type=="class")
    return(predclass)  else
      return(list("group.pred"=predclass,"prob.group"=classfinal,
                  "votes"=pred))
}
# @rdname predict.classif.adaboost
# @export 
predict.classif.adaboost=function(object, newdata, type = "class",...){
  
  # object<-out4a
  # newdata<-ltrain
  # newdata<-ltest
  #  type = "prob"
  
#   nldata<-length(newdata)
#   ny<-levels(object$group)
    fit<-object$list.fit
    mfinal<- length(object$list.fit)
    group<-object$group
    n <- length(group)
    nlev <- nlevels(group)
    lev <- levels(group)

    #numg=vardep=object$group
    
    pond <- object$alpha.boost
  
    numgr2<-nrow(newdata[[1]])
    classfinal <- array(0, c(numgr2, nlev))
    
#    newdata <- data.frame(newdata, pesos)
    pred <- data.frame(rep(0, numgr2))
    #ind<-1
    if (type=="probs") {
      for (m in 1:mfinal) {#binario
       # if (object$call[[1]]=="classif.ML")
       #   p<-predict.classif.ML(fit[[m]],newdata,type="probs")$prob.group[,2]
        #else                              
          #p<-predict.classif(fit[[m]],newdata,type="probs")$prob.group[,2]
          p<-predict(fit[[m]],newdata,type="probs")$prob.group[,2]
        if (m> 1) pred<-data.frame(pred,p)  else    pred <-p
      }
      pond <- pond/sum(pond)
      probclass <- as.numeric(as.matrix(pred)%*%matrix(pond))
      predclass <- factor(ifelse(probclass <.5,lev [1],lev[2]),levels=lev)
    }
    #else{ 
     
      for (m in 1:mfinal) {
        #p<-predict.classif(fit[[m]],newdata,type="probs")$group.pred #idem de:
        #if (object$call[[1]]=="classif.ML"){
          #p<-predict.classif.ML(fit[[m]],newdata,type="class")
         # p<-predict(fit[[m]],newdata,type="class")
        #}
        #else
          p<-predict(fit[[m]],newdata,type="class")
          #p<-predict.classif(fit[[m]],newdata,type="class")
        if (m> 1) pred<-data.frame(pred,p)  else    pred <-p
      }
      for (i in 1:nlev){
        classfinal[,i] <- matrix(as.numeric(pred==lev[i]),nrow=numgr2)%*%as.vector(pond)
      }
      #2019-04-10   
      #predclass<-apply(classfinal,1,FUN=select, vardep.summary=summary(vardep))
      predclass <-lev[apply(classfinal,1,which.max)] 
      predclass <- factor(predclass,levels=lev)
      if (type=="probs") 
        predclass<-list("group.pred"=predclass, "prob.group"=probclass)
      
  #}
return(predclass)
}    

# predict.classif.adaboost=function(object, newdata, type = "class",...){
#   
#   # object<-out4a
#   # newdata<-ltrain
#   # newdata<-ltest
#   #  type = "prob"
#   
#   #   nldata<-length(newdata)
#   #   ny<-levels(object$group)
#   fit<-object$list.fit
#   mfinal<- length(object$list.fit)
#   group<-object$group
#   n <- length(group)
#   nlev <- nlevels(group)
#   lev <- levels(group)
#   
#   #numg=vardep=object$group
#   
#   pond <- object$alpha.boost
#   
#   numgr2<-nrow(newdata[[1]])
#   classfinal <- array(0, c(numgr2, nlev))
#   
#   #    newdata <- data.frame(newdata, pesos)
#   pred <- data.frame(rep(0, numgr2))
#   #ind<-1
#   
#   
#   if (type=="probs") {
#     for (m in 1:mfinal) {#binario
#       # if (object$call[[1]]=="classif.ML")
#       #   p<-predict.classif.ML(fit[[m]],newdata,type="probs")$prob.group[,2]
#       #else                              
#       p<-predict.classif(fit[[m]],newdata,type="probs")$prob.group[,2]
#       if (m> 1) pred<-data.frame(pred,p)  else    pred <-p
#     }
#     classfinal <- as.numeric(as.matrix(pred)%*%matrix(pond))
#     predclass <- factor(ifelse(classfinal <.5,lev [1],lev[2]),levels=lev)
#   }
#   else{ 
#     
#     for (m in 1:mfinal) {
#       #p<-predict.classif(fit[[m]],newdata,type="probs")$group.pred #idem de:
#       if (object$call[[1]]=="classif.ML")
#         p<-predict.classif.ML(fit[[m]],newdata,type="class")
#       else
#         p<-predict.classif(fit[[m]],newdata,type="class")
#       if (m> 1) pred<-data.frame(pred,p)  else    pred <-p
#     }
#     for (i in 1:nlev){
#       classfinal[,i] <- matrix(as.numeric(pred==lev[i]),nrow=numgr2)%*%as.vector(pond)
#     }
#     #2019-04-10   
#     #predclass<-apply(classfinal,1,FUN=select, vardep.summary=summary(vardep))
#     predclass <-lev[apply(classfinal,1,which.max)] 
#     predclass <- factor(predclass,levels=lev)
#     predclass
#   }
#   return(predclass)
#   # ojo  solo valido para binomial       
#   
# }   