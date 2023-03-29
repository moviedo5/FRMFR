##############
gower.dist <- function (data.x, data.y = data.x, rngs = NULL, KR.corr = TRUE, 
          var.weights = NULL) 
{
  library(StatMatch)
  
  x1 <- as.logical(rbinom(10,1,0.5)) 
  x2 <- sample(letters, 10, replace=TRUE)
  x3 <- rnorm(10)
  x4 <- ordered(cut(x3, -4:4, include.lowest=TRUE))
  xx <- data.frame(x1, x2, x3, x4, stringsAsFactors = FALSE)
  
  # matrix of distances between observations in xx
  dx <- gower.dist(xx)
  head(dx)
  #  metric.dist(data.x[,1,drop=F],data.y[,1,drop=F])
  # funciona para variables logicas 0 o 1
  # si misma letra 0, si diferente 1
  
  # matrix of distances between first obs. in xx
  # and the remaining ones
  gower.dist(data.x=xx[1:6,], data.y=xx[7:10,], var.weights = c(1,2,5,2))
  d4<-gower.fcn(, xx[7:10,4], rng = NULL, KR.corr = TRUE)
  
# adapatar para que funciones con metric.lp y metric.dist si
# los factores estan ordenados

  gower.fcn <- function(x, y, rng = NULL, KR.corr = TRUE) {
    nx <- length(x)
    ny <- length(y)
    cx <- class(x)
    cy <- class(y)
    delta <- matrix(1, nx, ny)
    if (!identical(cx, cy)) 
      stop("the x and y object are of different type")
    if (is.logical(x)) {
      dd <- abs(outer(X = x, Y = y, FUN = "-"))
      delta[outer(x == FALSE, y == FALSE, FUN = "&")] <- 0
      delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
    }
    else if (is.character(x) || (is.factor(x) && !is.ordered(x))) {
      if (is.factor(x) && !identical(levels(x), levels(y))) 
        stop("x and y have different levels")
      dd <- 1 - outer(x, y, FUN = "==")
      delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
    }
    else if (is.ordered(x)) {
      if (KR.corr) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        if (is.null(rng) || is.na(rng)) 
          rng <- max(x, y, na.rm = TRUE) - 1
        if (rng == 0) {
          dd <- matrix(0, nx, ny)
          delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
        }
        else {
          zx <- (x - 1)/rng
          zy <- (y - 1)/rng
          dd <- abs(outer(X = zx, Y = zy, FUN = "-"))/(max(zx, 
                                                           zy, na.rm = TRUE) - min(zx, zy, na.rm = TRUE))
          delta[outer(is.na(zx), is.na(zy), FUN = "|")] <- 0
        }
      }
      else {
        x <- as.numeric(x)
        y <- as.numeric(y)
        if (is.null(rng) || is.na(rng)) 
          rng <- max(x, y, na.rm = TRUE) - 1
        if (rng == 0) 
          dd <- matrix(0, nx, ny)
        else dd <- abs(outer(X = x, Y = y, FUN = "-"))/rng
        delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
      }
    }
    else {
      if (is.null(rng) || is.na(rng)) 
        rng <- max(x, y, na.rm = TRUE) - min(x, y, na.rm = TRUE)
      if (rng == 0) 
        dd <- matrix(0, nx, ny)
      else dd <- abs(outer(X = x, Y = y, FUN = "-"))/rng
      delta[outer(is.na(x), is.na(y), FUN = "|")] <- 0
    }
    list(dist = dd, delta = delta)
  }

  
    if (is.null(dim(data.x)) && is.null(dim(data.y))) {
    out.gow <- gower.fcn(x = data.x, y = data.y, rng = rngs, 
                         KR.corr = KR.corr)
    out <- (out.gow$dist * out.gow$delta)/out.gow$delta
  }
  else if (is.null(dim(data.x)) && !is.null(dim(data.y))) {
    p <- ncol(data.y)
    if (length(data.x) != p) 
      stop("data.x should be of the same length of the no. of cols of data.y")
    num <- array(0, c(1, nrow(data.y)))
    den <- array(0, c(1, nrow(data.y)))
    if (is.null(var.weights)) 
      var.weights <- rep(1, p)
    for (k in 1:p) {
      if (is.null(rngs)) 
        rng.k <- NULL
      else rng.k <- rngs[k]
      w.k <- var.weights[k]
      out.gow <- gower.fcn(x = data.x[, k], y = data.y[, 
                                                       k], rng = rng.k, KR.corr = KR.corr)
      n <- out.gow$dist * out.gow$delta * w.k
      n[is.na(n)] <- 0
      num <- num + n
      d <- out.gow$delta * w.k
      d[is.na(d)] <- 0
      den <- den + d
    }
    out <- num/den
  }
  else {
    p <- ncol(data.y)
    
    if (ncol(data.x) != p) 
      stop("data.x and data.y must have the same no. of cols")
    num <- array(0, c(nrow(data.x), nrow(data.y)))
    den <- array(0, c(nrow(data.x), nrow(data.y)))
    if (is.null(var.weights)) 
      var.weights <- rep(1, p)
    for (k in 1:p) {
      if (is.null(rngs)) 
        rng.k <- NULL
      else rng.k <- rngs[k]
      w.k <- var.weights[k]
      out.gow <- gower.fcn(x = data.x[, k], 
                  y = data.y[, k], rng = rng.k, KR.corr = KR.corr)
      n <- out.gow$dist * out.gow$delta * w.k
      n[is.na(n)] <- 0
      num <- num + n
      d <- out.gow$delta * w.k
      d[is.na(d)] <- 0
      den <- den + d
    }
    out <- num/den
  }
  out
}

##############

# library(StatMatch) 
fact2dummy <- function (data, all = TRUE, lab = "x") 
{
  dum.fcn <- function(x, all = TRUE) {
    fine <- model.matrix(~x - 1)
    colnames(fine) <- levels(x)
    if (!all) 
      fine <- fine[, -ncol(fine), drop = FALSE]
    fine
  }
  if (is.null(dim(data))) {
    if (class(data)[1] == "numeric" || class(data)[1] == 
        "integer" || class(data)[1] == "logical") 
      oo <- cbind(1 * data)
    else {
      oo <- dum.fcn(data, all = all)
      colnames(oo) <- paste(lab, colnames(oo), sep = "")
    }
  }
  else {
    if (is.matrix(data)) 
      oo <- data
    else {
      p <- ncol(data)
      n <- nrow(data)
      vtype <- lapply(data, class)
      out <- as.list(rep(NA, p))
      lab <- names(data)
      for (i in 1:p) {
        if (vtype[[i]][1] == "numeric" || vtype[[i]][1] == 
            "integer" || vtype[[i]][1] == "logical") {
          aa <- matrix(1 * data[, i])
          colnames(aa) <- lab[i]
          out[[i]] <- aa
        }
        else {
          aa <- dum.fcn(data[, i], all = all)
          colnames(aa) <- paste(lab[i], 1:ncol(aa), sep = "")
          out[[i]] <- aa
        }
      }
      oo <- do.call("cbind", out)
    }
  }
  rownames(oo) <- rownames(data)
  oo
}
##############  
fac <- factor(sample(1:,20,replace=T))
fac
fac2<-ordered(fac)

metric.dist(fac)
metric.dist(fac2)


x1 <- as.logical(rbinom(10,1,0.5)) 
x2 <- sample(letters, 10, replace=TRUE)
x3 <- rnorm(10)
x4 <- ordered(cut(x3, -4:4, include.lowest=TRUE))
xx <- data.frame(x1, x2, x3, x4, stringsAsFactors = FALSE)

# matrix of distances between observations in xx
dx <- gower.dist(xx)
head(dx)

# matrix of distances between first obs. in xx
# and the remaining ones
gower.dist(data.x=xx[1:6,], data.y=xx[7:10,], var.weights = c(1,2,5,2))

