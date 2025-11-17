# -------------------------------------------------------------------
# BikeSharing dataset
# -------------------------------------------------------------------
# Packages
# library(fda.usc.devel) # old version (until 2025/11/14)
library(fda.frm)   # new version (since 2025/11/14)
library(FRegSigCom)
library(refund)

# -------------------------------------------------------------------
# Simulation 2018–19
# -------------------------------------------------------------------
# load("omel2018-19.rda")
data(omel2018_19)
names(omel2018_19)

ldat <- omel2018_19 
npx <- 3                                                        #4 
npy <- 3                                                        #4 

ttx <- ldat$En$argvals
tty <- ldat$Pr$argvals
nbasis.x <- nbasis.y <- s.n.basis <- t.n.basis <- x.n.basis <- 7 # 11
nk <- nbasis.y
nkk <- c(nbasis.x, nbasis.y)

set.seed(20030101)
nrep <- 10                                                 # nrep <- 20
rr <- matrix(NA, nrep, 7)
colnames(rr) <- c("FLMFR", "FSAMFR", "FKAMFR", "PFR", "FAMM", "LSC", "DISC")

for (i in 1:nrep) { # i <- 1
  print(i)
  ii <- sample(nrow(ldat$df), floor(nrow(ldat$df) * 0.7)) # 0.75
  ldatest <- ldat[ii, row = TRUE]
  ldatpred <- ldat[-ii, row = TRUE]
  
  b.x <- list("En" = create.pc.basis(ldatest$En, 1:npx))
  b.y <- create.pc.basis(ldatest$Pr, 1:npy)
  
  rest <- sum(norm.fdata(ldatpred$Pr - func.mean(ldatest$Pr))^2)
  
  reslm <- fregre.mlm.fr(Pr ~ En, ldatest, basis.x = b.x, basis.y = b.y)
  predlm <- predict(reslm, ldatpred)
  rr[i, 1] <- 1 - sum(norm.fdata(ldatpred$Pr - predlm)^2) / rest
  
  ressam <- fregre.sam.fr(Pr ~ s(En), ldatest,
                          basis.x = b.x, basis.y = b.y)
  predsam <- predict(ressam, ldatpred)
  rr[i, 2] <- 1 - sum(norm.fdata(ldatpred$Pr - predsam)^2) / rest
  
  par.np <- list()
  par.metric <- list(En = list(metric = metric.lp, lp = 2))
  reskam <- fregre.kam.fr(Pr ~ En, ldatest, 
                          par.np = par.np, 
                          par.metric = par.metric)
  predkam <- predict(reskam, ldatpred)
  rr[i, 3] <- 1 - sum(norm.fdata(ldatpred$Pr - predkam)^2) / rest
  
  ydat <- ldatest$Pr$data
  xdat <- ldatest$En$data
  tj <- ldatest$Pr$argvals
  si <- ldatest$En$argvals
  
  # time consuming
  respffpc <- pffr(ydat ~ ffpc(xdat, xind = si, npc.max = npx), yind = tj)
  predffpc <- fdata(predict(respffpc, list(xdat = ldatpred$En$data)), 
                    argvals = tty)
  rr[i, 4] <- 1 - sum(norm.fdata(ldatpred$Pr - predffpc)^2) / rest
  print("pffr1")
  respffnl <- pffr(
    ydat ~ sff(
      xdat,
      xind = si,
      splinepars = list(bs = "ps", m = c(2, 2, 2), k = nkk)
    ),
    yind = tj,
    bs.yindex = list(bs = "ps", k = nk, m = c(2, 1))
    # bs.int = list(bs = "ps", k = 21, m = c(2, 1))
  )
  print("pffr2")
  predffnl <- fdata(predict(respffnl, list(xdat = ldatpred$En$data)), 
                    argvals = tty)
  rr[i, 5] <- 1 - sum(norm.fdata(ldatpred$Pr - predffnl)^2) / rest
  
  XX <- list(Xen = ldatest$En$data)
  XXnew <- list(Xen = ldatpred$En$data)
  t.x.list <- list(Xen = si)
  
  resSCLin <- cv.sigcom(XX, ydat, t.x.list, tj, 
                        s.n.basis = s.n.basis, 
                        t.n.basis = t.n.basis)
  predSCLin <- fdata(pred.sigcom(resSCLin, XXnew), tty)
  rr[i, 6] <- 1 - sum(norm.fdata(ldatpred$Pr - predSCLin)^2) / rest
  print("cv.sigcom")
  resSCNL <- cv.nonlinear(XX, ydat, t.x.list, tj, 
                          s.n.basis = s.n.basis, 
                          x.n.basis = x.n.basis, t.n.basis = t.n.basis)
  predSCNL <- fdata(pred.nonlinear(resSCNL, XXnew), tty)
  rr[i, 7] <- 1 - sum(norm.fdata(ldatpred$Pr - predSCNL)^2) / rest
}

print(apply(rr, 2, mean, na.rm = TRUE))
print(apply(rtim, 2, mean, na.rm = TRUE))

 
# -------------------------------------------------------------------
# Simulation 2008–09
# -------------------------------------------------------------------
# load("omel2008-09.rda")
data("omel2018_19")
ldat <- omel2018_19
npx <- 3                                                        #4 
npy <- 3                                                        #4 

nn <- nrow(ldat$En)
nlag <- 7
nl <- (nlag + 1):nn

ldatm <- ldata(
  df    = ldat$df[nl, ],
  Pr = ldat$Pr[nl],
  ener   = ldat$En[nl],
  ener1  = ldat$En[nl - 1],
  ener7  = ldat$En[nl - nlag],
  price1  = ldat$Pr[nl - 1],
  price7  = ldat$Pr[nl - nlag]
)

ttx <- ldatm$ener$argvals
tty <- ldatm$Pr$argvals
nbasis.x <- nbasis.y <- s.n.basis <- t.n.basis <- x.n.basis <- 7 # 11
nk <- nbasis.y
nkk <- c(nbasis.x, nbasis.y)

set.seed(20030101)
nrep <- 10
rr <- matrix(NA, nrep, 7)
rtim <- matrix(NA, nrep, 7)
colnames(rr) <- c("FLMFR", "FSAMFR", "FKAMFR", "PFR", "FAMM", "LSC", "DISC")

for (i in 1:nrep) {
  print(paste("Repetition:", i))
  
  ii <- sample(nrow(ldatm$df), floor(nrow(ldatm$df) * 0.7)) # 0.75
  ldatest <- ldatm[ii, row = TRUE]
  ldatpred <- ldatm[-ii, row = TRUE]
  
  itime <- Sys.time()
  
  b.x <- list(
    "ener"  = create.pc.basis(ldatest$ener, 1:npx),
    "ener1" = create.pc.basis(ldatest$ener1, 1:npx),
    "ener7" = create.pc.basis(ldatest$ener7, 1:npx),
    "price1" = create.pc.basis(ldatest$price1, 1:npy),
    "price7" = create.pc.basis(ldatest$price7, 1:npy)
  )
  b.y <- create.pc.basis(ldatest$Pr, 1:npy)
  
  rest <- sum(norm.fdata(ldatpred$Pr - func.mean(ldatest$Pr))^2)
  
  # reslm <- fregre.mlm.fr(Pr ~ ener1 + ener7, ldatest, basis.x = b.x, basis.y = b.y)
  reslm <- fregre.mlm.fr(Pr ~ price1 + price7, ldatest, 
                         basis.x = b.x, basis.y = b.y)
  predlm <- predict(reslm, ldatpred)
  rr[i, 1] <- 1 - sum(norm.fdata(ldatpred$Pr - predlm)^2) / rest
  
  itime2 <- Sys.time()
  print(paste("FLMFR:", round(difftime(itime2, itime, units = "mins"), 2)))
  rtim[i, 1] <- difftime(itime2, itime, units = "mins")
  
  # ressam <- fregre.sam.fr(Pr ~ s(ener1) + s(ener7), ldatest, basis.x = b.x, basis.y = b.y)
  ressam <- fregre.sam.fr(Pr ~ s(price1) + s(price7), ldatest, 
                          basis.x = b.x, basis.y = b.y)
  predsam <- predict(ressam, ldatpred)
  rr[i, 2] <- 1 - sum(norm.fdata(ldatpred$Pr - predsam)^2) / rest
  
  itime3 <- Sys.time()
  print(paste("FSAMFR:", round(difftime(itime3, itime2, units = "mins"), 2)))
  rtim[i, 2] <- difftime(itime3, itime2, units = "mins")
  
  # par.np <- list(ener = list(Ker = AKer.norm, type.S = "S.NW"),
  #                ener1 = list(Ker = AKer.norm, type.S = "S.NW"),
  #                ener7 = list(Ker = AKer.norm, type.S = "S.NW"))
  # par.metric <- list(ener = list(metric = metric.lp, lp = 2),
  #                    ener1 = list(metric = metric.lp, lp = 2),
  #                    ener7 = list(metric = metric.lp, lp = 2))
  # reskam <- fregre.kam.fr(Pr ~ ener1 + ener7, ldatest, 
  #                         par.np = par.np, par.metric = par.metric)
  
  par.np <- list(
    Pr = list(Ker = AKer.norm, type.S = "S.NW"),
    price1  = list(Ker = AKer.norm, type.S = "S.NW"),
    price7  = list(Ker = AKer.norm, type.S = "S.NW")
  )
  par.metric <- list(
    Pr = list(metric = metric.lp, lp = 2),
    price1  = list(metric = metric.lp, lp = 2),
    price7  = list(metric = metric.lp, lp = 2)
  )
  reskam <- fregre.kam.fr(Pr ~ price1 + price7, ldatest,
                          par.np = par.np, par.metric = par.metric)
  predkam <- predict(reskam, ldatpred)
  rr[i, 3] <- 1 - sum(norm.fdata(ldatpred$Pr - predkam)^2) / rest
  
  itime4 <- Sys.time()
  print(paste("FKAMFR:", round(difftime(itime4, itime3, units = "mins"), 2)))
  rtim[i, 3] <- difftime(itime4, itime3, units = "mins")
  
  ydat  <- ldatest$Pr$data
  xdat  <- ldatest$ener$data
  xdat1 <- ldatest$ener1$data
  xdat7 <- ldatest$ener7$data
  ydat1 <- ldatest$price1$data
  ydat7 <- ldatest$price7$data
  tj <- ldatest$Pr$argvals
  si <- ldatest$ener$argvals
  
  # respffpc <- pffr(ydat ~ ffpc(xdat1, xind = si, npc.max = npx) +
  #                        ffpc(xdat7, xind = si, npc.max = npx), yind = tj)
  # predffpc <- fdata(
  #   predict(respffpc, list(xdat1 = ldatpred$ener1$data, xdat7 = ldatpred$ener7$data)),
  #   argvals = tty
  # )
  
  respffpc <- pffr(ydat ~ ffpc(ydat1, xind = tj, npc.max = npy) +
                     ffpc(ydat7, xind = tj, npc.max = npy), yind = tj)
  predffpc <- fdata(predict(respffpc, list(ydat1 = ldatpred$price1$data, 
                                           ydat7 = ldatpred$price7$data)), 
                    argvals = tty)
  rr[i, 4] <- 1 - sum(norm.fdata(ldatpred$Pr - predffpc)^2) / rest
  
  itime5 <- Sys.time()
  print(paste("PFR:", round(difftime(itime5, itime4, units = "mins"), 2)))
  rtim[i, 4] <- difftime(itime5, itime4, units = "mins")
  # llist <- list(bs = "ps", m = c(2, 2, 2), k = nkk)
  # respffnl <- pffr(ydat ~ sff(xdat1, xind = si, splinepars = llist) +
  #                        sff(xdat7, xind = si, splinepars = llist)),
  #                  yind = tj, bs.yindex = list(bs = "ps", k = nk, m = c(2, 1)))
  # predffnl <- fdata(
  #   predict(respffnl, list(xdat1 = ldatpred$ener1$data, xdat7 = ldatpred$ener7$data)),
  #   argvals = tty
  # )
  lparam <- list(bs = "ps", m = c(2, 2, 2), k = c(nk, nkk))
  respffnl <- pffr(
    ydat ~ sff(ydat1, xind = tj, splinepars = lparam) +
      sff(ydat7, xind = tj, splinepars = lparam), yind = tj,
      bs.yindex = list(bs = "ps", k = nk, m = c(2, 1))
    # bs.int = list(bs = "ps", k = 21, m = c(2, 1))
  )
  predffnl <- fdata(predict(respffnl, list(ydat1 = ldatpred$price1$data,
                                           ydat7 = ldatpred$price7$data)),
                    argvals = tty)
  rr[i, 5] <- 1 - sum(norm.fdata(ldatpred$Pr - predffnl)^2) / rest
  
  itime6 <- Sys.time()
  print(paste("FAMM:", round(difftime(itime6, itime5, units = "mins"), 2)))
  rtim[i, 5] <- difftime(itime6, itime5, units = "mins")
  
  # XX <- list(Xen1 = ldatest$ener1$data, Xen7 = ldatest$ener7$data)
  # XXnew <- list(Xen1 = ldatpred$ener1$data, Xen7 = ldatpred$ener7$data)
  # t.x.list <- list(Xen1 = si, Xen7 = si)
  
  XX <- list(Yen1 = ldatest$price1$data, Yen7 = ldatest$price7$data)
  XXnew <- list(Yen1 = ldatpred$price1$data, Yen7 = ldatpred$price7$data)
  t.x.list <- list(Yen1 = tj, Yen7 = tj)
  
  resSCLin <- cv.sigcom(XX, ydat, t.x.list, tj, 
                        s.n.basis = s.n.basis, t.n.basis = t.n.basis)
  predSCLin <- fdata(pred.sigcom(resSCLin, XXnew), tty)
  rr[i, 6] <- 1 - sum(norm.fdata(ldatpred$Pr - predSCLin)^2) / rest
  
  itime7 <- Sys.time()
  print(paste("LSC:", round(difftime(itime7, itime6, units = "mins"), 2)))
  rtim[i, 6] <- difftime(itime7, itime6, units = "mins")
  
  resSCNL <- cv.nonlinear(XX, ydat, t.x.list, tj, 
                          s.n.basis = s.n.basis, x.n.basis = x.n.basis, 
                          t.n.basis = t.n.basis)
  predSCNL <- fdata(pred.nonlinear(resSCNL, XXnew), tty)
  rr[i, 7] <- 1 - sum(norm.fdata(ldatpred$Pr - predSCNL)^2) / rest
  
  itime8 <- Sys.time()
  print(paste("DISC:", round(difftime(itime8, itime7, units = "mins"), 2)))
  rtim[i, 7] <- difftime(itime8, itime7, units = "mins")
  
  print(apply(rr, 2, mean, na.rm = TRUE))
  print(apply(rtim, 2, mean, na.rm = TRUE))
}

# Results
print(apply(rr, 2, mean, na.rm = TRUE))
print(apply(rtim, 2, mean, na.rm = TRUE))