# ejemplo 2 suavizado
install.packages("wavelets")
install.packages("waveslim")
library(wavelets)
library(waveslim)

jj<-5
J<-4
x  <- seq(0, 1, len=2^(jj+1))
n<-length(x)
tt <- seq(0, 1, len=2^(jj+1))
ntt<-length(tt)
fit<-rep(0,len=ntt)
y<-sin(2*pi*tt)+rnorm(ntt,sd=.1)

library(fda.usc)
fy<-fdata(y,tt)   
plot(fy)
## Transformar un argvals normal en uno floor(log2(n))
#a<-approx(x, x, method = "constant",n=n)
for (j in 0:J){
  K<-2^j-1    
  for (k in 0:K){
    fi<-sapply(x, function(x) haar(x,j,k))
    ffi<-fdata(fi,x) #cambiarle la escala!!!
    #      bet<-inprod.fdata(fy,ffi)
    bet=mean(y*fi)
    newfit<-drop(bet)*fi
    fit<-fit+newfit
  }
}   
dady<-padre(x)
#  dady<-padre.mex(2*pi*x)
fpadre<-fdata(dady,x)
alfa<-inprod.fdata(fy,fpadre)
# alfa<-mean(y*dady)
#  fit<-fit+drop(alfa)*fpadre   

fit<-fit+drop(alfa)*fpadre   


# como guardar las bases ortogonales
# range val
plot(fy)
lines(fit,col=4,lwd=2)