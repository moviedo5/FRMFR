bcorrdistbin=function(resx,resy){
if (length(resx$Adiag)!=length(resy$Adiag)) stop("The number of elements is not the same for X and Y")
n=length(resx$Adiag)
Uxy=0;Ux=0;Uy=0
for (i in 1:length(resx$A)){
	for (j in 1:length(resy$A)){
	aa=intersect(rownames(resx$A[[i]]$x),rownames(resy$A[[j]]$x))
#	nM[i,j]=length(aa)
#	Mcross[i,j]=length(aa)*resx$A[[i]]$A*resy$A[[j]]$A
	Uxy=Uxy+length(aa)*resx$A[[i]]$A*resy$A[[j]]$A
	}
}
Uxy=Uxy-(n/(n-2))*sum(resx$Adiag*resy$Adiag)
#Vxy=sum(Mcross)
for (i in 1:length(resx$A)){Ux=Ux+resx$A[[i]]$n*resx$A[[i]]$A^2}
Ux=Ux-(n/(n-2))*sum(resx$Adiag^2)
for (i in 1:length(resy$A)){Uy=Uy+resy$A[[i]]$n*resy$A[[i]]$A^2}
Uy=Uy-(n/(n-2))*sum(resy$Adiag^2)
R=Uxy/sqrt(Ux*Uy)
M=n*(n-3)/2
df=M-1
tstat=sqrt(M-1)*R/sqrt(1-R^2)
names(tstat)="T"
names(R)="Bias corrected dcor"
pval=1-pt(tstat,df=df)
method="Binned dcor t-test of independence"
rval=list(statistic=tstat,parameter=df,p.value=pval,estimate=R,method=method)
class(rval)="htest"
return(rval)
}


Abinx=function(x,fdist=metric.lp,delta=-.005,nlim=10000,...){
if (is.vector(x)) x=matrix(x,ncol=1)
if (is.data.frame(x)) x=as.matrix(x)
if (is.matrix(x)) {n=nrow(x);ismat=TRUE} else {n=length(x);ismat=FALSE}
cte=n/(n-1)
dm=numeric(n)
k=ceiling(n/nlim)
bord=1+pmin((0:k)*nlim,n)
lim1=bord[-(k+1)]
lim2=bord[-1]-1
dbord=diff(bord)
dm=rep(0,n)
ndm=rep(0,n)
	if (k==1){
		D=fdist(x,...)
		dm=apply(D,1,mean,na.rm=TRUE)
		ndm=apply(!is.na(D),1,sum)
			} else {
		for (i in 1:length(lim1)) {
			for (j in i:length(lim1)){
				if (ismat){
				D=fdist(x[lim1[i]:lim2[i],,drop=FALSE],x[lim1[j]:lim2[j],,drop=FALSE],...)
				} else {
				D=fdist(x[lim1[i]:lim2[i]],x[lim1[j]:lim2[j]],...)
				}
				if (i==j){
				dm[lim1[j]:lim2[j]]=dm[lim1[j]:lim2[j]]*ndm[lim1[j]:lim2[j]]/(ndm[lim1[j]:lim2[j]]+dbord[i])+apply(D,2,sum,na.rm=TRUE)/(ndm[lim1[j]:lim2[j]]+dbord[i])
				ndm[lim1[j]:lim2[j]]=ndm[lim1[j]:lim2[j]]+dbord[i]
				}else{
				dm[lim1[j]:lim2[j]]=dm[lim1[j]:lim2[j]]*ndm[lim1[j]:lim2[j]]/(ndm[lim1[j]:lim2[j]]+dbord[i])+apply(D,2,sum,na.rm=TRUE)/(ndm[lim1[j]:lim2[j]]+dbord[i])
				ndm[lim1[j]:lim2[j]]=ndm[lim1[j]:lim2[j]]+dbord[i]
				dm[lim1[i]:lim2[i]]=dm[lim1[i]:lim2[i]]*ndm[lim1[i]:lim2[i]]/(ndm[lim1[i]:lim2[i]]+dbord[j])+apply(D,1,sum,na.rm=TRUE)/(ndm[lim1[i]:lim2[i]]+dbord[j])
				ndm[lim1[i]:lim2[i]]=ndm[lim1[i]:lim2[i]]+dbord[j]
				}
			}
		}
}

#dmm=mean(dm)
dmm=sum(dm*ndm)/sum(ndm)
maxA=1.05*max(abs(dmm-2*dm))
if (delta<0) delta=abs(delta)*2*maxA
nbinx=ceiling(1+2*maxA/delta)
lisx=vector("list",length=nbinx)
for (i in 1:nbinx){lisx[[i]]$n=0;lisx[[i]]$A=0}
if (k==1){
#	D=fdist(x,...)  # Computed before?
 	Aij=cte*(sweep(sweep(D,1,dm,"-"),2,dm,"-")+dmm-D/n)
	idx=1+floor((Aij+maxA)/delta)
	tt=table(idx)
	for (l in 1:length(tt)){
	iaux=as.numeric(names(tt)[l])
	lisx[[iaux]]$A=mean(Aij[idx==iaux])
	lisx[[iaux]]$n=tt[l]
	lisx[[iaux]]$x=which(idx==iaux,arr.ind=TRUE)
					}
		}else{
		for (i in 1:length(lim1)){
			for (j in i:length(lim1)){
#			print(paste(lim1[i],lim2[i],lim1[j],lim2[j]))
				if (ismat) {
					D=fdist(x[lim1[i]:lim2[i],,drop=FALSE],x[lim1[j]:lim2[j],,drop=FALSE],...)
						} else {
					D=fdist(x[lim1[i]:lim2[i]],x[lim1[j]:lim2[j]],...)
				}
			Aij=cte*(sweep(sweep(D,1,dm[lim1[i]:lim2[i]],"-"),2,dm[lim1[j]:lim2[j]],"-")+dmm-D/n)
			idx=1+floor((Aij+maxA)/delta)
			tt=table(idx)
			for (l in 1:length(tt)){
			iaux=as.numeric(names(tt)[l])
			xaux=which(idx==iaux,arr.ind=TRUE)
			xaux=sweep(xaux,2,c(lim1[i]-1,lim1[j]-1),"+")
			if (i==j){
		lisx[[iaux]]$A=(lisx[[iaux]]$A*lisx[[iaux]]$n+sum(Aij[idx==iaux]))/(lisx[[iaux]]$n+tt[l])
		lisx[[iaux]]$n=lisx[[iaux]]$n+tt[l]
		lisx[[iaux]]$x=rbind(lisx[[iaux]]$x,xaux)
			} else {
		lisx[[iaux]]$A=(lisx[[iaux]]$A*lisx[[iaux]]$n+2*sum(Aij[idx==iaux]))/(lisx[[iaux]]$n+2*tt[l])
		lisx[[iaux]]$n=lisx[[iaux]]$n+2*tt[l]
		lisx[[iaux]]$x=rbind(lisx[[iaux]]$x,xaux,xaux[,c(2,1)])
			}				
		}
	}
}
}
vn=unlist(sapply(lisx,"[","n"))
nb=length(vn[vn>0])
lisx=lisx[vn>0]
for (idx in 1:length(lisx)){
rownames(lisx[[idx]]$x)=paste(lisx[[idx]]$x[,1],lisx[[idx]]$x[,2],sep="-")
}
return(list(A=lisx,Adiag=cte*(dm-dmm),dm=dm,maxA=maxA))
}



bcordistBD=function(x,y,fdistx=metric.lp,fdisty=metric.lp,par.distx=NULL,par.disty=NULL,nlim=10000){
if (is.vector(x)) x=matrix(x,ncol=1)
if (is.data.frame(x)) x=as.matrix(x)
if (is.matrix(x)) {nx=nrow(x);ismatx=TRUE} else {nx=length(x);ismatx=FALSE}
if (is.vector(y)) y=matrix(y,ncol=1)
if (is.data.frame(y)) y=as.matrix(y)
if (is.matrix(y)) {ny=nrow(y);ismaty=TRUE} else {ny=length(y);ismaty=FALSE}
if (nx!=ny) stop("length(x)!=length(y)")
aax=formals(fdistx)[-c(1:2)]
if (length(par.distx)>0){
	namx=names(par.distx)
	for (i in 1:length(par.distx)){aax[[namx[i]]]=par.distx[[i]]}
	}
aay=formals(fdisty)[-c(1:2)]
if (length(par.disty)>0){
	namy=names(par.disty)
	for (i in 1:length(par.disty)){aax[[namy[i]]]=par.disty[[i]]}
	}

n=nx
cte=n/(n-1)
dmx=rep(0,n)
dmy=rep(0,n)
ndmx=rep(0,n)
ndmy=rep(0,n)

k=ceiling(n/nlim)
bord=1+pmin((0:k)*nlim,n)
lim1=bord[-(k+1)]
lim2=bord[-1]-1

dbord=diff(bord)
	if (k==1){
	    Dx=do.call(fdistx,c(alist(x),aax))
#		Dx=fdistx(x,...=par.distx)
		dmx=apply(Dx,1,mean,na.rm=TRUE)
		ndmx=apply(!is.na(Dx),1,sum)
		Dy=do.call(fdisty,c(alist(y),aay))
#		Dy=fdisty(y,...=par.disty)
		dmy=apply(Dy,1,mean,na.rm=TRUE)
		ndmy=apply(!is.na(Dy),1,sum)
			} else {
		for (i in 1:length(lim1)) {
			for (j in i:length(lim1)){
				if (ismatx){
				Dx=do.call(fdistx,c(alist(x[lim1[i]:lim2[i],,drop=FALSE],x[lim1[j]:lim2[j],,drop=FALSE]),aax))
#				Dx=fdistx(x[lim1[i]:lim2[i],,drop=FALSE],x[lim1[j]:lim2[j],,drop=FALSE],par.distx)
				} else {
				Dx=do.call(fdistx,c(alist(x[lim1[i]:lim2[i]],x[lim1[j]:lim2[j]]),aax))
#				Dx=fdistx(x[lim1[i]:lim2[i]],x[lim1[j]:lim2[j]],par.distx)
				}
				if (ismaty){
				Dy=do.call(fdisty,c(alist(y[lim1[i]:lim2[i],,drop=FALSE],y[lim1[j]:lim2[j],,drop=FALSE]),aay))
#				Dy=fdisty(y[lim1[i]:lim2[i],,drop=FALSE],y[lim1[j]:lim2[j],,drop=FALSE],par.disty)
				} else {
				Dy=do.call(fdisty,c(alist(y[lim1[i]:lim2[i]],y[lim1[j]:lim2[j]]),aay))
#				Dy=fdisty(y[lim1[i]:lim2[i]],y[lim1[j]:lim2[j]],par.disty)
				}
				if (i==j){
				dmx[lim1[j]:lim2[j]]=dmx[lim1[j]:lim2[j]]*ndmx[lim1[j]:lim2[j]]/(ndmx[lim1[j]:lim2[j]]+dbord[i])+apply(Dx,2,sum,na.rm=TRUE)/(ndmx[lim1[j]:lim2[j]]+dbord[i])
				ndmx[lim1[j]:lim2[j]]=ndmx[lim1[j]:lim2[j]]+dbord[i]
				dmy[lim1[j]:lim2[j]]=dmy[lim1[j]:lim2[j]]*ndmy[lim1[j]:lim2[j]]/(ndmy[lim1[j]:lim2[j]]+dbord[i])+apply(Dy,2,sum,na.rm=TRUE)/(ndmy[lim1[j]:lim2[j]]+dbord[i])
				ndmy[lim1[j]:lim2[j]]=ndmy[lim1[j]:lim2[j]]+dbord[i]
				}else{
				dmx[lim1[j]:lim2[j]]=dmx[lim1[j]:lim2[j]]*ndmx[lim1[j]:lim2[j]]/(ndmx[lim1[j]:lim2[j]]+dbord[i])+apply(Dx,2,sum,na.rm=TRUE)/(ndmx[lim1[j]:lim2[j]]+dbord[i])
				ndmx[lim1[j]:lim2[j]]=ndmx[lim1[j]:lim2[j]]+dbord[i]
				dmx[lim1[i]:lim2[i]]=dmx[lim1[i]:lim2[i]]*ndmx[lim1[i]:lim2[i]]/(ndmx[lim1[i]:lim2[i]]+dbord[j])+apply(Dx,1,sum,na.rm=TRUE)/(ndmx[lim1[i]:lim2[i]]+dbord[j])
				ndmx[lim1[i]:lim2[i]]=ndmx[lim1[i]:lim2[i]]+dbord[j]
				dmy[lim1[j]:lim2[j]]=dmy[lim1[j]:lim2[j]]*ndmy[lim1[j]:lim2[j]]/(ndmy[lim1[j]:lim2[j]]+dbord[i])+apply(Dy,2,sum,na.rm=TRUE)/(ndmy[lim1[j]:lim2[j]]+dbord[i])
				ndmy[lim1[j]:lim2[j]]=ndmy[lim1[j]:lim2[j]]+dbord[i]
				dmy[lim1[i]:lim2[i]]=dmy[lim1[i]:lim2[i]]*ndmy[lim1[i]:lim2[i]]/(ndmy[lim1[i]:lim2[i]]+dbord[j])+apply(Dy,1,sum,na.rm=TRUE)/(ndmy[lim1[i]:lim2[i]]+dbord[j])
				ndmy[lim1[i]:lim2[i]]=ndmy[lim1[i]:lim2[i]]+dbord[j]
				}
			}
		}
}

dmmx=sum(dmx*ndmx)/sum(ndmx)
dmmy=sum(dmy*ndmy)/sum(ndmy)
Uxy=0;Ux=0;Uy=0
if (k==1){
#	Dx=fdistx(x,par.distx)  # Computed before?
#	Dy=fdisty(y,par.disty)  # Computed before?
 	Aij=sweep(sweep(Dx,1,dmx,"-"),2,dmx,"-")+dmmx-Dx/n
 	Bij=sweep(sweep(Dy,1,dmy,"-"),2,dmy,"-")+dmmy-Dy/n
	diag(Aij)=dmx-dmmx
	diag(Bij)=dmy-dmmy
	Aij=cte*Aij;Bij=cte*Bij
	Uxy=sum(Aij*Bij)-(n/(n-2))*sum(diag(Aij)*diag(Bij))
	Ux=sum(Aij^2)-(n/(n-2))*sum(diag(Aij^2))
	Uy=sum(Bij^2)-(n/(n-2))*sum(diag(Bij^2))
		}else{
		for (i in 1:length(lim1)){
			for (j in i:length(lim1)){
#			print(paste(lim1[i],lim2[i],lim1[j],lim2[j]))
				if (ismatx) {
					Dx=do.call(fdistx,c(alist(x[lim1[i]:lim2[i],,drop=FALSE],x[lim1[j]:lim2[j],,drop=FALSE]),aax))
#					Dx=fdistx(x[lim1[i]:lim2[i],,drop=FALSE],x[lim1[j]:lim2[j],,drop=FALSE],par.distx)
						} else {
					Dx=do.call(fdistx,c(alist(x[lim1[i]:lim2[i]],x[lim1[j]:lim2[j]]),aax))
#					Dx=fdistx(x[lim1[i]:lim2[i]],x[lim1[j]:lim2[j]],par.distx)
				}
				if (ismaty) {
					Dy=do.call(fdisty,c(alist(y[lim1[i]:lim2[i],,drop=FALSE],y[lim1[j]:lim2[j],,drop=FALSE]),aay))
#					Dy=fdisty(y[lim1[i]:lim2[i],,drop=FALSE],y[lim1[j]:lim2[j],,drop=FALSE],par.disty)
						} else {
					Dy=do.call(fdisty,c(alist(y[lim1[i]:lim2[i]],y[lim1[j]:lim2[j]]),aay))
#					Dy=fdisty(y[lim1[i]:lim2[i]],y[lim1[j]:lim2[j]],par.disty)
				}
			Aij=cte*(sweep(sweep(Dx,1,dmx[lim1[i]:lim2[i]],"-"),2,dmx[lim1[j]:lim2[j]],"-")+dmmx-Dx/n)
			Bij=cte*(sweep(sweep(Dy,1,dmy[lim1[i]:lim2[i]],"-"),2,dmy[lim1[j]:lim2[j]],"-")+dmmy-Dy/n)
			if (i==j){
			diag(Aij)=dmx[lim1[i]:lim2[i]]-dmmx
			diag(Bij)=dmy[lim1[i]:lim2[i]]-dmmy
			Ux=Ux+sum(Aij^2)-(n/(n-2))*sum(diag(Aij*Aij))
			Uy=Uy+sum(Bij^2)-(n/(n-2))*sum(diag(Bij*Bij))
			Uxy=Uxy+sum(Aij*Bij)-(n/(n-2))*sum(diag(Aij*Bij))
			} else {
			Ux=Ux+2*sum(Aij^2)
			Uy=Uy+2*sum(Bij^2)
			Uxy=Uxy+2*sum(Aij*Bij)
			}				
		}
	}
}

R=Uxy/sqrt(Ux*Uy)
M=n*(n-3)/2
df=M-1
tstat=sqrt(M-1)*R/sqrt(1-R^2)
names(tstat)="T"
names(R)="Bias corrected dcor"
pval=1-pt(tstat,df=df)
method="Binned dcor t-test of independence"
rval=list(statistic=tstat,parameter=df,p.value=pval,estimate=R,method=method)
class(rval)="htest"
return(rval)
}


