

skal = function(x,y){ apply(x*y,1,sum)}
skalt = function(x,y){ apply(x*y,2,sum)}

plot.stairs = function(pareto,...)
{
	points(pareto,...)
	lines(c(rep(pareto[,1],each=2),pareto[nrow(pareto),1]+10),c(pareto[1,2]+10,rep(pareto[,2],each=2)),...)
}

plot.stairs3d = function(pareto,z=0,...)
{
	points3d(pareto[,1],pareto[,2],z,...)
	lines3d(c(rep(pareto[,1],each=2),pareto[nrow(pareto),1]+0.1),c(pareto[1,2]+0.1,rep(pareto[,2],each=2)),0,...)
}


EI = function(a,x,sigma)
{
	pnorm(a,mean=x,sd=sigma)*(x-a) - dnorm(a,mean=x,sd=sigma)*(sigma**2)
}


myEI_old = function (pareto,x,sds,minx,miny)
{
	f = function(r) { pnorm(r,mean=x[,1], sd=sds[,1]) }
	pn = t(apply(pareto[,1,drop=F],1,f))
	pn = rbind(pn[1,],diff(pn),1-pn[nrow(pn),])
	f = function(r) { EI(r, x[,2], sds[,2]) }
	ei = t(apply(rbind(minx,pareto[,2,drop=F]),1,f))
	re = skalt(pn,ei)
	f = function(r) { pnorm(r,mean=x[,2], sd=sds[,2]) }
	pn = t(apply(pareto[,2,drop=F],1,f))
	pn = rbind(1-pn[1,],-diff(pn),pn[nrow(pn),])
	f = function(r) { EI(r, x[,1], sds[,1]) }
	ei = t(apply(rbind(pareto[,1,drop=F],miny),1,f))
	cbind(skalt(pn,ei),re)
}


myEI = function (pareto,x,sds,minx,miny)
{
	f = function(r) { pnorm(r,mean=x[,1], sd=sds[,1]) }
	pn = t(apply(rbind(pareto[,1,drop=F],minx),1,f))
	pn = rbind(pn[1,],diff(pn))
	f = function(r) { EI(r, x[,2], sds[,2]) }
	ei = t(apply(rbind(miny,pareto[,2,drop=F]),1,f))
	re = skalt(pn,ei)
	f = function(r) { pnorm(r,mean=x[,2], sd=sds[,2]) }
	pn = t(apply(rbind(miny,pareto[,2,drop=F]),1,f))
	pn = rbind(-diff(pn),pn[nrow(pn),])
	f = function(r) { EI(r, x[,1], sds[,1]) }
	ei = t(apply(rbind(pareto[,1,drop=F],minx),1,f))
	cbind(skalt(pn,ei),re)
}

myEI1 = function (pareto,x,sds,minx,miny)
{
	pn = pnorm(pareto[,1],mean=x[1],sd=sds[1])
	pn = c(pn[1],diff(pn),1-pn[length(pn)])
	ei = EI(c(minx,pareto[,2]),x[2],sds[2])
	re = crossprod(pn,ei)
	pn = pnorm(pareto[,2],mean=x[2],sd=sds[2])
	pn = c(1-pn[1],-diff(pn),pn[length(pn)])
	ei = EI(c(pareto[,2],miny),x[1],sds[1])
	c(crossprod(pn,ei),re)
}



myEHV = function (pareto,x,sds,minx,miny)
{
	f = function(r) { EI(r, x[,1], sds[,1]) }
	pn = t(apply(rbind(pareto[,1,drop=F],minx),1,f))
	pn = rbind(pn[1,,drop=F],diff(pn))
	f = function(r) { EI(r, x[,2], sds[,2]) }
	ei = t(apply(rbind(miny,pareto[,2,drop=F]),1,f))
	skalt(pn,ei)	
}

myEHV1 = function (pareto,x,sds,minx,miny)
{

	pn = EI(c(pareto[,1,drop=F],minx), x[,1], sds[,1])
	pn = c(0,diff(pn))
	f = function(r) { EI(r, x[,2], sds[,2]) }
	ei = EI(c(miny,pareto[,2,drop=F]), x[,2], sds[,2])
	crossprod(pn,ei)	
}

myPI = function (pareto,x,sds,minx,miny)
{
	f = function(r) { pnorm(r,mean=x[,1], sd=sds[,1]) }
	pn = t(apply(rbind(pareto[,1,drop=F],minx),1,f))
	pn = rbind(pn[1,,drop=F],diff(pn))
	f = function(r) { pnorm(r,mean=x[,2], sd=sds[,2]) }
	ei = t(apply(rbind(miny,pareto[,2,drop=F]),1,f))
	skalt(pn,ei)	
}

predict.se.lm = function(obj, x, ...)
{
	rep(sqrt(deviance(obj)/df.residual(obj)),nrow(x))
}

EHVI = function(mins, pareto, x, sd)
{
	minsEI   = -EI(mins, x, sd)
	x  = matrix(x ,nrow=nrow(pareto),ncol=ncol(pareto),byrow=T)
	sd = matrix(sd,nrow=nrow(pareto),ncol=ncol(pareto),byrow=T)
	paretoEI = -EI(pareto, x, sd)
	prod(minsEI)-dominatedHypervolume(paretoEI,minsEI)
}

lower.n = function(x)
{
        n = nrow(x)
        x=apply(x,2,range)
        x=x[2,]+(x[2,]-x[1,])/n
        x
}
