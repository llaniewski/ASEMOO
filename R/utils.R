#!/bin/echo some utils. file:

# Some utils

source.to.list = function(file, ...)
{ local( {
   source(file, local=T, ...)
   as.list(environment())
  } )
}


my.case = function(x) {
	case = list()
	case$parameters = data.frame(
		name=as.character(x$labin),
		lower=x$bounds[,1],
		upper=x$bounds[,2],
		value=NA,
                stringsAsFactors=F
	)
	case$name = x$name
	case$objectives = data.frame(
		name=as.character(x$labout),
		threshold=1,
		minimize=T,
		stringsAsFactors=F
	)
	case$constraints = data.frame(
	        name=c(),
	        above=c(),
	        threshold=c(),
                stringsAsFactors=F
	)
	case$fun=x$fun
	case
}

points.sd = function(x,y,x.sd,y.sd,lty=2,col=1,...) {
	points(x,y,col=col,...)
	a=seq(0,2*pi,len=31)[-1]
	for (i in 1:length(x)) {
		polygon(x[i]+sin(a)*x.sd[i],y[i]+cos(a)*y.sd[i],lty=lty, border=col, ...)
	}
}

smart_scale = function(n,k,doe) {
	1 - pnorm(n, mean=pmax(doe*1.7,(doe+k)/2), sd=doe/3)
}

progress.bar = function() {
	progress.len = 1
	progress.max = 20
	progress.maxlen = 40
	progress.ticks = 0
	clock = "-\\|/"
	
	list(
		bar = function(n)
		{
			progress.len <<- n
			progress.max <<- n
			if (progress.len > progress.maxlen) progress.len <<- progress.maxlen
			n = progress.len+1
			progress.ticks <<- 0
			cat(paste(c("[",rep("-",n),"]",rep("\b",n+1)),collapse=""))
		},
		tick = function() {
			w1 = progress.ticks / progress.max * progress.len
			progress.ticks <<- progress.ticks + 1
			w2 = progress.ticks / progress.max * progress.len
			if ((floor(w1) < floor(w2)) && (w2 <= progress.len)) {
				 cat("#");
			} 
			i = (progress.ticks %% 4) + 1
			cat(substr(clock, i, i), "\b", sep="")
		},
		fail = function() cat("X"),
		finish = function() cat("#\n")
	}
}	
