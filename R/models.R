#!/bin/echo Run throu newdesign.R. This file's name is 


DiceKrigingModels = list(#----------------------------- DiceKriging Library based model ---------------------------------

name="Kriging model from DiceKriging package"

,acr="DiceKriging"

,init = function(...) {
	library(DiceKriging)
}

,makemodel = function(X,Y, ...) {
	cat("DiceKriging: fitting model: ");
#                Kr = km(~., X, Y, covtype="gauss", lower=c(rep(0.1,ncol(X))), upper=c(rep(20,ncol(X))), nugget.est=T)
		Krls  = list()
		j=0;
		form = as.formula(paste("~.",paste("I(",names(X),"^2)",sep="",collapse="+"),sep="+"))
		print(form)
		progress.bar(30)
		for (i in 1:30) {
#	                Kr = try(km(~., X, Y, nugget.est=F))
			f = textConnection("output_dump","w"); sink(f)
	                Kr = try(km(form, X, Y, nugget.est=F))
#	                Kr = try(km(~., X, Y, nugget.est=F))
			sink(); close(f)
			if (!("try-error" %in% class(Kr))) {
				progress.tick()
#				if (all(Kr@covariance@range.val > Kr@lower) & all(Kr@covariance@range.val < Kr@upper)) {
				if (all(Kr@covariance@range.val > Kr@lower)) {
					j = j + 1;
					Krls[[j]] = Kr;
					if (j >= 1) break;
				} else {
					if (i > 20) {
						j = j + 1;
						Krls[[j]] = Kr;
						break;
					}
				}
			} else {
				progress.fail()
			}
		}
		progress.finish()
		if (j == 0) stop("Fitting the model failed (",i,"try)")
#		nuggets = sapply(Krls, function(Kr) { Kr@covariance@nugget } )
#		sel=which.min(nuggets)
		sel=1
		cat("Selecting",sel,"fit, because of the lowest nugget\n");
		Krls[[sel]]
}

,oldmodel = function(Kr, X,Y, ...) {
	cat("---------------------- DiceKriging: old model --\n");
		form = as.formula(paste("~.",paste("I(",names(X),"^2)",sep="",collapse="+"),sep="+"))
		Kr = km(form, X, Y,
	        	coef.cov=Kr@covariance@range.val,
		        coef.trend=Kr@trend.coef,
		        coef.var=Kr@covariance@sd2,
			nugget=Kr@covariance@nugget)
	Kr
}

,predict = function(obj, X, ...) {
	predict(obj, X, type="UK")
}

,getrows = function(obj) {
	nrow(obj@X)
}
)#-----------------------------------------------------------------------------------------------------------------------


KrigModels = list(#----------------------------- DiceKriging Library based model ---------------------------------

name="Kriging model from fields package"

,acr="Krig"

,init = function(...) {
	library(fields)
}

,makemodel = function(X,Y, ...) {
	cat("---------------------- Krig (fields): make model --\n");
                Kr = Krig(as.matrix(X), as.vector(Y))
	Kr
}

,oldmodel = function(Kr, X,Y, ...) {
	cat("---------------------- Krig (fields): old model --\n");
		print(dim(X))
		print(length(Y))
		print(dim(as.matrix(X)))
                Kr = Krig(as.matrix(X), as.vector(Y))
	Kr
}

,predict = function(obj, X, ...) {
	list( mean = predict(obj,X), sd = predict.se(obj, X) )
}

,getrows = function(obj) {
	nrow(obj$x)
}
)#-----------------------------------------------------------------------------------------------------------------------



GGM.model = function(X, Y, parameters, ...) {
	print(parameters)
	X = as.matrix(X)
	Y = as.vector(Y)
	cat("czesc\n");
	Y.scale = sd(Y)
	mymodel = Kriging(X,Y,
		Mixed(~., X,Y,
			ModelSum(
				Noise(parameters, scale=Y.scale),
				SpatialGauss(parameters, scale=Y.scale, theta.min=log(1e-5), theta.max=log(1e+5))
			)
		)
	)
	print(mymodel$parameters)
	mymodel
}


GGModels = list(#------------------------------------------- aKrig based model------------------------------------------------

name="Kriging model with General Gaussian Models 2"

,acr="GGM2"

,init = function(path=NULL,...) {
#	source("~/Rwork/lib/OR2.R")
	source(paste(path,"models.R",sep="/"))
	source(paste(path,"kriging.R",sep="/"))
	cat("------------- Using aKrig ------------------------\n");
}

,makemodel = function(...) {
	mymodel=GGM.model(...)
	cat("------------- Fitting model ----------------------\n");
	mymodel=mymodel$fit()
	mymodel=mymodel$recalculate()
	print(mymodel$parameters)
	mymodel
}

,oldmodel = function(obj,newfit=F,...) {
	mymodel=GGM.model(...)
	cat("------------- Using old model --------------------\n");
	mymodel$setparameters(obj$parameters$value)
	if (newfit) mymodel=mymodel$fit(p=obj$parameters$value)
	mymodel=mymodel$recalculate()
	mymodel
}

,predict = function(obj, X, ...) {
	predict(obj, X, ...)
}

,getrows = function(obj) {
	nrow(obj$X)
}

)#-----------------------------------------------------------------------------------------------------------------------