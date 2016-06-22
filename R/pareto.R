require(mco)

my.paretoInfo = function(x) {
	x = as.matrix(x)
	sapply(1:nrow(x), function(j) {
		r.all = rep(TRUE,nrow(x))
		r.any = rep(FALSE,nrow(x))
		for (i in 1:ncol(x)) {
			r.all = r.all & (x[j,i] >= x[,i])
			r.any = r.any | (x[j,i] > x[,i])
		}
		min(which(r.all & r.any))
	})
}

my.paretoFilter = function(x) {
	w = my.paretoInfo(x)
	x[is.infinite(w),]
}

my.paretoApply = function(x, fun, simplify=T, ...)
{
	w = my.paretoInfo(x)
	sel = w>1:length(w)
	sel = which(sel)
	nsel = rep(FALSE,length(w))
	ret = sapply(sel, function(i) {
		nsel[1:i] = w[1:i] > i
		fun(x[nsel,,drop=F], ...)
	},simplify=simplify)
	nsel = rep(0,length(w))
	nsel[sel]=1
	nsel = cumsum(nsel)
	ret[nsel]
}

my.HV = function(x,ref) {
	if (nrow(x) > 0) {
		dominatedHypervolume(x,ref)
	} else {
		0
	}
}

my.progressiveHV = function(x,ref) {
	my.paretoApply(x, my.HV, ref=ref)
}
