nsga2_optim = function(popsize=0, generations=0)
{
	function(case) {
		if (popsize == 0) { popsize = nrow(case$parameters) * 12 }
		if (generations == 0) { generations = nrow(case$objectives) * 100 }
		cat("Running NSGA II (",generations,"x",popsize,"): ",sep="")
		progress.bar(generations+2)
		ret = nsga2(
			function(x) {
				ret=case$fun(x);
				progress.tick();
				t(ret)
			},
			idim=nrow(case$parameters),
			odim=nrow(case$objectives),
			cdim=0,
			lower.bounds=case$parameters$lower,
			upper.bounds=case$parameters$upper,
			popsize=popsize,
			generations=generations,
			vectorized = TRUE
		)
		progress.finish()
		ret
	}
}

nsga2_optim_full = function(case) {
	cat("Running NSGA II: ")
	progress.bar(102)
	totpop = NULL
	ret = nsga2_vec(
		function(x) {
			ret=case$fun(x);
			totpop <<- rbind(totpop,cbind(x,ret))
			progress.tick();
			ret
		},
		idim=nrow(case$parameters),
		odim=nrow(case$objectives),
		cdim=0,
		lower.bounds=case$parameters$lower,
		upper.bounds=case$parameters$upper
	)
	progress.finish()
	tab = data.frame(totpop)
        names(tab) = c(case$parameters$name, case$objectives$name)
        w = my.paretoInfo(tab[,case$objectives$name])
        pareto = tab[is.infinite(w),]
        ret = list(
                par=as.matrix(pareto[,case$parameters$name]),
                value=as.matrix(pareto[,case$objectives$name]),
                tab=tab
        )
        ret
}
