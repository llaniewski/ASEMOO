#!/home/llaniewski/bin/RS

####################################################################
#                                                                  #
#  ASEMOO: ASynchronious Efficient Multi-Objctive Optimization     #
#                                                                  #
#  optim.R : R script to generate new designs                      #
#                                                                  #
####################################################################


#----------------------- Libraries --------------------------------#

library(optparse,quietly = T)
library(mco,quietly = T)
#library(DiceKriging,quietly = T)
library(fields,quietly = T)
library(lhs,quietly = T)

#----------------------- Options ----------------------------------#

source(paste(opt$libpath,"opt_fun.R",sep="/"))
source(paste(opt$libpath,"nsga2_vec/mco_vec.R",sep="/"),chdir=T)
source(paste(opt$libpath,"utils.R",sep="/"))
source(paste(opt$libpath,"models.R",sep="/"))
source(paste(opt$libpath,"makemodels.R",sep="/"))
source(paste(opt$libpath,"sampling.R",sep="/"))
source(paste(opt$libpath,"gen.new.R",sep="/"))
source(paste(opt$libpath,"plots.R",sep="/"))
source(paste(opt$libpath,"optims.R",sep="/"))
source(paste(opt$libpath,"pareto.R",sep="/"))

#models=GGModels
#models$init(paste(opt$libpath,"GGM2",sep="/"))

#----------------------- Preprocessing ----------------------------#
lg = data.frame(matrix(NA,nrow=1,ncol=0))

optim.asemoo = function(case, fun, opt=list(), maxiter = 19, models = DiceKrigingModels, sub.optimize=nsga2_optim()) {

	models$init()

	for (n in c("constraints", "objectives", "parameters"))
	{
		case[[n]]$name = as.character(case[[n]]$name)
		rownames(case[[n]]) = case[[n]]$name
	}

	tab = data.frame(design=c())

	#----------------------- newdesign generation ---------------------#
	comment = "Good"
	opt$count = 1
	while (nrow(tab) < maxiter) {
		doe = nrow(case$parameters)*2 + 2
		opt$sd_scale = smart_scale(nrow(tab), maxiter, doe)

		newdesign = try(gen.new(tab=tab, case=case, models=models, opt=opt, doe=doe, sub.optimize=sub.optimize))
		if ("try-error" %in% class(newdesign)) {
			warning("Encountered error in gen.new. quiting optimization")
			comment="Error"
			break;
		}
		newdesign[,c(case$objectives$name, case$constraints$name)] = as.numeric(NA)
		newdesign$state = "running"
		if (nrow(tab) > 0) {
			tab = rbind(tab,newdesign)
		} else {
			tab = newdesign
		}
		sel = apply( is.na(tab[,c(case$objectives$name, case$constraints$name),drop=F] ), 1, any)
		tab[sel, case$objectives$name] = fun(tab[sel,case$parameters$name,drop=F])
		tab$state = "done"
		save(tab, file=paste(opt$actual,paste(opt$casename, "Save", "Rdata", sep="."), sep="/"))
#		Sys.sleep(5)
	}
	list(tab=tab, comment=comment)
}
