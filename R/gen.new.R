#!/home/llaniewski/bin/RS

####################################################################
#                                                                  #
#  ASEMOO: ASynchronious Efficient Multi-Objctive Optimization     #
#                                                                  #
#  gen.new.R : R script to generate new designs                    #
#                                                                  #
####################################################################


gen.new = function(tab, case, models, opt, ...)
{
	parsel = is.na(case$parameters$value)
	parameters  = case$parameters$name[parsel]
	d = sum(parsel)
	cat("Number of observations:", nrow(tab), "\n");
	old_krig = F
	lg=cbind(lg,old_krig=old_krig)

	new.case = case
	new.case$parameters = case$parameters[parsel,]
	newdesign = gen.new.parsel(tab, new.case, models, opt, ...)

	for ( i in which(!parsel) )
		newdesign[,case$parameters$name[i]]=case$parameters$value[i]

	newdesign
}


gen.new.parsel = function(tab, case, models, opt, doe.min, sub.optimize=nsga2_optim())
{
	d = nrow(case$parameters)
	absolute.min = d*2+2
	if (missing(doe.min)) {	doe.min = absolute.min }
	if (doe.min < absolute.min) { warning("DoE should be above 2*dim+2"); doe.min=absolute.min; }

	newdesign=NULL

	if (nrow(tab) == 0) {
		ndone = 0
	} else {
		ndone = sum(!apply(is.na(tab[,case$objectives$name,drop=F]),1,any))
	}
	cat(paste(nrow(tab),"designs in table.",ndone,"done.\n"))

	#----------------------- Sampling -------------------------------#
	if (ndone >= absolute.min) {
#		pdf(file=paste(opt$actual,"plots.pdf",sep="/"))
		for (n.design in 1:opt$count)
		{ 
			newparams = EHVI_sampling(tab=tab, case=case, models=models, opt=opt, sub.optimize=sub.optimize)
			newparams = data.frame(newparams)
			newdesign = rbind(newdesign,newparams)
			newparams$state="running"
			tab = merge(tab,newparams,all=T)
		}
		newdesign$comment="ASEMOO generated"
	} else {
		cat("Number of (done) designs below absolute minimum:", ndone, "<", absolute.min,"\n");
	}

	#----------------------- DOE --------------------------------------#
	if (is.null(newdesign)) {
		cat("Starting DOE ...\n")
		k=doe.min
		if (k < opt$count) k=opt$count
		if (nrow(tab) == 0)
		{
			p = improvedLHS(k, d)
		} else {
			p = t(as.matrix(tab[,case$parameters$name]))
			p = (p-case$parameters$lower)/(case$parameters$upper-case$parameters$lower)
			p = t(p)
			if (ndone != 0) {
				ratio = ndone/nrow(tab)
				nk=ceiling((absolute.min-ndone)/ratio)
				cat("We need additional",absolute.min-ndone,"designs. Our succress rate is",sprintf("%2.0f%%",100*ratio),", so we are rxtending DoE by",nk,"designs\n")
				if (nk<1) stop("Something's wrong");
				p = augmentLHS(p,nk)
				p = p[nrow(tab)+1:nk,,drop=F]
			} else {
				stop("No designs from DoE calculated")
			}
		}

		cat("LHS table:\n")
		print(p)
		newdesign = t(p) * ( case$parameters$upper-case$parameters$lower)  + case$parameters$lower
		newdesign = data.frame(t(newdesign))
		names(newdesign)=case$parameters$name
		print(newdesign);
		cat("End of DOE\n")
		newdesign=data.frame(newdesign)
		newdesign$comment="Design of Experiment"
	}

	if (is.null(newdesign)) stop("No new designes genereted. Wierd.")

	#----------------------- Plot results -----------------------------#

#	pdf("newdesign_plots.pdf")
#	if (nrow(tab) > 0) plot(tab[,case$parameters$name],pch=ifelse(tab$state=="running", 16, 1), col=ifelse(tab$state=="running",3,1))
#	plot(newdesign[,case$parameters$name])
#	dev.off()

	newdesign
}
