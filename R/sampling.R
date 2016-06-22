#!/bin/echo Run throu newdesign.R. This file's name is 

EHVI_sampling = function(tab, case, models, opt, sub.optimize=nsga2_optim())
{
       #----------------------- Selection and transformation -------------#

        for (i in 1:nrow(case$objectives))
        {
                x = case$objectives[i,]
                if (x$name %in% names(tab))
                {
                        if (! x$minimize)
			{
                                nn=paste("neg.",x$name,sep="")
                                tab[,nn] = - tab[,x$name]
 				case$objectives$name[i] = nn
 				case$objectives$threshold[i] = -x$threshold
 				case$objectives$minimize[i] = TRUE
                        }
                } else {
                        stop("Objective ",x$name," not in data (",opt$data,")\n")
                }
        }

        if (nrow(case$constraints) > 0) {
        for (i in 1:nrow(case$constraints))
        {
                x = case$constraints[i,]
                if (x$name %in% names(tab))
                {
                         nn=paste("const.",x$name,sep="")
                         tab[,nn] = ifelse(x$above,1,-1)*(tab[,x$name] - x$threshold)
 			case$constraints$name[i] = nn
 			case$constraints$threshold[i] = 0
 			case$constraints$above[i] = TRUE
                } else {
                        stop("Constraint ",x$name," not in data (",opt$data,")\n")
                }
        }}

        K = make.models(
		tab = tab,
		to.model=c(case$objectives$name, case$constraints$name),
		models=models,
		parameters=case$parameters$name,
		case=case,
		opt=opt
	)

        selected = (tab$state == "done") & (!apply(is.na(tab[,case$objectives$name,drop=F]), 1, any))

        X = tab[selected,case$parameters$name,drop=F]
        Y = tab[selected,case$objectives$name,drop=F]

	if (nrow(case$constraints) >= 1) {
		YC = tab[selected,case$constraints$name,drop=F]
		if ("margin" %in% names(case$constraints)) {
			YC.m = matrix(case$constraints$margin,nrow(YC),ncol(YC),byrow=T)
			YC.m[is.na(YC.m)]=0
		} else {
			YC.m = 0
		}
		YC.s = apply(YC >= -YC.m,1,all)
		X = X[YC.s,,drop=F]
		Y = Y[YC.s,,drop=F]
	}
	mins = case$objectives$threshold
	if (nrow(Y) > 1) {
		Y = Y[,,drop=T]
		pareto=paretoFilter(as.matrix(Y))
        	if (! is.matrix(pareto)) pareto=t(pareto);
	} else {
		pareto = Y
	}
        cat("Actual sample Pareto front\n");
        print(pareto)

        if (is.null(opt$mins_correction)) opt$mins_correction = F;
        cat("Actual mins:",mins,"\n");
        if (opt$mins_correction && (nrow(pareto)>1)) {
		newmins = lower.n(pareto)
		mins = pmin(mins,newmins)
		rm(newmins)
		cat("Corrected mins:",mins,"\n")
        } else cat("No mins correction\n");

	if (nrow(pareto) >= 1) {
	        pareto.sel = pareto > matrix(mins,nrow(pareto),ncol(pareto),byrow=T)
	        pareto.sel = apply(pareto.sel,1,any)
	        cat("Pareto points in threshold:",sum(!pareto.sel),"in,",sum(pareto.sel),"outside\n")
	        pareto = pareto[!pareto.sel,,drop=F]
	}

        if (nrow(pareto) < 1) { warning("Empty pareto - can be an error\n"); pareto = t(mins); }

        tmp=data.frame(t(mins))
        names(tmp) = paste("mins",case$objectives$name,sep=".")
        lg=cbind(lg, mins_correction=opt$mins_correction,tmp)

        ehvi_debug_i = 0
	EHVI_crit = function(nX,verb=F)
        {
         ehvi_debug_i <<- ehvi_debug_i + 1
	 k = nrow(case$objectives)
         nS = sapply(c(case$objectives$name, case$constraints$name), function(i){
          pred = models$predict(K[[i]], nX)
          cbind(pred$mean, pred$sd*opt$sd_scale)
         })
         nY = nS[1:nrow(nX),,drop=F]
         nS = nS[1:nrow(nX)+nrow(nX),,drop=F]
         ehvi= 
          apply(
           cbind(
            nY[,case$objectives$name],
            nS[,case$objectives$name]
           ),
           1,
           function(x){
            EHVI(
             mins,
             pareto,
             x[1:k],
             x[1:k+k]
            )
           }
          );
         if (length(case$constraints$name)>0) {
                ps = nY[,case$constraints$name,drop=F];
                ps[]=0;
                ps = pnorm(nY[,case$constraints$name,drop=F],mean=ps,sd=nS[,case$constraints$name,drop=F])
                psp = apply(ps,1,prod)
         } else {
                psp = 1
         }
         if (verb) {
                colnames(nY) = paste(colnames(nY),"mean",sep=".")
                colnames(nS) = paste(colnames(nS),"sd",sep=".")
                if (exists("ps")){
                        colnames(ps) = paste(colnames(ps),"prob",sep=".")
                        ret=cbind(nY,nS,ps,ehvi,crit=-psp*ehvi)
                } else {
                        ret=cbind(nY,nS,ehvi,crit=-psp*ehvi)
                }
         } else {
                ret=cbind(nY[,case$objectives$name],-psp*ehvi)
#                points(ret)
         }
         cat("EHVI step",ehvi_debug_i,"min:",min(ret),"\n");
         ret
        }
	EHVI_crit_ = function(nX,verb=F)
        {
                r = lapply(1:nrow(case$objectives), function(i) {
                        w = rbind(mins,pareto)[,i]
                        ret = models$predict(K[[case$objectives$name[i]]], nX)
                        list(ei=-EI(
                                matrix(w,nrow(nX),length(w),byrow=T),
                                matrix(ret$mean,nrow(nX),length(w)),
                                matrix(ret$sd,nrow(nX),length(w))
                        ), mean=ret$mean, sd=ret$sd)
                })
                rn = names(r[[1]])
                r = lapply(seq_along(r[[1]]), function(var) {
                        do.call(cbind,
                                lapply(seq_along(r), function(obj) {
                                        r[[obj]][[var]]
                                })
                        )
                })
                names(r) = rn
                nY = r$mean
                ehvi = sapply(1:nrow(nX), function(i) {
                        np = matrix(r$ei[i,],nrow(pareto)+1)
                        prod(np[1,])-dominatedHypervolume(np[-1,,drop=F],ref=np[1,])		
                })
                psp=1
                if (verb) {
			nS = r$sd
                        colnames(nY) = paste(case$objectives$name,"mean",sep=".")
                        colnames(nS) = paste(case$objectives$name,"sd",sep=".")
                        if (exists("ps")){
                                colnames(ps) = paste(colnames(ps),"prob",sep=".")
                                ret=cbind(nY,nS,ps,ehvi,crit=-psp*ehvi)
                        } else {
                                ret=cbind(nY,nS,ehvi,crit=-psp*ehvi)
                        }
                 } else {
                        ret=cbind(nY,-psp*ehvi)
                 }
                 ret
        }

	if (is.null(opt$sd_scale)) opt$sd_scale = 1;
        cat("Using sd_scale:",opt$sd_scale,"\n")
        lg=cbind(lg,sd_scale=opt$sd_scale)

        cat("Running NSGA II for EHVI\n")

#        plot(NA,  
#                type="s",
#                main="Pareto fronts (ehvi-sample)",
#                xlab=case$objectives$name[1],  
#                ylab=case$objectives$name[2],  
#                xlim=range(0.60,0.62),
#                ylim=range(-0.01,0))

#        r = nsga2_vec(
#                EHVI_crit,
#                idim=nrow(case$parameters),
#                odim=nrow(case$objectives)+1,
#                cdim=0,
#                lower.bounds=case$parameters$lower,
#                upper.bounds=case$parameters$upper
#        )
	ncase = case
	ncase$fun = EHVI_crit
	ncase$objectives = rbind(ncase$objectives, data.frame(
		name="EHVI",
		minimize=T,
		threshold=0
	))
	ncase$constraints = data.frame()
#	save(K, EHVI_crit, file="test.Rdata")
		
	r = sub.optimize(ncase)
        colnames(r$par) = case$parameters$name

#        plot_pareto(case=case, pareto=pareto, mins=mins, estpareto=r$val, Y=Y)

        #cat("NSGA II result:\n")
        #print(r)
        #cat("---------------\n")

        ind = which.min(r$val[,nrow(case$objectives)+1])
        cat("Selected index from set:",ind,"\n")
        cat("Parameters of new design:\n")
        print(r$par[ind,])
        cat("Calulating full EHVI values:\n")
        ret = EHVI_crit(r$par, verb=T)
# --------------------- plots -----------------------
	k = nrow(case$objectives)
if (k>1) {
	par(mfrow=c(k-1,k-1))
	for (j in 1:(k-1)) {
		for (i in 2:k) {
			if (j < i) {

				plot(ret[,c(i,j)],
					xlim=range(mins[i], min(ret[,i], ret[,i]-ret[,i+k], mins[i], pareto[,i])),
					ylim=range(mins[j], min(ret[,j], ret[,j]-ret[,j+k], mins[j], pareto[,j])),
					pch=4,cex=0.5
				)
#	save(ret, file="ret.Rdata")
				points.sd(ret[,i],ret[,j],ret[,i+k],ret[,j+k],pch=4,cex=0.5)
			        maxs=par("usr")[c(1,3)]
			        tmp=rbind(c(maxs[1],mins[j]),pareto[order(pareto[,i]),c(i,j)],c(mins[i],maxs[2]))
			        lines(tmp,type="s",col=3, lwd=3)
			        points(pareto[,c(i,j)],col=3)
				points.sd(ret[ind,i],ret[ind,j],ret[ind,i+k],ret[ind,j+k],pch=16,col=2,cex=1,lwd=3,lty=1)
			} else {
				plot(0)
			}
		}
	}
}
        cat("Estimated values in new desing:\n")
        print(ret[ind,])

        cat("EHVI in new desing:",r$val[ind,nrow(case$objectives)+1],"\n")

        design_log_file = paste(opt$actual,"newdesign_log.csv",sep="/")
        if (file.exists(design_log_file))
        {
#                design_log = read.csv(design_log_file)
        } else {
                design_log = NULL
        }
       #tmp= data.frame(t(r$val[ind,1:nrow(case$objectives)]))
       #names(tmp) = paste("estim",case$objectives$name,sep=".")

        tmp = data.frame(t(ret[ind,]))
        lg=cbind(lg, ehvi.selected=ind, tmp)
#        cat("Newdesign log:\n")
#        print(lg)
        design_log <<- merge(design_log,lg,all=T)
        write.csv(design_log, file=design_log_file,row.names=F)

        newparams=r$par[ind,,drop=F]
	newparams
}
