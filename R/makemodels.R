#!/bin/echo Run throu newdesign.R. This file's name is 

make.models = function(tab, parameters, to.model, models, case, opt) {
        K = list()
	to.model.list = as.list(to.model)
	names(to.model.list) = to.model

	pdist = pmax(
			case$parameters$lower - t(tab[,case$parameters$name]),
			0,
			t(tab[,case$parameters$name]) - case$parameters$upper
		)/(case$parameters$upper-case$parameters$lower)
	pdist = sqrt(colSums(pdist**2))
	min.o = nrow(case$parameters)*3;
	K = lapply(to.model.list, function(i) {
                cat("Modeling function",i,":\n");
                selected.o = (tab$state == "done") & (!is.na(tab[,i]))
		pdist.o = pdist[selected.o]
		sel.p = pdist.o == 0
		if (sum(sel.p) < min.o) {
			sel.p[order(pdist.o)[1:min.o]] = T
		} 
		selected = selected.o
		selected[]=F
		selected[selected.o] = sel.p
		selected.r = pdist == 0
                X = tab[selected,parameters,drop=F]
                Y = tab[selected,i,drop=T]
                cat("   Number of accepted observations:",sum(selected),"\n");
                Krig.file = paste(opt$actual,paste(opt$casename,models$acr,i,"Rdata",sep="."),sep="/")
                old_krig = (nrow(X)>400)
                old_start = FALSE
#                if (file.exists(Krig.file))
#                {
#                        cat("Loading old model ...\n");
#                        load(Krig.file)
#                        if (nrow(X) == models$getrows(Kr))
#                        {
#                                cat("same number of cases. I'll keep the model\n");
#                                old_krig=T;
#                        }
#                        old_start = FALSE
#                }
                if (old_krig & file.exists(Krig.file))
                {
                        cat("   Applying old Kriging model\n");
                        #load(Krig.file);
                        Kr = models$oldmodel(Kr,X=X,Y=Y,parameters=case$parameters[parameters,])
                } else {
                        cat("   Constricting Kriging model:\n");
                        if (old_start)
                        { cat("   Starting from the old model\n");
                          Kr = models$oldmodel(Kr,X=X,Y=Y,parameters=case$parameters[parameters,],newfit=T)
                        } else {
                          Kr = models$makemodel(X=X,Y=Y,parameters=case$parameters[parameters,])
                        }
                        save(Kr, file=Krig.file)
                }

                if (any(!selected & selected.r)) {
                        cat("   Adding missing observations:",sum(!selected & selected.r),"\n")
                        nX = tab[!selected & selected.r, parameters,drop=F]
                        X  = rbind(X,nX)
                        nY = models$predict(Kr, nX)$mean
                        tab[!selected  & selected.r,i] = nY
                        Y = c(Y,nY)
                        cat("   Sample set was augmented reconstructing Kriging model\n");
                        print(Y)
                        Kr = models$oldmodel(Kr,X=X,Y=Y,parameters=case$parameters[parameters,])
                }
                Kr;
        })
	K
}
