#!/bin/echo Run throu newdesign.R. This file's name is 

plot_pareto = function(case, estpareto, pareto, mins, Y)
{
        estpareto=estpareto[order(estpareto[,1]),]  
        if (ncol(estpareto)==3) 
        {
          plot(estpareto,  
        #        type="s",
                main="Pareto fronts (ehvi-sample)",
                xlab=case$objectives$name[1],  
                ylab=case$objectives$name[2],  
                xlim=range(pareto[,1],estpareto[,1],mins[1]),
                ylim=range(pareto[,2],estpareto[,2],mins[2]))
          maxs=par("usr")[c(1,3)]
          points(estpareto[,1:2])
          tmp=rbind(c(maxs[1],mins[2]),pareto[order(pareto[,1]),1:2],c(mins[1],maxs[2]))
          lines(tmp,type="s",col=3)
          points(pareto,col=3)
          plot(estpareto[,3],
                type="p",
                main="EHVI",
                xlab="Pareto points",
                ylab="-EHVI")

          plot(Y,  
        #        type="p",
                main="Sample Points",
                xlab=case$objectives$name[1],
                ylab=case$objectives$name[2],
                xlim=range(pareto[,1],Y[,1],mins[1]),
                ylim=range(pareto[,2],Y[,2],mins[2]))
          lines(tmp,type="s",col=3)
          points(pareto,col=3)
        }
        if (ncol(estpareto)==2)
        {
          plot(estpareto,
                type="s",
                main="Pareto fronts (ehvi-sample)",
                xlab=case$objectives$name[1],
                ylab="EHVI",
                xlim=range(pareto[,1],estpareto[,1],mins[1]),
                ylim=range(estpareto[,2]))
          maxs=par("usr")[c(1,3)]
          points(estpareto[,1:2])
        #  tmp=rbind(c(maxs[1],mins[2]),pareto[order(pareto[,1]),1:2],c(mins[1],maxs[2])
        #  lines(tmp,type="s",col=3)
          abline(v=pareto[1])
        }
}
