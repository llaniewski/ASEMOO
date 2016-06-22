#!/home/llaniewski/bin/RS

####################################################################
#                                                                  #
#  ASEMOO: ASynchronious Efficient Multi-Objctive Optimization     #
#                                                                  #
#  newdesign.R : R script to generate new designs                  #
#                                                                  #
####################################################################

library(optparse,quietly = T)

#----------------------- Options ----------------------------------#

options <- list(
	make_option(c("-d","--data"), "store", default="", help="CSV data file", type="character"),
	make_option(c("-s","--settings"), "store", default="", help="Settings file (.R)", type="character"),
	make_option(c("-o","--output"), "store", default="", help="CSV file with new design", type="character"),
	make_option(c("-a","--actual"), "store", default="actual", help="Where to keep Kriging models and related data", type="character"),
	make_option(c("-n","--count"), "store", default="1", help="Number of designs", type="character"),
	make_option(c("-l","--log"), "store", default="newdesign.log", help="Log file", type="character"),
	make_option(c("-p","--libpath"), "store", default=".", help="Path to newdesign (libs, utils etc)", type="character"),
	make_option(c("-i","--design"), "store", default="", help="Index of the first design to generate", type="character")
)

opt <- parse_args(OptionParser(usage="Usage: newdesing -d datafile.csv -s settingfile.R [-l log] [-o output]",options))
if (opt$data == "") stop("ERROR: No data file");
if (opt$settings == "") stop("ERROR: No settings file");
opt$count = as.numeric(opt$count)
opt$design = as.numeric(opt$design)

#----------------------- Libraries --------------------------------#

library(optparse,quietly = T)
library(mco,quietly = T)
#library(DiceKriging,quietly = T)
library(fields,quietly = T)
library(lhs,quietly = T)

source(paste(opt$libpath,"opt_fun.R",sep="/"))
source(paste(opt$libpath,"nsga2_vec/mco_vec.R",sep="/"),chdir=T)
source(paste(opt$libpath,"utils.R",sep="/"))
source(paste(opt$libpath,"models.R",sep="/"))
source(paste(opt$libpath,"makemodels.R",sep="/"))
source(paste(opt$libpath,"sampling.R",sep="/"))
source(paste(opt$libpath,"gen.new.R",sep="/"))
source(paste(opt$libpath,"plots.R",sep="/"))
source(paste(opt$libpath,"optims.R",sep="/"))

#models=GGModels
#models$init(paste(opt$libpath,"GGM2",sep="/"))
models=DiceKrigingModels
models$init()

#----------------------- Preprocessing ----------------------------#

case=source.to.list(opt$settings)
for (n in c("constraints", "objectives", "parameters"))
{
	case[[n]]$name = as.character(case[[n]]$name)
	rownames(case[[n]]) = case[[n]]$name
}



lg = data.frame(matrix(NA,nrow=1,ncol=0))
tab = try(read.csv(opt$data),TRUE)
if (class(tab) == "try-error") stop("ERROR: Invalid data file: ",opt$data);

if (!("state" %in% names(tab)))
{	cat("Adding \"state\" column to data\n");
	if (nrow(tab) > 0) tab$state="done"
} else {
	tab$state = as.character(tab$state)
}

#----------------------- newdesign generation ---------------------#

newdesign = gen.new(tab=tab, case=case, models=models, opt=opt)

#----------------------- Saving results ---------------------------#

newdesign=data.frame(newdesign)

if (!is.na(opt$design)) {
	cat("Design option set, adding design index to table:", paste(opt$design,opt$design+nrow(newdesign)-1,sep=":"),"\n");
	if (opt$count != nrow(newdesign)) {
		cat(" Note: number of rows different then count set in options (count:",opt$count,"nrow(tab):",nrow(newdesign),") Maybe DoE\n");
	}
	newdesign$design = 1:nrow(newdesign) + opt$design - 1
}

if (opt$output != "")
{
	print(newdesign)
	write.csv(newdesign, file=opt$output, row.names=F)
}
