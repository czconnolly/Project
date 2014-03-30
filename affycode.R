# script to calculate differential expression and plot raw/normalised data using LIMMA on the microarray data

#install and require the necessary packages
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
biocLite("hgu133plus2.db")
library(ggplot2)
library(reshape2)
library(GEOquery)
library(affy)
library(affyPLM)
library(arrayQualityMetrics)
library(annotate)
library(chicken.db)
library(chickenprobe)
library(genefilter)
library(limma)

.ggboxIt <- function(dfin,pdfname)
{	# makes a boxplot of the data from .summariseIt function
	gp <- ggplot(dfin, aes(x = name, ymin = lwisk, lower = lbox,
													 middle = mid, upper = ubox,
													 ymax = uwisk, fill = name))
  gp <- gp + geom_boxplot(stat="identity",show_guide=F) + coord_flip() 
	gp <- gp + theme(panel.background = element_blank())
	gp <- gp + labs(y = "Log2(Intensity)", x = "")
  ggsave(pdfname) 
}

.summariseIt <- function(dfraw,phenoData)
{ # creates quantile data from a dataframe of depths
  
	dfin <- melt(dfraw)
	colnames(dfin) <- c("probeId","sample","value")
	dfin$group<-factor(dfin$sample)
  stats <- as.data.frame(t(sapply(levels(dfin$sample),
            function(x){
              quantile(log2(as.numeric(dfin$value))[which(dfin$group==x)],
											 prob=c(0.05,0.25,0.5,0.75,0.95),na.rm=T)
            })))
	colnames(stats)<-c("lwisk","lbox","mid","ubox","uwisk")
  stats$name <- rownames(stats)
  rownames(stats) <- NULL
 	return(stats)
}

.getData <- function()
{	#sets working directory
	baseDir <- "/homes/czconnolly/BS32010/Project"
	workingDir=paste0(baseDir)
	# input file names 
	filenames <- c("ROS1-_9.CEL", "ROS1+_10.CEL", "ROS1-_11.CEL",
							"ROS1+_12.CEL", "ROS1-_13.CEL", "ROS1+_14.CEL",
							"ROS1-_15.CEL", "ROS1+_16.CEL")
	# input sample names
	samplenames <- c("ROS1-_9", "ROS1+_10", "ROS1-_11",
  						"ROS1+_12", "ROS1-_13", "ROS1+_14",
							"ROS1-_15", "ROS1+_16")
	#input the treatments to match the samples
	targets <- c("minus", "plus", "minus", "plus", "minus", "plus",
						"minus", "plus")
	#create a table of the filenames,samplenames, and targets
	phenodata<-as.data.frame(cbind(filenames,samplenames,targets))
	write.table(phenodata,paste(workingDir,"phenodata.txt",sep="/")
							,quote=F,row.name=F)
	celRAW <- ReadAffy(celfile.path=workingDir,compress=T,
	
}
#plot the probe intensities
.plotDensity <- function(exps,filename)
{
	pdf(filename)
	# Plot a density vs log intensity histogram for the unnormalised data
	d <- apply(exps,2,
				function(x){
					density(x)
				})
	xmax <- max(sapply(d,function(x)max(x$x)))
	xmin <- min(sapply(d,function(x)min(x$x)))
	ymax <- max(sapply(d,function(x)max(x$y)))
	plot(0,pch='',ylab='',xlab='',
			 xlim=c(xmin,round(xmax+1)),ylim=c(0,ymax))
	lapply(1:length(d),function(x) lines(d[[x]],col=x))
	dev.off()
}

# Perform probe-level metric calculations on the CEL files:
.doPLM <- function(celRAW)
{
	pdf("celRAWqc.pdf")
	#affyPLM is required to interrogate celRMA
	celRAWqc <- fitPLM(celRAW)

	# Create an image of GSM24662.CEL:
	image(celRAWqc, which=1, add.legend=TRUE)

	# Create an image of GSM524665.CEL
	# There is a spatial artifact present
	image(celRAWqc, which=4, add.legend=TRUE)

	# affyPLM also provides more informative boxplots
	# RLE (Relative Log Expression) plots should have
	# values close to zero. GSM524665.CEL is an outlier
	RLE(celRAWqc, main="RLE")

	# We can also use NUSE (Normalised Unscaled Standard Errors).
	# The median standard error should be 1 for most genes.
	# GSM524665.CEL appears to be an outlier on this plot too
	NUSE(celRAWqc, main="NUSE")
	dev.off()
}

.doCluster <- function(celRMA)
{
#clusters the data, and makes a boxplot

	eset <- exprs(celRMA)
	distance <- dist(t(eset),method="maximum")
	clusters <- hclust(distance)
	plot(clusters)
}
#filters the data
.doFilter <- function(celRMA)
{
	celfiles.filtered <- nsFilter(celRMA, 
															require.entrez=FALSE, 
															remove.dupEntrez=FALSE)
}


.doDE <- function(eset)
{
	samples <- eset$targets
# check the results of this
# convert into factors
	samples <- as.factor(samples)
# set up the experimental design

	design <- model.matrix(~0 + samples)
	colnames(design) <- c("minus", "plus")

# fit the linear model to the filtered expression set
	fit <- lmFit(exprs(eset), design)

# compute moderated t-statistics, moderated F-statistic, and log-odds of 
#differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
	plus_ebFit <- eBayes(fit)
# return the top 10 results for any given contrast
	ttab <- topTable(plus_ebFit, number=30000, coef=1)

	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=5))
	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=4))
	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=3))
	nrow(topTable(plus_ebFit, coef=1, number=30000, lfc=2))
# Get a list for probesets with a four fold change or more
	tTable <- topTable(plus_ebFit, number=30000, p.value=0.05)
#annotate the results
	annotation <- as.data.frame(select(chicken.db,
																		 rownames(tTable), 
																		 c("ENSEMBL","SYMBOL")))
	colnames(annotation) <- c("probeId","ensemblId","geneSymbol")
	results <- merge(annotation, tTable,by.x="probeId",by.y="row.names")

	head(results)
	write.table(results, "results.txt", sep="\t", quote=FALSE)
	return(results)
}
#a loop to carry out limma analysis if celResults exists
if(!exists("celResults"))
{
	celRAW <- .getData()
	eset<-exprs(celRAW)
	celGCRMA <- gcrma(celRAW)
	celRMA <- rma(celRAW)
	.ggboxIt(.summariseIt(log2(exprs(celRAW))),"sumRAW.pdf")
	.ggboxIt(.summariseIt(exprs(celRMA)),"sumRMA.pdf")
	.plotDensity(log2(exprs(celRAW)),"densityRAW.pdf")
	.plotDensity(log2(exprs(celRMA)),"densityRMA.pdf")
	celFilt <- .doFilter(celRMA)
	celResults <- .doDE(celFilt$eset)

}
