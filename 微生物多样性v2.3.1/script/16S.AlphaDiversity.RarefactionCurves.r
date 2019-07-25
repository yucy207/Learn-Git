rm(list=ls())

args=commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== Rarefaction ===========
	Usage: 
	Rscript 16S.AlphaDiversity.RarefactionCurves.r inputpath N outputpath
		parameters ->
			inputpath: [path -> contain .rarefaction files derived from mothur];
				N: [number -> smooth coefficient, always set 40];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

inputpath = normalizePath(args[1])
N = as.numeric(args[2])
outputpath = normalizePath(args[3])

rareffile = list.files(path = inputpath, pattern = ".rarefaction")

if(length(rareffile) == 0){      # check files

	stop("No Required .rarefaction files Under This Path , please check it!")

}else{

	rarefName = unlist(strsplit(rareffile, ".rarefaction.*"))
	rarefName = gsub("sample.","",rarefName)
	
	numsampled = list()
	numOTUs = list()
	
	setwd(inputpath)
	for(i in 1:length(rarefName)){

		A = read.table(rareffile[i], header=TRUE, check.name=F)
		numsampled[[i]] = A[,1]
		numOTUs[[i]] = A[,2]
		
	}
	
	xmax = max(unlist(numsampled))
	ymax = max(unlist(numOTUs))
	
	library(RColorBrewer)
	setwd(outputpath)
	
	#ColorsInR = grDevices::colors()[grep("gr(a|e)y|light|skyblue|oldlace|linen|powderblue|rosybrown1|sea|papayawhip|wheat|pink|moccasin|rosybrown|cornsilk|azure|bisque|burlywood|white|lavender|gainsboro|honeydew|lemonchiffon|tan|ivory|mistyrose|snow|thistle",grDevices::colors(),invert=T)]
	#color = sample(ColorsInR, length(rarefName))      # randomly selected colors
	mycol = c(119,132,147,454,553,89,123,463,552,28,117,54,84,100,43,652,31,610,477,529,588,614,99,81,491,503,104,562,404,76,96,495,372,108,75,622,631,471,653)
	color = colors()[rep(mycol,20)][1:length(rarefName)]
	
	pdf("Samples_Rarefaction_Curves.pdf", width=12, height=11)
	par(mar = c(3.5,3.5,4.5,1), mgp = c(2,0.5,0))
	
	cex = ifelse(17/length(rarefName) > 1.4, 1.4, 17/length(rarefName))
	
	plot(spline(numsampled[[1]], numOTUs[[1]], method = "hyman", n = N), lwd = 4, col = color[1], type = "l", xlim = c(1,xmax+xmax/5), ylim = c(1,ymax+ymax/3), xlab = "Number of Reads Sampled", ylab = "Number of OTUs", main = "Samples Rarefaction Curves", cex.lab = 1.3, cex.main = 1.6)
	# text(x = max(numsampled[[1]]), y = max(numOTUs[[1]]), cex = cex, labels=rarefName[1], pos=4, col=color[1])
	
	if( length(rarefName) > 1 ){
	
		for(j in 2:length(rarefName)){
	
			lines(spline(numsampled[[j]], numOTUs[[j]], method="hyman", n = N), lwd = 4, col = color[j])
			# text(x = max(numsampled[[j]]), y = max(numOTUs[[j]]), cex = cex, labels = rarefName[j], pos = 4, col = color[j])
		
		}
	
	}
	
	legend("topright", legend = rarefName, col = color, lty = 1 , cex = 0.9, bty = "n", ncol = ceiling(length(rarefName)/12))
	dev.off()

}
