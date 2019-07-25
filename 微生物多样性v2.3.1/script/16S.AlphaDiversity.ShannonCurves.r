rm(list=ls())

args=commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== shannon ===========
	Usage: 
	Rscript 16S.AlphaDiversity.ShannonCurves.r inputpath N outputpath
		parameters ->
			inputpath: [path -> contain .r_shannon files derived from mothur];
				N: [number -> smooth coefficient, always set 40];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

inputpath = normalizePath(args[1])
N = as.numeric(args[2])
outputpath = normalizePath(args[3])

shannonfile = list.files(path = inputpath, pattern = ".r_shannon")

if(length(shannonfile) == 0){      # check files

	stop("No .r_shannon files Under This Path , pleae check it!")

}else{

	shannonName = unlist(strsplit(shannonfile, ".r_shannon.*"))
	shannonName = gsub("sample.","",shannonName)
	
	numsampled = list()
	shannnonIndex = list()
	
	setwd(inputpath)
	for(i in 1:length(shannonName)){

		A = read.table(shannonfile[i], header = TRUE, check.name = F)
		A[1,2:4] = rep(0, 3)
		numsampled[[i]] = A[,1]
		shannnonIndex[[i]] = A[2]
		
	}
	
	xmax = max(unlist(numsampled))
	ymax = max(unlist(shannnonIndex))+1.5
	
	library(RColorBrewer)
	setwd(outputpath)
	
	mycol = c(119,132,147,454,553,89,123,463,552,28,117,54,84,100,43,652,31,610,477,529,588,614,99,81,491,503,104,562,404,76,96,495,372,108,75,622,631,471,653)
	color = colors()[rep(mycol,20)][1:length(shannonName)]
	
	pdf("Samples_Shannon_Curves.pdf", width = 11, height = 11)
	par(mar = c(3.5,3.5,4.5,1), mgp = c(2,0.5,0))
	
	cex = ifelse(17/length(shannonfile) > 1.4, 1.4, 17/length(shannonfile))
	
	plot(spline(numsampled[[1]], shannnonIndex[[1]][[1]], n = N), lwd = 2, col = color[1], type = "l", xlim = c(0, xmax+xmax/5), ylim = c(0,ymax), xlab = "Number of Reads Sampled", ylab = "Shannon Index", main = "Samples Shannon Curves", cex.lab = 1.3,cex.main = 1.6)
	# text(x = max(numsampled[[1]]), y = max(shannnonIndex[[1]]), cex = cex, labels = shannonName[1], pos = 4, col = color[1])
	
	if( length(shannonName) > 1 ){
	
		for(j in 2:length(shannonName)){
	
			lines(spline(numsampled[[j]], shannnonIndex[[j]][[1]], n = N), lwd = 2, col = color[j])
			# text(x = max(numsampled[[j]]), y = max(shannnonIndex[[j]]), cex = cex, labels = shannonName[j], pos = 4, col = color[j])
		
		}
		
	}
	legend("topright", legend = shannonName, col = color, lty = 1 , cex = 0.9, bty = "n", ncol = ceiling(length(shannonName)/12))
	dev.off()

}
