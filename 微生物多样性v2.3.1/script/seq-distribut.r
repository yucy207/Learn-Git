rm(list=ls())

args = commandArgs(trailingOnly=T) 

if(length(args) != 2){

	cat(
	"Usage: 
	Rscript seq-distribut.r inputfile outputpath
		parameters ->
		    inputfile: [file -> always *length.distribution.xls derived from seq-distribut.pl];
		   outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

inputfile = normalizePath(args[1])
outputpath = normalizePath(args[2])

name = unlist(strsplit(inputfile,"/|.xls"))[length(unlist(strsplit(inputfile,"/|.xls")))]

setwd(outputpath)
stat = read.table(inputfile, header = FALSE, sep = "\t", row.names = 1)

# 0.01%
index = which(stat[,1] == max(stat[,1]))
rig = ifelse(index + 10 < nrow(stat), index + 10, nrow(stat))
lef = ifelse(index - 10 > 0, index - 10, 1)

data = stat[c(lef:rig),1]
names(data) = rownames(stat)[c(lef:rig)]

pdf(paste(name, ".pdf", sep = ""), width = 12, height = 9)

if (max(nchar(names(data))) < 4){

	par(mar = c(4.5, 4.5, 3, 1), mgp = c(3, 0.8, 0))
	bp = barplot(data, col = "#4F81BD", border = FALSE, main = "Sequence length distribution", xlab = "Length", ylab = "Number", ylim = c(0, max(data)*1.2), cex.lab = 1.1)

	
}else{
	
	par(mar = c(6, 6, 3, 1), mgp = c(4.8, 0.8, 0))
	bp = barplot(data, col = "#4F81BD", border = FALSE, axisnames = FALSE, main = "Sequence length distribution", xlab = "Length", ylab = "Number", ylim = c(0, max(data)*1.2), cex.lab = 1.1)
	text(x = bp, y = -0.01*max(data), srt = 70, adj = 1, labels = names(data), xpd = T, cex = 1)

}

dev.off()
