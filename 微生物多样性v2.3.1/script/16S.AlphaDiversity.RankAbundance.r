rm(list=ls())

args=commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== plot OTU Rank Abundane ===========
	Usage: 
	Rscript 16S.AlphaDiversity.RankAbundance.r otufile outputpath
		parameters ->
			otufile: [file -> always otu.tax.0.03.xls];
		      groupfile: [file -> always sample.groups without header and only two column];
		     outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

# ====================== prepare data =======================
SampleInfo = read.table(groupfile, row.names = 1, header = F, sep = "\t",quote="")

otudata = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
names = rownames(otudata)

otudata = as.matrix(otudata[,rownames(SampleInfo)])
rownames(otudata) = names
colnames(otudata) = rownames(SampleInfo)

# sort and rank Relative Abundance
A = apply(otudata, 2, function(x) sort(x, decreasing = T))
RA = apply(A, 2, function(d) d/sum(as.numeric(d)))*100

# get ranks that Relative Abundance >0 in all samples for plot
index = which(apply(RA > 0, 1, any)==T)
data_for_plot = as.matrix(RA[index,])
colnames(data_for_plot) = rownames(SampleInfo)

# prepare set
mycol = c(119,132,147,454,553,89,404,123,529,463,552,28,117,653,54,84,100,558,43,652,31,610,477,588,614,99,81,491,503,104,562,76,96,495,372,108,75,622,631,471)
color = colors()[rep(mycol,20)][1:ncol(data_for_plot)]
ylim = -3:2

# x is rank; y is Relative Abundance
pdf(paste(outputpath, "Rank.Abundance.pdf", sep = "/"))
par(mar = c(2.5, 2.5, 3, 1), mgp = c(1.5, 0.4, 0))
plot(c(1:nrow(data_for_plot)), log10(data_for_plot[,1]), type = "l", xlab="OTU Rank", ylab="Relative Abundance (%)", main="Rank-abundance distribution curve", cex.main = 1.2, cex.axis=0.7, yaxt='n', ylim = c(-3,2), col = color[1])

if (ncol(data_for_plot) > 1){

	for(i in 2:ncol(data_for_plot)){
	
		lines(c(1:nrow(data_for_plot)), log10(data_for_plot[,i]), col = color[i])
		
	}
	
}

axis(side = 2, at = ylim, labels = as.character(10^ylim), cex.axis = 0.7) 
legend("topright", legend = colnames(data_for_plot), col = color, lty = 1 , cex = 0.7, bty = "n", ncol = ceiling(ncol(data_for_plot)/26))
dev.off()
