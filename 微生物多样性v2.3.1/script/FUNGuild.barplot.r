rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"Usage: 
	Rscript 16S.Funguild.Average.r inputfile groupfile outputpath
		parameters ->
			  inputfile: [file -> always *.taxon.Abundance.xls with row taxon and column samples];
			  groupfile: [file -> always sample.groups without header and only two column];
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

inputfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

inputdata = read.table(inputfile, header = TRUE, row.names = 1, check.name = F, 
					   comment.char = "", quote = "", sep = "\t", fill = T)

SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Group = as.character(SampleInfo[,2])
Groups = unique(as.character(SampleInfo$Group))

Data = inputdata[,which(SampleInfo[,2]==Groups[1])]
res = as.matrix(apply(Data,1,mean))
colnames(res)[ncol(res)] = Groups[1]
for (i in 2:length(Groups)){
	
	res1= as.matrix(apply(inputdata[,which(SampleInfo[,2]==Groups[i])],1,mean))
	colnames(res1)[ncol(res1)] = Groups[i]	
	res = cbind(res, res1)
	
}

avera = cbind(rownames(res),res)
colnames(avera)[1] = "Guild"
#filename = paste(outputpath, "guilds.average.xls", sep = "/")
#write.table(avera, filename, sep = "\t", row.names = F, quote = F)

# get relative.abundance
data_relative_abundance = sapply(1:ncol(res),function(x) res[,x]/sum(res[,x]))

#小于30个全画，大于等于30个按丰度排序，排序30名及之后的归为Others
data_for_plot = data_relative_abundance
colnames(data_for_plot) = colnames(res)
data_for_plot = as.matrix(data_for_plot[order(as.matrix(apply(data_for_plot,1,sum)),decreasing = T),])

if(nrow(data_for_plot) > 30){

	A = data_for_plot[30:nrow(data_for_plot),]
	A = apply(A,2,sum)
	data_for_plot = rbind(data_for_plot[1:29,],A)
	rownames(data_for_plot)[nrow(data_for_plot)] = "Others"

}

mycol = c(119,132,147,454,89,404,123,529,463,552,28,54,84,100,558,43,652,31,610,477,256,588,99,81,503,104,562,76,96,495)
mycol = colors()[rep(mycol,20)]
color = mycol[1:nrow(data_for_plot)]

# set width and height
mr = max(nchar(rownames(data_for_plot)))
if (mr < 20) mr = 20
mc = max(nchar(colnames(data_for_plot)))
nr = 0.15*ncol(data_for_plot)  

if(0.2*mr < 16){
	width = 16
}else{
	width = 0.2*mr
}

if(nr < 15){

	if(ncol(data_for_plot) < 8){
		cutpoint = c(7.5-6*nr,12*nr,7.5-6*nr)      
	}else{
		cutpoint = c(0.1,15,0.1)         
	}

}else{

	cutpoint = c(0.1,nr,0.1)         
}

#ncolu = ceiling(width*6.5/mr)        
nrowu = nrow(data_for_plot)/2     
if (nrowu < 5) nrowu = 5
height = 8 + 0.15*mc + 0.12*nrowu 
#height = 0.75*width  

pdf(paste(outputpath, "guilds.Barplot.pdf", sep="/"), width = width, height = height)
layout(rbind(c(0,1,0),c(2,2,2)), heights = c(height-0.28*nrowu, 0.28*nrowu), widths = cutpoint, respect = FALSE)    # ·ÖÃæ
par(mar = c(0.75*mc, 5.5, 2.5, 2.5), mgp = c(3,0.8,0), xpd=TRUE)
barplot(data_for_plot*100, space = 1,    
					  beside = FALSE,   
				   cex.names = 1.5,    
				   	  border = NA,     
				   	     las = 2,       
        horiz = FALSE, density = NULL, col = color, ylab = "Relative abundance (%)", axes = TRUE, cex.lab = 2.5, cex.axis = 1.2, xaxs = "i", ylim = c(0,100))
par(mar = c(2,4,3,1), xpd = TRUE)
plot.new()
legend("center", legend = rownames(data_for_plot), fill = color, ncol = 2, bty = "n", cex = 1.5)
dev.off()


