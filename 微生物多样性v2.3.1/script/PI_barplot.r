rm(list = ls())

args = commandArgs(trailingOnly = T)
cogfile = normalizePath(args[1])
outputpath = normalizePath(args[2])


data = read.table(cogfile, fill = T, header = T,sep = "\t",row.names = 1,
				  check.names = F,stringsAsFactors = F)

rownames(data) = gsub("\\[", "", rownames(data))
rownames(data) = gsub("\\]", ":", rownames(data))

la = LETTERS[seq(1,23)]
labels = c(la,"Y","Z")
mycol = rainbow(25)

dir = paste(outputpath,"cog_barplot",sep = "/")
system(paste("mkdir -p ", dir, sep = ""))
setwd(dir)

lapply(1:ncol(data),function(t){

	name = colnames(data)[t]
	pdf(paste(name,".pdf", sep=""),width = 14,height = 8)
	layout(matrix(c(1,2), nrow = 1, byrow = T),widths = c(1.3,1))
	
	par(mar = c(6.5,7,5.5,2) , mgp = c(5,0.8,0))
	options(scipen=999)  #不用科学计数法
	b = barplot(data[,t], ylab = "Number of Sequences", 
		col = mycol, 
		las = 2,
		main = "COG Function Classification" ,
		ylim = c(0, max(data[,t])*1.2))

	title(xlab = "Function Class",line =2 )
	text( x = b[,1], y = (par("usr")[3]-(par("usr")[4] - par("usr")[3])*0.03), labels = labels,xpd = TRUE)

	par(mar = c(0,0,0,0), xpd = TRUE)
	plot.new()
	legend("left", legend = rownames(data),cex = 0.9 ,bty = "n", y.intersp = 1.3)

	dev.off()

})