rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 2){

    cat(
    "Usage: 
    Rscript RNA_barplot.r inputfile outputpath
        parameters ->
            inputfile: [file -> always RNA_seQc.xls];
            outputpath: [path -> path for output]; \n")
    options("show.error.messages" = F) 
    stop()
	
}

inputfile = normalizePath(args[1])
outputpath = normalizePath(args[2])

inputdata = read.table(inputfile, header = TRUE, row.names = 1, check.name = F, comment.char = "", quote = "", sep = "\t", fill = T)
data_for_plot = t(inputdata[,2:4])

# set color
mycol = c("#F8766D", "#00BA38", "#619CFF")  
# set width and height
mc = max(nchar(colnames(data_for_plot)))
if(ncol(data_for_plot) < 6){

	pdf(paste(outputpath, "Map.region.distribution.pdf", sep = "/"), width = 8, height = 6)
	layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(1.5,1))

}else if(ncol(data_for_plot) >= 6 && ncol(data_for_plot)< 12){

	pdf(paste(outputpath, "Map.region.distribution.pdf", sep = "/"), width = 12, height = 8)
	layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(3,1))

}else{

	pdf(paste(outputpath, "Map.region.distribution.pdf", sep = "/"), width = ncol(data_for_plot), height = 8)
	layout(matrix(c(1,2), nrow = 1, byrow = T), widths = c(3,1))

}

if(mc <= 10){

	par(mar = c(6, 5.5, 2.5, 0), mgp = c(3,0.8,0), xpd = TRUE)

}else{

	par(mar = c(0.6*mc, 5.5, 2.5, 0), mgp = c(3,0.8,0), xpd = TRUE)
	
}

#plot
barplot(data_for_plot, space = 0.5,    
		beside = FALSE, 
		cex.names = 1, cex.lab = 1.3, cex.axis = 1,   
		border = NA,     
		las = 2,       
        horiz = FALSE, density = NULL, col = mycol, 
        ylab = "Fraction of mapped bases", axes = TRUE,  xaxs = "i", ylim = c(0,1))

plot.new()
par(mar = c(4,0,2.5,0))
legend("left",legend = rownames(data_for_plot), fill = mycol, bty = "n", cex = 1)
dev.off()
	

