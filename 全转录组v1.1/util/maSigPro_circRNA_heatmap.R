rm(list=ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
    "Usage: 
    Rscript maSigPro_heatmap.R diff_file outputpath
        parameters ->
            diff_file: [file -> always all.diff.xls];
            outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()
	
}

library(gplots)

diff_file  = normalizePath(args[1])
outputpath = normalizePath(args[2])

diffdata = read.table(diff_file, header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "", quote = "", fill = T)

#heatmap
myheatcol    = colorpanel(75, 'green','black','red')
index        = which(colnames(diffdata) == "beta0")
data         = data.matrix(diffdata[diffdata$type != "Not DEG", 3:index-1])
heatmap_file = paste(outputpath, "/heatmap.pdf", sep="")
pdf(heatmap_file , width = 12, height = 12)
heatmap.2( data,
           dendrogram   = 'row', 
           Colv         = FALSE,
           col          = myheatcol, 
           scale        = "row", 
           margins      = c(5,10),
           density.info = "none",
           trace        = "none",
           key          = TRUE,
           keysize      = 0.6,
           cexRow       = 0.01,
           cexCol       = 0.8,
           srtCol       = 90
)
dev.off()

#Top50
res = diffdata[diffdata$type != 'Not DEG', ]
if( nrow(res) < 50 ){

    data = data.matrix(res[ , 3:index-1])          # 小于50,全画

}else{

    res  = res[order(res$"p-value"), ]
    data = data.matrix(res[1:50, 3:index-1])       #  相对丰度最高的50个

}

heatmap_file = paste(outputpath, "/Top50.heatmap.pdf", sep="")
pdf(heatmap_file , width = 12, height = 12)
heatmap.2( data,
           dendrogram   = 'both', 
           #Colv         = FALSE,
           col          = myheatcol, 
           scale        = "row", 
           margins      = c(5,12),
           density.info = "none",
           trace        = "none",
           key          = TRUE,
           keysize      = 0.6,
           cexRow       = 0.8,
           cexCol       = 0.8,
           srtCol       = 90
)

dev.off()
