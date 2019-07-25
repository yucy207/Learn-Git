# prepare lefse result for plot 
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== prepare lefse result for plot  ===========
	Usage: 
	Rscript lefse4plot.r resfile TopN outputpath
		parameters ->
			resfile: [file -> always taxon.stat.res.xls derived from run_lefse.py];
			    TopN: [number -> N most significant items to plot, always 50];
		     outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

resfile = normalizePath(args[1])
TopN = as.numeric(args[2])
outputpath = normalizePath(args[3])
setwd(outputpath)

res = read.table(resfile, header = F, sep = "\t", na.strings = "")

clean_res = res[which(as.matrix(res[,5]) != "-"), ]
resout = clean_res[which(!is.na(as.matrix(clean_res[,4]))), c(1,3,4,5)]
colnames(resout) = c("Taxonomy", "Group", "LDA", "Pvalue")
write.table(resout, "taxon.lefse.result.xls", row.names = F, sep = "\t", quote = F, na = "")

groups = unique(res[,3])[which(!is.na(unique(res[,3])))]

if (nrow(clean_res) > 280){

	res4plot = res[order(res[,5]),]
	
	for (i in 1:length(groups)){
	
		if (i == 1){
		
			if( length(which(res4plot[,3] == groups[i])) < TopN ){
			
				dataforplot = res4plot[which(res4plot[,3] == groups[i]),]
				
			}else{
			
				dataforplot = res4plot[which(res4plot[,3] == groups[i]),][1:TopN,]
				
			}
		
		}else{
		
			if( length(which(res4plot[,3] == groups[i])) < TopN ){
			
				dataforplot = rbind(dataforplot, res4plot[which(res4plot[,3] == groups[i]),])
				
			}else{
			
				dataforplot = rbind(dataforplot, res4plot[which(res4plot[,3] == groups[i]),][1:TopN,])
				
			}
		
		}
		
	}
	
	index = which((res4plot[,1] %in% dataforplot[,1]) == "FALSE")
	res4plot[index, 5] = "-"
	res4plot[which(res4plot[,5] == "-"), c(3, 4)] = NA

}else{

	res4plot = res

}

write.table(res4plot, "taxon.stat.res4plot.xls", row.names = F, col.names = F, sep = "\t", quote = F, na = "")
