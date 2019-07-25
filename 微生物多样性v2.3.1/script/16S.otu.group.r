#!/usr/bin/Rscript
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat("Usage: 
	Rscript 16S.otu.group.r otufile groupfile outputpath
		parameters ->
			  otufile: [file -> always subsample_otu.tax.0.03.xls];
			  groupfile: [file -> always sample.groups without header and only two column];
		      outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])
setwd(outputpath)

# READ¸÷·ÖÀàË®Æ½·á¶È±í¼°·Ö×éÐÅÏ¢±í£¨È·¶¨Ñù±¾ÊýÄ¿£©
otudata = read.table(otufile, header = T, fill = T, check.names = F, sep = "\t", row.names = 1,comment.char = "")
SampleInfo = read.table(groupfile, header = F, check.names = F, sep = "\t", fill = T)

########################## Ìí¼ÓgroupÁÐ #######################

Info = unique(SampleInfo[,2])

if (length(Info) != 1 && nrow(SampleInfo) > length(Info)){
	
	res = sapply(1:length(Info), function(m) {

		group = Info[m]
		newSample = SampleInfo[which(SampleInfo[,2] %in% group), ]
		data = otudata[which(colnames(otudata) %in% newSample[,1])]

		sum = as.matrix(rowSums(data))
		colnames(sum) = paste("group",as.character(Info[m]),sep = "_")
		sum
	}
	)
	colnames(res) = Info
	rownames(res) = rownames(otudata)
	write.table(res,"otu.group.xls", sep ="\t",row.names = T, col.names = NA, quote = F)
	sample.group = data.frame(Info,Info) 
	write.table(sample.group,"sample.group", sep ="\t",row.names = F,col.names = F, quote = F)

}
#############################################################
