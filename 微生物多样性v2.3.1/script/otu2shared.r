rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	"============== prepare shared file for mothur ===========
	Usage: 
	Rscript otu2shared.r otufile groupfile dissimi outputpath
		parameters ->
			otufile: [file -> always otu.tax.0.03.xls];
			groupfile: [file -> always sample.groups without header and only two column];
			dissimi: [number -> always 0.03];
		    outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
dissimi = as.numeric(args[3])
outputpath = normalizePath(args[4])

SampleInfo = read.table(groupfile, header = F, sep="\t",quote="")
otudata = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
rowname = rownames(otudata)
otudata = as.matrix(otudata[,as.character(SampleInfo[,1])])
colnames(otudata) = as.character(SampleInfo[,1])
rownames(otudata) = rowname

newdata = cbind(rep(dissimi,ncol(otudata)),colnames(otudata),rep(nrow(otudata),ncol(otudata)),t(otudata))
write.table(newdata, paste(outputpath, "sample.shared", sep = "/"), sep = "\t", quote = F, row.names = F, col.names = F)
