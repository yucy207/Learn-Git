rm(list = ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== do betaNTI ===========
	Usage: 
	Rscript format.r otufile groupfile fasta outputpath
            parameters ->
                otufile: [file -> always always subsample_otu.tax.0.03.xls];
                groupfile: [file -> always sample.groups];
                outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])
setwd(outputpath)

SampleInfo = read.table(groupfile, header = F, sep = "\t",quote="")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Groups = unique(SampleInfo$Group)

otudata = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote = "", fill = T)

dataforplot = data.frame()

for(i in 1:length(Groups)){

	newSample = as.character(SampleInfo$Sample[which(SampleInfo$Group %in% Groups[i])])
	otu_data = as.matrix(otudata[,newSample])

	res = data.frame()
	for(m in 1:ncol(otu_data)){

		df = data.frame(sample = rep(colnames(otu_data)[m],times = nrow(otu_data)), value = otu_data[,m], OTU = rownames(otu_data))
		df1 = df[-which(df$value == 0),]
		res = rbind(res,df1)

	}

	write.table(res,paste(Groups[i],"format.xls",sep ="."),sep = "\t",row.names = F,col.names = F,quote=F)
	
}
