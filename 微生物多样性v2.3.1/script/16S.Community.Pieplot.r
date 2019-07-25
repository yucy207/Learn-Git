#!/usr/bin/Rscript
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== Pieplot ===========
	Usage: 
	Rscript 16S.Community.Pieplot.r taxonfile groupfile outputpath
		parameters ->
			  taxonfile: [file -> always *.taxon.Abundance.xls];
			  groupfile: [file -> always sample.groups without header and only two column];
		      outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

taxonfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])


# READ������ˮƽ��ȱ�������Ϣ��ȷ��������Ŀ��
taxondata = read.table(taxonfile, header = T, fill = T, check.names = F, sep = "\t", row.names = 1,comment.char = "")
SampleInfo = read.table(groupfile, header = F, check.names = F, sep = "\t", fill = T)

#ȡ��������ˮƽ�а�������������
rownames = rownames(taxondata)
colnames = as.character(SampleInfo[,1])
taxondata = as.matrix(taxondata[,as.character(SampleInfo[,1])]) ###ֻ��һ��ʱ
rownames(taxondata) = rownames                                  ###ֻ��һ��ʱ
colnames(taxondata) = colnames                                  ###ֻ��һ��ʱ

#�ҵ�����all�ַ�������
index = which(rownames(taxondata) == "All")
if(length(index) > 0) {
	taxondata = as.matrix(taxondata[-index,])
	colnames(taxondata) = colnames                              ###ֻ��һ��ʱ
	if (nrow(taxondata) == 1) rownames(taxondata) = rownames[which(rownames != "All")]
}

########################## ���group�� #######################
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("Sample", "Group")
Info = unique(SampleInfo[[2]])
if (length(Info) != 1 && nrow(SampleInfo) > length(Info)){
	for(m in 1:(length(Info))){

		Sample = Info[m]
		newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
		data = as.matrix(taxondata[,as.character(newSample[,1])])
		#taxondata = cbind(taxondata, apply(data, 1, sum))
		#colnames(taxondata)[ncol(taxondata)] = as.character(Info[m])

		sum = as.matrix(rowSums(data))
		colnames(sum) = paste("group",as.character(Info[m]),sep = "_")
		taxondata = cbind(taxondata,sum)
	}	
}
#############################################################

names = strsplit(strsplit(taxonfile,"/")[[1]][length(strsplit(taxonfile,"/")[[1]])],"[.]")[[1]][1]
# dir = paste(outputpath, names, sep = "/")
# system(paste("mkdir -p ", dir, sep = ""))
# setwd(dir)

# plot��ͼ
mycol = c(119,132,147,454,404,123,529,463,552,28,54,84,100,558,89,43,652,31,610,477,256,588,99,81,503,104,562,76,96,495)
mycol = colors()[rep(mycol,20)]

for (i in 1:ncol(taxondata)) {
	
	A = as.matrix(taxondata[,i])
	#rownames(A) = rownames(taxondata)
	#colnames(A) = clonames(taxondata)[i]
	#A = as.matrix(A[which(A[,1]!=0),])
	
	# ��Է��
	B = round(A/sum(A)*100,2)
	
	loc = which(B[,1] < 1)
	if(length(loc) > 0){	
		C = rbind(as.matrix(B[-loc,]), as.matrix(sum(B[loc,])))
		rownames(C)[nrow(C)] = "Others"
		C = as.matrix(C[order(C[,1],decreasing = T),])
	}else {
		C = as.matrix(B[order(B[,1],decreasing = T),])
		if (nrow(C) == 1) rownames(C)[nrow(C)] = rownames(taxondata)
	}
	
	#plot pie
	color = mycol[1:nrow(C)] 
	lbls = rownames(C)
	lbls2 = paste(lbls, "(", C, "%", ")", sep = "")
	mr = max(nchar(lbls2))
	
	pdf(paste(outputpath, paste(colnames(taxondata)[i],names,"pieplot.pdf", sep = "."), sep  ="/"), width = 9, height = 7)
	
	#���ñ�ͼ��ͼ�е�λ�ã��£����ϣ��ң� 
	par(mar = c(0.5, 8, 0.5, 8))
	pie(C, labels = lbls2, col = color, cex = 8/mr+0.3)
	dev.off()
}
