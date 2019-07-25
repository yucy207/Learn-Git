rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"Usage: 
	Rscript 16S.Community.Barplot.r inputfile outputpath
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


inputdata = read.table(inputfile, header = T, fill = T, check.names = F, sep = "\t", row.names = 1,comment.char = "") 
SampleInfo = read.table(groupfile, header = F, check.names = F, sep = "\t", fill = T)
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("Sample", "Group")
Info = unique(SampleInfo[[2]])

B = inputdata
#分组
for(m in 1:(length(Info))){

    Sample = Info[m]
    newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
    newinputdata = as.matrix(inputdata[,as.character(newSample[,1])])
    datamean = rowSums(newinputdata)/ncol(newinputdata)
    Data = as.matrix(datamean)
    colnames(Data)= Info[m]
    rownames(Data)=rownames(inputdata)

    B = cbind(B,Data)
}

z = length(Info)
n = ncol(B)
s = n - z 
s = s + 1

Data = as.matrix(B[s:n])
colnames(Data) = colnames(B)[s:n]

# prepare data
rownames(Data) = gsub("\"","",rownames(Data))     # È¥³ýÃû×ÖÖÐµÄË«ÒýºÅ,È«¾ÖÌæ»»
colnames(Data) = gsub("\"","",colnames(Data))
data1 = as.matrix(Data)
rownames(data1) = as.character(strsplit(rownames(data1),split="{.*}",perl=T))

loc = which(rownames(data1) == "All")       # È¥³ý All ÐÐ
if(length(loc) > 0) data1 = as.matrix(data1[-loc, ])

if( ncol(data1) == 1){

	data1 = as.matrix(data1[data1 != 0, ])

}else{

	data1 = as.matrix(data1[rowSums(data1) != 0, ])

}

if( !is.null(nrow(data1)) ){     # Ö»ÓÐÒ»ÀàµÄ,²»»­

	# get relative.abundance
	data_relative_abundance = sapply(1:ncol(data1),function(x) data1[,x]/sum(data1[,x]))

	# ÖÖÀàÉÙÓÚ30¸öµÄÈ«»­£¬¶àÓÚ30¸öµÄ£¬½«Ïà¶Ô·á¶ÈµÍµÄ½øÐÐºÏ²¢ÎªOthers½øÐÐÕ¹Ê¾£¡
	if(nrow(data_relative_abundance) < 30){

		data_for_plot = data_relative_abundance
		colnames(data_for_plot) = colnames(Data)
		data_for_plot = as.matrix(data_for_plot[order(as.matrix(data_for_plot[,1]),decreasing = T),])

	}else{

		# ·ÇÔÚËùÓÐÑù±¾ÖÐÏà¶Ô·á¶È¶¼´óÓÚ1%µÄ¾Í¹éÎªothers
		index = which(apply(data_relative_abundance > 0.01, 1, any)==T)
		A = as.matrix(data_relative_abundance[index,])
		A = cbind(A,as.matrix(apply(A,1,sum)))
		A = as.matrix(A[order(A[,ncol(A)],decreasing = T),])
		A = A[,-ncol(A)]
                A = as.matrix(A)
		data_for_plot = rbind(A,apply(as.matrix(data_relative_abundance[-index,]),2,sum))
		rownames(data_for_plot)[nrow(data_for_plot)] = "Others"

	}

	if(nrow(data_for_plot) > 2){
		index = grep("No_Rank",rownames(data_for_plot))
		if(length(index) > 0) {
			data_for_plot = rbind(as.matrix(data_for_plot[-index,]),data_for_plot[index,])
			rownames(data_for_plot)[nrow(data_for_plot)] = "No_Rank"
		}
	}
	
	data_for_plot = as.matrix(data_for_plot)
	colnames(data_for_plot) = colnames(Data)
	
	# set color
	mycol = c(119,132,147,454,89,404,123,529,463,552,28,54,84,100,558,43,652,31,610,477,256,588,99,81,503,104,562,76,96,495)
	mycol = colors()[rep(mycol,20)]
	color = mycol[1:nrow(data_for_plot)] 

	# set width and height
	mr = max(nchar(rownames(data_for_plot)))
	if (mr < 20) mr = 20
	mc = max(nchar(colnames(data_for_plot)))

	nr = 0.15*ncol(data_for_plot)    # Ö»¸ù¾ÝÑù±¾ÊýÉè¶¨Í¼ÐÎ¿í¶È

	if(nr < 15){

		width = 15
		if(ncol(data_for_plot) < 8){
			cutpoint = c(7.5-6*nr,12*nr,7.5-6*nr)       # Ñù±¾Á¿ÉÙÓÚ8¸ö£¬ÔòÑ¹Ëõ»æÍ¼ÇøÓò£¬À©´óÁ½²àÇøÓò
		}else{
			cutpoint = c(0.1,width,0.1)          # ÈôÑù±¾Á¿½Ï´ó£¬Ôò¿ÉÒÔ°ÑÁ½²à¿Õ°×ÇøÓòÑ¹Ëõ
		}

	}else{

		width = nr
		cutpoint = c(0.1,width,0.1)          # ÈôÑù±¾Á¿½Ï´ó£¬Ôò¿ÉÒÔ°ÑÁ½²à¿Õ°×ÇøÓòÑ¹Ëõ
	}

	ncolu = ceiling(width*6.5/mr)         # ¸ù¾ÝÍ¼ÐÎ¿í¶È£¬ÏÞÖÆ±êÇ©ÁÐÊý
	nrowu = nrow(data_for_plot)/ncolu      # ±êÇ©ÐÐÊý
	if (nrowu < 5) nrowu = 5
	height = 8 + 0.15*mc + 0.2*nrowu       # ¸ù¾Ý bar labels×Ö·û³¤¶ÈºÍ±êÇ©ÐÐÊýÏÞ¶¨¸ß¶È

	name = unlist(strsplit(inputfile,"/"))
	pdf(paste(outputpath, paste(strsplit(name[length(name)], "\\.\\w+\\.xls"), "Community.Barplot.group.pdf", sep = "."), sep="/"), width = width, height = height)
	layout(rbind(c(0,1,0),c(2,2,2)), heights = c(height-0.35*nrowu, 0.35*nrowu), widths = cutpoint, respect = FALSE)    # ·ÖÃæ
	par(mar = c(0.75*mc, 5.5, 2.5, 2.5), mgp = c(3,0.8,0), xpd=TRUE)
	barplot(data_for_plot*100, space = 1,     # ÌõÐÎ¼ä¾à
				  beside = FALSE,   # beside ¶Ñ»ýÌõÐÎÍ¼
			   cex.names = 1.5,     # bar labels ×ÖÌå´óÐ¡
				  border = NA,      # ÎÞÌõÐÎ±ß¿ò
					 las = 2,       # las±êÇ©ÊúÏò
		 horiz = FALSE, density = NULL, col = color, ylab = "Relative abundance (%)", axes = TRUE, cex.lab = 1.8, cex.axis = 1.2, xaxs = "i", ylim = c(0,100))
	par(mar = c(1,2,3,1), xpd = TRUE)
	plot.new()
	legend("center", legend = rownames(data_for_plot), fill = color, ncol = ncolu, bty = "n", cex = 1.8)
	dev.off()
}
	

