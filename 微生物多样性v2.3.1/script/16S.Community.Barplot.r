rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"Usage: 
	Rscript 16S.Community.Barplot.r inputfile outputpath
		parameters ->
			  inputfile: [file -> always *.taxon.Abundance.xls with row taxon and column samples];
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

inputfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

# inputfile = normalizePath("/home/zhengxq/workdir/meta_genomics/16C1025A/ITS2/report/Community/Community_Structure/class/class.taxon.Abundance.xls")
# outputpath = normalizePath("/home/zhengxq/workdir/meta_genomics/16C1025A/ITS2/report/Community/Barplot")

SampleInfo = read.table(groupfile, header = F, check.names = F, sep = "\t", fill = T)
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("Sample", "Group")
sample = rownames(SampleInfo)

inputdata = read.table(inputfile, header = TRUE, row.names = 1, check.name = F, comment.char = "", quote = "", sep = "\t", fill = T)    # ע�� quote="" 
# ������ڷ���ע����Ϣ,��ȥ��
loc = grep("OTUsize|Abundance", colnames(inputdata))

if(length(loc) > 0){

	Data = as.matrix(inputdata[,1:(loc-1)])
	colnames(Data) = colnames(inputdata)[1:(loc-1)]

}

Data = Data[,sample]
rownames(Data) = rownames(inputdata)

# prepare data
rownames(Data) = gsub("\"","",rownames(Data))     # ȥ�������е�˫����,ȫ���滻
colnames(Data) = gsub("\"","",colnames(Data))
data1 = as.matrix(Data)
rownames(data1) = as.character(strsplit(rownames(data1),split="{.*}",perl=T))

loc = which(rownames(data1) == "All")       # ȥ�� All ��
if(length(loc) > 0) data1 = as.matrix(data1[-loc, ])

if( ncol(data1) == 1){

	data1 = as.matrix(data1[data1 != 0, ])

}else{

	data1 = as.matrix(data1[rowSums(data1) != 0, ])

}

if( !is.null(nrow(data1)) ){     # ֻ��һ���,����

	# get relative.abundance
	data_relative_abundance = sapply(1:ncol(data1),function(x) data1[,x]/sum(data1[,x]))
	
	# ��������30����ȫ��������30���ģ�����Է�ȵ͵Ľ��кϲ�ΪOthers����չʾ��
	if(nrow(data_relative_abundance) < 30){

		data_for_plot = data_relative_abundance
		colnames(data_for_plot) = colnames(Data)
		data_for_plot = as.matrix(data_for_plot[order(as.matrix(data_for_plot[,1]),decreasing = T),])

	}else{

		# ����������������Է�ȶ�����1%�ľ͹�Ϊothers
		index = which(apply(data_relative_abundance > 0.01, 1, any)==T)
		A = as.matrix(data_relative_abundance[index,])
		A = cbind(A,as.matrix(apply(A,1,sum)))
		A = as.matrix(A[order(A[,ncol(A)],decreasing = T),])
		A = A[,-ncol(A)]
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

	nr = 0.15*ncol(data_for_plot)    # ֻ�����������趨ͼ�ο��

	if(nr < 15){

		width = 15
		if(ncol(data_for_plot) < 8){
			cutpoint = c(7.5-6*nr,12*nr,7.5-6*nr)       # ����������8������ѹ����ͼ����������������
		}else{
			cutpoint = c(0.1,width,0.1)          # ���������ϴ�����԰�����հ�����ѹ��
		}

	}else{

		width = nr
		cutpoint = c(0.1,width,0.1)          # ���������ϴ�����԰�����հ�����ѹ��
	}

	ncolu = ceiling(width*6.5/mr)         # ����ͼ�ο�ȣ����Ʊ�ǩ����
	nrowu = nrow(data_for_plot)/ncolu      # ��ǩ����
	if (nrowu < 5) nrowu = 5
	height = 8 + 0.15*mc + 0.2*nrowu       # ���� bar labels�ַ����Ⱥͱ�ǩ�����޶��߶�

	name = unlist(strsplit(inputfile,"/"))
	pdf(paste(outputpath, paste(strsplit(name[length(name)], "\\.\\w+\\.xls"), "Community.Barplot.pdf", sep = "."), sep="/"), width = width, height = height)
	layout(rbind(c(0,1,0),c(2,2,2)), heights = c(height-0.35*nrowu, 0.35*nrowu), widths = cutpoint, respect = FALSE)    # ����
	par(mar = c(0.75*mc, 5.5, 2.5, 2.5), mgp = c(3,0.8,0), xpd=TRUE)
	barplot(data_for_plot*100, space = 1,     # ���μ��
					  beside = FALSE,   # beside �ѻ�����ͼ
				   cex.names = 1.5,     # bar labels �����С
				   	  border = NA,      # �����α߿�
				   	     las = 2,       # las��ǩ����
        horiz = FALSE, density = NULL, col = color, ylab = "Relative abundance (%)", axes = TRUE, cex.lab = 1.8, cex.axis = 1.2, xaxs = "i", ylim = c(0,100))
	par(mar = c(1,2,3,1), xpd = TRUE)
	plot.new()
	legend("center", legend = rownames(data_for_plot), fill = color, ncol = ncolu, bty = "n", cex = 1.8)
	dev.off()
	
}
