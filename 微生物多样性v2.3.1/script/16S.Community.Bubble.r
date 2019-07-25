rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	"Usage: 
	Rscript 16S.Community.Bubble.r inputfile outputpath
		parameters ->
			  inputfile: [file -> always *.taxon.Abundance.xls with row taxon and column samples];
			  groupfile: [file -> always sample.groups without header and only two column];
			  N: [number -> always N most abundant taxa to plot Ternary];
			  outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

library(ggplot2)

inputfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
N = as.numeric(args[3])
outputpath = normalizePath(args[4])

setwd(outputpath)

# ====================== prepare data =======================
SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Input_file = read.table(inputfile, header = TRUE, row.names = 1, check.name = F, comment.char = "", quote = "", sep = "\t", fill = T)

###只有一个样本时需要重新赋值行名，列名
input_data = as.matrix(Input_file[,rownames(SampleInfo)])
colnames(input_data) = SampleInfo[,1]
rownames(input_data) = rownames(Input_file)

loc = which(rownames(input_data)=="All")

###只有一个物种时不画图
if(length(loc) > 0 & nrow(input_data) == 2)  stop(" taxon number less than 2 ! Please Cheack It! ***\n")

if(length(loc) > 0 & nrow(input_data) > 2) {
	input_data = as.matrix(input_data[-loc,])
	colnames(input_data) = SampleInfo[,1]
	rownames(input_data) = rownames(Input_file)[-loc]
}

input_data = as.matrix(input_data[order(apply(input_data,1,sum),decreasing = T),])
colnames(input_data) = SampleInfo[,1]

# ====================== prepare data =======================

# most abundant taxa
if(nrow(input_data) < N){
	DataForPlot = as.matrix(input_data)
}else{
	DataForPlot = as.matrix(input_data[1:N,])
}
colnames(DataForPlot) = colnames(input_data)

# relative abundance
re = sapply(1:ncol(DataForPlot), function(x) DataForPlot[,x]/sum(DataForPlot[,x]))
rownames(re) = rownames(DataForPlot)
colnames(re) = colnames(DataForPlot)

# prepare data
data = stack(as.data.frame(re))
data$ylabel = rep(rownames(DataForPlot), times = ncol(DataForPlot))
data$groups = rep(as.character(SampleInfo[,2]), each = nrow(DataForPlot))
colnames(data) <- c("Abundance", "xlabel", "ylabel", "Groups")

# sort
data$xlabel = factor(data$xlabel, levels = unique(data$xlabel))

# plot
mycol = c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol = colors()[mycol[1:length(unique(data$Groups))]]

names = strsplit(unlist(strsplit(inputfile,"/"))[length(unlist(strsplit(inputfile,"/")))],".",fixed=T)[[1]][1]

pdf(paste(outputpath, paste(names, "taxon.Community.Bubble.pdf", sep = "."), sep="/"), width = 6 + 0.6*ncol(re), height = 4 + 0.1*nrow(re))

ggplot(data, aes(x = xlabel, y = ylabel, size = Abundance, colour = Groups)) + 
#guides(colour = guide_legend()) +		                      # combine size and color legend
geom_point() + theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +    # xÖá±êÇ©ÇãÐ±
theme(panel.border = element_blank()) +                       # È¥³ý×îÍâ²ã±ß¿ò
scale_size(range = c(0,6)) + labs(x = "Sample", y = "") +
scale_color_manual(values = mycol)		                      # custom colors

dev.off()


########################## Ìí¼ÓgroupÁÐ #######################
Info = unique(SampleInfo[,2])
###样本无重复时不画分组的图
if ( length(Info) == 1 || nrow(SampleInfo) == length(Info) ) stop(" sample repeat number less than 2 ! Please Cheack It! ***\n")

res = lapply(1:length(Info), function(t){
	Sample = Info[t]
	newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
	data = as.matrix(input_data[,as.character(newSample[,1])])
	sum = as.matrix(rowSums(data))
	colnames(sum) = Sample
	sum})

new_data <- as.matrix(do.call(cbind, res))

# most abundant taxa
if(nrow(new_data) < N){
	DataForPlot = as.matrix(new_data)
	rownames(DataForPlot) = rownames(new_data)
}else{
	DataForPlot = as.matrix(new_data[1:N,])
	rownames(DataForPlot) = rownames(new_data)[1:N]
}
colnames(DataForPlot) = colnames(new_data)

# relative abundance
re = sapply(1:ncol(DataForPlot), function(x) DataForPlot[,x]/sum(DataForPlot[,x]))
rownames(re) = rownames(DataForPlot)
colnames(re) = colnames(DataForPlot)

# prepare data
data = stack(as.data.frame(re))
data$ylabel = rep(rownames(DataForPlot), times = ncol(DataForPlot))
data$groups = rep(colnames(DataForPlot), each = nrow(DataForPlot))
colnames(data) <- c("Abundance", "xlabel", "ylabel", "Groups")

# sort
data$xlabel = factor(data$xlabel, levels = unique(data$xlabel))

# plot
mycol = c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol = colors()[mycol[1:length(unique(data$Groups))]]

names = strsplit(unlist(strsplit(inputfile,"/"))[length(unlist(strsplit(inputfile,"/")))],".",fixed=T)[[1]][1]
pdf(paste(outputpath, paste(names, "taxon.Community.Bubble.group.pdf", sep = "."), sep="/"), width = 6 + 0.6*ncol(re), height = 4 + 0.1*nrow(re))

ggplot(data, aes(x = xlabel, y = ylabel, size = Abundance, colour = Groups)) + 
#guides(colour = guide_legend()) +		                      # combine size and color legend
geom_point() + theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +    # xÖá±êÇ©ÇãÐ±
theme(panel.border = element_blank()) +                       # È¥³ý×îÍâ²ã±ß¿ò
scale_size(range = c(0,6)) + labs(x = "Groups", y = "") +
scale_color_manual(values = mycol)		                      # custom colors

dev.off()

