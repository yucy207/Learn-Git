rm(list=ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== Venn ===========
	Usage: 
	Rscript 16S.Community.Venn.r taxonfile groupfile outputpath
		parameters ->
			taxonfile: [file -> always *.taxon.Abundance.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(plyr)
library(VennDiagram)

taxonfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

setwd(outputpath)

# Ñù±¾ÐÅÏ¢
SampleInfo = read.table(groupfile, header = F, sep="\t")
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("Sample", "Group")
Info = unique(SampleInfo[[2]])

# ÎïÖÖÊý¾Ý
taxondata = read.table(taxonfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
count = taxondata[, rownames(SampleInfo)] 
# prepare data
rownames(count) = gsub("\"","",rownames(count))     
colnames(count) = gsub("\"","",colnames(count))
rownames(count) = as.character(strsplit(rownames(count),split="{.*}",perl=T))
# È¥³ý All ºÍ No_Rank ÐÐ
count = count[setdiff(rownames(count), c("All", "No_Rank")),]

venn_plot <- function(count,Group,taxonname){

	GroupTaxon = lapply(unique(Group), function(g) { A = as.matrix(count[, as.factor(Group) %in% g]); index = which(apply(A, 1, sum) != 0); B = as.matrix(A[index,]); rownames(B) = rownames(count)[index]; B} )
	names(GroupTaxon) = unique(Group)
	# Ã¿×éÑù±¾µÄÎïÖÖÁÐ±í
	TaxonInGroup = lapply(GroupTaxon, function(g) data.frame( t(as.matrix(rownames(g)))))
	DataForPlot = lapply(GroupTaxon, function(g) rownames(g))
	

	# Êä³öµ½ÎÄ¼þ
	result = t(do.call(rbind.fill, TaxonInGroup))      
	colnames(result) = paste(taxonname, "In", names(GroupTaxon), sep = ".")
	result = as.data.frame(result)
	if (nrow(result) > 0)  write.csv(result, paste(taxonname, "InGroups.csv", sep = "."), row.names = F, quote = F, na = "") 

	# Ã¿×éÑù±¾ÌØÓÐµÄÎïÖÖÁÐ±í
	AllTaxon = unlist(DataForPlot)
	TaxonOnlyInOneGroup = lapply(DataForPlot, function(x) {T = table(AllTaxon)[x]; names(which(T==1))})
	names = names(TaxonOnlyInOneGroup)

	if( length(names) > 1){

		result = t(do.call(rbind.fill, lapply(TaxonOnlyInOneGroup, function(x) as.data.frame(t(as.matrix(x)))))) 
		if (nrow(result) != 1) result = as.matrix(apply(result, 2, function(r) sort(as.character(r), na.last = T)))      
	
		if (nrow(result) > 0){
	
			colnames(result) = paste(taxonname, "OnlyIn", names, sep = ".")
			result = as.data.frame(result)
			write.csv(result, paste(taxonname, "OnlyInOneGroup.csv", sep = "."), row.names = F, quote = F, na = "")     
	
		}

		# Ã¿×éÑù±¾ÌØÓÐµÄÎïÖÖÊý¾Ý
		for (t in 1:length(TaxonOnlyInOneGroup)){
			output = taxondata[TaxonOnlyInOneGroup[[t]],]
			if ( nrow(output) > 0 ) write.csv(output, paste(taxonname, "DataOnlyInGroup", names[t], "csv", sep = "."), quote = F, na = "")
		}
	
		# ×é¼äÑù±¾¹²ÓÐµÄÎïÖÖ
		common = Reduce(intersect, DataForPlot)
		common = as.matrix(common)
		colnames(common) = paste(taxonname, "Common", sep = ".")
		write.csv(common, paste(taxonname, "Common.csv", sep = "."), row.names = F, quote = F)      
		commonData = taxondata[common,]
		write.csv(commonData, paste(taxonname, "DataCommon.csv", sep = "."), row.names = T, quote = F)

		# Venn plot
		mycol = c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
		mycol = colors()[rep(mycol,20)]
		#str(DataForPlot)
	
		if( length(DataForPlot) < 6 ) {
			source("/home/zhengy/bin/modules/script/venn.diagram_new.r")
			venn.diagram.new(DataForPlot, file = paste(taxonname, "vennDiagram.pdf", sep = "."), fill = mycol[1:length(DataForPlot)])
		}else{
			source("/home/zhengy/bin/modules/script/flower_plot.r")
			a = lengths(TaxonOnlyInOneGroup)
			pdf(paste(taxonname, ".flower.pdf", sep = ""), width = 10, height = 10)
			flower_plot(names(a), as.character(a), 90, 1, 2 ,circle_text = length(common), clockwise = T, ellipse_col = mycol[1:length(DataForPlot)])
		}


		# delete log files
		logfile = list.files(pattern = "*.log", full.names = T)
		for (i in logfile) system(paste("rm -rf", i, sep = " "))   


	}else{

		stop("ERROR: *** Groups number less than 2 ! Please Cheack It! ***\n")

	}

}

if(length(Info) > 5){

	name = unlist(strsplit(taxonfile,"/"))
	taxonname = unlist(strsplit(name[length(name)], "\\."))[1]
	Group = as.character(SampleInfo$Group)
	venn_plot(count,Group,taxonname)

}else{

	for(m in 2:(length(Info))){

		xg = combn(Info, m)
		for (n in 1:ncol(xg)){

			Sample = xg[,n]
			newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
			Group = as.character(newSample[,2])
			#####创建差异分析的文件夹###########
			taxonname = paste(unique(Group), collapse = "_vs_")
			dir = paste(outputpath, taxonname ,sep = "/")
			system(paste("mkdir -p ", dir, sep = ""))
			setwd(dir)
			#venn_plot(count,Group,taxonname)
			countdata = count[, as.character(newSample[,1])]
			GroupTaxon = lapply(unique(Group), function(g) { A = as.matrix(countdata[, as.factor(Group) %in% g]); index = which(apply(A, 1, sum) != 0); B = as.matrix(A[index,]); rownames(B) = rownames(countdata)[index]; B} )
			names(GroupTaxon) = unique(Group)
			# Ã¿×éÑù±¾µÄÎïÖÖÁÐ±í
			TaxonInGroup = lapply(GroupTaxon, function(g) data.frame( t(as.matrix(rownames(g)))))
			DataForPlot = lapply(GroupTaxon, function(g) rownames(g))
			

			# Êä³öµ½ÎÄ¼þ
			result = t(do.call(rbind.fill, TaxonInGroup))      
			colnames(result) = paste(taxonname, "In", names(GroupTaxon), sep = ".")
			result = as.data.frame(result)
			if (nrow(result) > 0)  write.csv(result, paste(taxonname, "InGroups.csv", sep = "."), row.names = F, quote = F, na = "") 

			# Ã¿×éÑù±¾ÌØÓÐµÄÎïÖÖÁÐ±í
			AllTaxon = unlist(DataForPlot)
			TaxonOnlyInOneGroup = lapply(DataForPlot, function(x) {T = table(AllTaxon)[x]; names(which(T==1))})
			names = names(TaxonOnlyInOneGroup)

			if( length(names) > 1){

				result = t(do.call(rbind.fill, lapply(TaxonOnlyInOneGroup, function(x) as.data.frame(t(as.matrix(x)))))) 
				if (nrow(result) != 1) result = as.matrix(apply(result, 2, function(r) sort(as.character(r), na.last = T)))      
			
				if (nrow(result) > 0){
			
					colnames(result) = paste(taxonname, "OnlyIn", names, sep = ".")
					result = as.data.frame(result)
					write.csv(result, paste(taxonname, "OnlyInOneGroup.csv", sep = "."), row.names = F, quote = F, na = "")     
			
				}

				# Ã¿×éÑù±¾ÌØÓÐµÄÎïÖÖÊý¾Ý
				for (t in 1:length(TaxonOnlyInOneGroup)){
					output = taxondata[TaxonOnlyInOneGroup[[t]],]
					if ( nrow(output) > 0 ) write.csv(output, paste(taxonname, "DataOnlyInGroup", names[t], "csv", sep = "."), quote = F, na = "")
				}
			
				# ×é¼äÑù±¾¹²ÓÐµÄÎïÖÖ
				common = Reduce(intersect, DataForPlot)
				common = as.matrix(common)
				colnames(common) = paste(taxonname, "Common", sep = ".")
				write.csv(common, paste(taxonname, "Common.csv", sep = "."), row.names = F, quote = F)      
				commonData = taxondata[common,]
				write.csv(commonData, paste(taxonname, "DataCommon.csv", sep = "."), row.names = T, quote = F)

				# Venn plot
				mycol = c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
				mycol = colors()[rep(mycol,20)]
				#str(DataForPlot)
			
				if( length(DataForPlot) < 6 ) {
					source("/home/zhengy/bin/modules/script/venn.diagram_new.r")
					venn.diagram.new(DataForPlot, file = paste(taxonname, "vennDiagram.pdf", sep = "."), fill = mycol[1:length(DataForPlot)])
				}else{
					source("/home/zhengy/bin/modules/script/flower_plot.r")
					a = lengths(TaxonOnlyInOneGroup)
					pdf(paste(taxonname, ".flower.pdf", sep = ""), width = 10, height = 10)
					flower_plot(names(a), as.character(a), 90, 1, 2 ,circle_text = length(common), clockwise = T, ellipse_col = mycol[1:length(DataForPlot)])
				}


				# delete log files
				logfile = list.files(pattern = "*.log", full.names = T)
				for (i in logfile) system(paste("rm -rf", i, sep = " "))   


			}else{

				stop("ERROR: *** Groups number less than 2 ! Please Cheack It! ***\n")

			}
		}
	}	
}



