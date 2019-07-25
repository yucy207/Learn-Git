rm(list=ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	"============== Venn ===========
	Usage: 
	Rscript 16S.BetaDiversity.Venn.r taxonfile groupfile outputpath scriptdir
		parameters ->
			taxonfile: [file -> always subsample_otu.tax.0.03.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		    outputpath: [path -> path for output];
		    scriptdir：[path -> path of sub_script, such as :/home/genesky/pipeline/absolute_quantification_16s/v1.1/sub_script];\n")
	options("show.error.messages" = F) 
	stop()

}

library(plyr)
library(VennDiagram)

taxonfile <- normalizePath(args[1])
groupfile <- normalizePath(args[2])
# select    <- args[3]
outputpath <- normalizePath(args[3])
scriptdir <- args[4]

setwd(outputpath)

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

venn_data <- function(count,Group,taxonname){

	GroupTaxon <- lapply(unique(Group), function(g) { A = as.matrix(count[, factor(Group, levels = unique(Group)) %in% g]); index = which(apply(A, 1, sum) != 0); B = as.matrix(A[index,]); rownames(B) = rownames(count)[index]; B} )
	names(GroupTaxon) <- unique(Group)

	TaxonInGroup <- lapply(GroupTaxon, function(g) data.frame( t(as.matrix(rownames(g)))))
	DataForPlot <<- lapply(GroupTaxon, function(g) rownames(g))	

	result <- t(do.call(rbind.fill, TaxonInGroup))      
	colnames(result) <- paste(taxonname, "In", names(GroupTaxon), sep = ".")
	result <- as.data.frame(result)
	if (nrow(result) > 0)  write.csv(result, paste(taxonname, "InGroups.csv", sep = "."), row.names = F, quote = F, na = "") 

	AllTaxon <- unlist(DataForPlot)
	TaxonOnlyInOneGroup <- lapply(DataForPlot, function(x) {T = table(AllTaxon)[x]; names(which(T==1))})
	names <- names(TaxonOnlyInOneGroup)

	result <- t(do.call(rbind.fill, lapply(TaxonOnlyInOneGroup, function(x) as.data.frame(t(as.matrix(x)))))) 
	if (nrow(result) != 1) result <- as.matrix(apply(result, 2, function(r) sort(as.character(r), na.last = T)))      	
	if (nrow(result) > 0){
	
		colnames(result) <- paste(taxonname, "OnlyIn", names, sep = ".")
		result <- as.data.frame(result)
		write.csv(result, paste(taxonname, "OnlyInOneGroup.csv", sep = "."), row.names = F, quote = F, na = "")     
	
	}

	for (t in 1:length(TaxonOnlyInOneGroup)){
		output <- taxondata[TaxonOnlyInOneGroup[[t]],]
		if ( nrow(output) > 0 ) write.csv(output, paste(taxonname, "DataOnlyInGroup", names[t], "csv", sep = "."), quote = F, na = "")
	}
	
	common <- Reduce(intersect, DataForPlot)
	common <- as.matrix(common)
	colnames(common) <- paste(taxonname, "Common", sep = ".")
	write.csv(common, paste(taxonname, "Common.csv", sep = "."), row.names = F, quote = F)      
	commonData <- taxondata[common,]
	write.csv(commonData, paste(taxonname, "DataCommon.csv", sep = "."), row.names = T, quote = F)

	Result <- list(TaxonOnlyInOneGroup, DataForPlot,common)
	return(Result)
	
}


mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol <- colors()[rep(mycol,20)]
group2col = data.frame(unique(as.character(SampleInfo[[2]])),mycol,stringsAsFactors=FALSE)
colnames(group2col) = c("gro","col")


# if(nchar(select) > 0){

# 	##指定两两分析
# 	pair = unlist(strsplit(select,";"))
# 	# print(pair)
# 	for (i in 1:length(pair)) {
# 		Sample = unlist(strsplit(pair[i], ","))
# 		newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
# 		newGroup = as.character(newSample[,2])
# 		newmycol = group2col$col[group2col$gro %in% unique(newGroup)]
# 		#####创建差异分析的文件夹###########
# 		taxonname <- paste(unique(newGroup), collapse = "_vs_")
# 		dir <- paste(outputpath, taxonname ,sep = "/")
# 		system(paste("mkdir -p ", dir, sep = ""))
# 		setwd(dir)

# 		countdata <- count[, as.character(newSample[,1])] 
# 		Result <- venn_data(countdata,newGroup,taxonname)
	
# 		source(paste(scriptdir,"venn.diagram_new.r", sep = "/"))
# 		venn.diagram.new(Result[[2]], file = paste(taxonname, "vennDiagram.pdf", sep = "."), fill = newmycol)	

# 	}	



# 	##所有样本一起分析
# 	Group = as.character(SampleInfo[,2])
# 	col = group2col$col[group2col$gro %in% unique(Group)]
# 	taxonname = "all_group"
# 	dir = paste(outputpath, "all_group", sep = "/")
# 	system(paste("mkdir -p ", dir, sep = ""))
# 	setwd(dir)	

# 	Result <- venn_data(count,Group,taxonname)

# 	if(length(Info) > 5){
# 		source(paste(scriptdir,"flower_plot.r", sep = "/"))
# 		a <- lengths(Result[[1]])
# 		pdf(paste(taxonname, ".flower.pdf", sep = ""), width = 10, height = 10)
# 		flower_plot(names(a), as.character(a), 90, 1, 2 ,circle_text = length(Result[[3]]), clockwise = T, ellipse_col = col)	

# 	}else{

# 		source(paste(scriptdir,"venn.diagram_new.r", sep = "/"))
# 		venn.diagram.new(Result[[2]], file = paste(taxonname, "vennDiagram.pdf", sep = "."), fill = col)
# 	}


# }else{


	if(length(Info) > 5){

		name = unlist(strsplit(taxonfile,"/"))
		taxonname = unlist(strsplit(name[length(name)], "\\."))[1]
		Group = as.character(SampleInfo$Group)
		Result <- venn_data(count,Group,taxonname)

		source(paste(scriptdir,"flower_plot.r", sep = "/"))
		a <- lengths(Result[[1]])
		pdf(paste(taxonname, ".flower.pdf", sep = ""), width = 10, height = 10)
		flower_plot(names(a), as.character(a), 90, 1, 2 ,circle_text = length(Result[[3]]), clockwise = T, ellipse_col = mycol[1:length(Result[[2]])])

	}else{

		for(m in 2:(length(Info))){

			xg = combn(Info, m)
			for (n in 1:ncol(xg)){

				Sample = xg[,n]
				newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
				Group = as.character(newSample[,2])

				if(length(unique(Group)) == 4){
					col = group2col$col[group2col$gro %in% unique(Group)]
					# col = col[c(3,4,1,2)]
					col = col[c(1,4,2,3)]
				}else{
					col = group2col$col[group2col$gro %in% unique(Group)]
				}			

				#####创建差异分析的文件夹###########
				taxonname <- paste(unique(Group), collapse = "_vs_")
				dir <- paste(outputpath, taxonname ,sep = "/")
				system(paste("mkdir -p ", dir, sep = ""))
				setwd(dir)

				countdata <- count[, as.character(newSample[,1])] 
				Result <- venn_data(countdata,Group,taxonname)
			
				source(paste(scriptdir,"venn.diagram_new.r", sep = "/"))
				venn.diagram.new(Result[[2]], file = paste(taxonname, "vennDiagram.pdf", sep = "."), fill = col)


			}
		}	
	}


# }


