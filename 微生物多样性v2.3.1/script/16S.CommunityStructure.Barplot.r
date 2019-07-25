# 16S.CommunityStructurePlot.r
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
	"Usage: 
	Rscript 16S.CommunityStructure.Barplot.r taxonfile outputpath
		parameters ->
			  taxonfile: [file -> always *.taxon.Abundance.xls with row taxon and column samples];
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

taxonfile = normalizePath(args[1])
outputpath = normalizePath(args[2])

# taxonfile = normalizePath("/home/yangkj/work/metagenomics/16B1102C/report/Community/Community_Structure/class/class.taxon.Abundance.xls")
# outputpath = normalizePath("/home/yangkj/work/metagenomics/16B1102C/report/Community/Community_Structure/class")

inputdata = read.table(taxonfile, header=TRUE, row.names=1, check.name=F, comment.char="", quote="", sep="\t", fill=T)    # 注意 quote="" 

rownames(inputdata) = gsub("\\(", "", rownames(inputdata))
rownames(inputdata) = gsub("\\)", "", rownames(inputdata))
rownames(inputdata) = gsub("\\[", "", rownames(inputdata))
rownames(inputdata) = gsub("\\]", "", rownames(inputdata))
rownames(inputdata) = gsub("\\+", "", rownames(inputdata))
rownames(inputdata) = gsub("\\^", "", rownames(inputdata))

# 如果存在分类注释信息,则去除
index = grep("OTUsize|Abundance", colnames(inputdata))
if(length(index) > 0){

	Data = as.matrix(inputdata[,1:(index-1)])
	colnames(Data) = colnames(inputdata)[1:(index-1)]

}
rownames(Data) = rownames(inputdata)

# prepare data
rownames(Data) = gsub("\"","",rownames(Data))     # 去除名字中的双引号,全局替换
colnames(Data) = gsub("\"","",colnames(Data))
Data = as.matrix(Data)
rownames(Data) = as.character(strsplit(rownames(Data),split="{.*}",perl=T))

loc = which(rownames(Data) == "All")       # 去除 All 行

if(length(loc) > 0) Data1 = as.matrix(Data[-loc, ])

if( ncol(Data1) == 1){

	Data1 = as.matrix(Data1[Data1 != 0, ])

}else{

	Data1 = as.matrix(Data1[rowSums(Data1) != 0, ])

}
if( !is.null(nrow(Data1)) ){     # 只有一类的,不画

	# get relative.abundance
	su = apply(Data1, 1, sum)
	ra = su/sum(su)
	loc = which(rownames(Data1) %in% c("No_Rank","Unclassified"))     # 不画 No_Rank Unclassified
	if(length(loc) > 0) {
		Data1 = Data1[-loc,]
		ra = ra[-loc]
	}

	# 种类少于30个的全画，多于30个的，取丰度最高的30个进行展示！
	# set color
	mycol = c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
	mycol = colors()[rep(mycol,20)]
	
	# 不同水平的legend
	if (length(index) < ncol(inputdata) && length(index) !=0 ){       #  有上层注释

		for(n in (index+1):ncol(inputdata)){      # 添加不同水平的注释
	
			if(length(ra) < 30){
				data_for_plot = ra
			}else{
				data_for_plot = ra[1:30]
			}
			
			chars = inputdata[names(data_for_plot), n]
			chars = inputdata[names(data_for_plot), n]
			
			chars = gsub("\\(", "", chars)
			chars = gsub("\\)", "", chars)
			chars = gsub("\\[", "", chars)
			chars = gsub("\\]", "", chars)
			chars = gsub("\\+", "", chars)
			chars = gsub("\\^", "", chars)
			
			label = unique(chars)
			
			# %
			data_for_plot = data_for_plot*100
		
			color = mycol[1 : length(label)] 
			colorss = rep(0, length(data_for_plot))
			for ( i in 1:length(color) )   colorss[grep(label[i], as.character(chars))] = color[i]    # 设置颜色

			# set width and height
			mn = max(nchar(names(data_for_plot)))
	
			width = ifelse(length(data_for_plot)*0.3 < 6, 6, length(data_for_plot)*0.3)
			height = 4 + 0.1*mn 

			# ymax = ifelse(max(data_for_plot)+5 < 100,  max(data_for_plot) + 5, 100)
			# ymax = ceiling(max(data_for_plot)/10)*10
			ymax = ifelse(ceiling(max(data_for_plot)%%10) < 5, (ceiling(max(data_for_plot)/10)*10)-5, ceiling(max(data_for_plot)/10)*10)

			name = unlist(strsplit(taxonfile,"/"))
			
			# 不同上层水平注释（与legend搭配使用！）
			pdf(paste(outputpath, paste(strsplit(name[length(name)], "\\.\\w+\\.xls"), "Barplot.Legend.With", colnames(inputdata)[n], "pdf", sep = "."), sep="/"), width = width, height = height)
			# 不加注释
			#pdf(paste(args[2], paste(strsplit(name[length(name)], "\\.\\w+\\.xls"), "Barplot.pdf", sep = "."), sep="/"), width = width, height = height)
			n = ifelse(0.2*nchar(names(data_for_plot)[1]) < 2.5, 2.5, 0.2*nchar(names(data_for_plot)[1]))
			par(mar=c(0.2*mn, n, 1, 1), mgp=c(1.5, 0.5, 0))		
			x = barplot(data_for_plot, space = 1,     # 条形间距
					  beside = FALSE,   # beside 堆积条形图
				   	  border = NA,      # 无条形边框
				   axisnames = FALSE,   # 不加坐标上的label
					   horiz = FALSE, density = NULL, col = colorss, xlab = "", ylab = "Relative abundance (%)", axes = TRUE, cex.main = 1, cex.lab = 0.9, cex.axis = 0.8, ylim=c(0, ymax))
			text(x, y = -ymax/40, labels = names(data_for_plot), adj = 1, srt = 35, xpd = TRUE, cex = 0.6)
			legend("topright", legend = label, fill = color, bty = "n", cex = 0.6)
			dev.off()
		
		}
	
	}else{          # superkingdom 无上层注释
	
		if(length(ra) < 30){
			data_for_plot = ra
		}else{
			data_for_plot = ra[1:30]
		}
		
		# %
		data_for_plot = data_for_plot*100
	
		color = mycol[1 : length(data_for_plot)] 

		# set width and height
		mn = max(nchar(names(data_for_plot)))

		width = ifelse(length(data_for_plot)*0.3 < 5, 5, length(data_for_plot)*0.3)
		height = 5 + 0.1*mn 

		# ymax = ifelse( max(data_for_plot)+5 < 100,  max(data_for_plot) + 5, 100)
		# ymax = ceiling(max(data_for_plot)/10)*10
		ymax = ifelse(ceiling(max(data_for_plot)%%10) < 5, (ceiling(max(data_for_plot)/10)*10)-5, ceiling(max(data_for_plot)/10)*10)

		name = unlist(strsplit(taxonfile,"/"))
		pdf(paste(args[2], paste(strsplit(name[length(name)], "\\.\\w+\\.xls"), "Barplot.pdf", sep = "."), sep="/"), width = width, height = height)
		n = ifelse(0.2*nchar(names(data_for_plot)[1]) < 2.5, 2.5, 0.2*nchar(names(data_for_plot)[1]))    # 到左边界的距离与第一个字符串长度有关
		par(mar=c(0.2*mn, n, 1, 1), mgp=c(1.5, 0.5, 0))	
		x = barplot(data_for_plot, space = 1,     # 条形间距
					beside = FALSE,   # beside 堆积条形图
				   	border = NA,      # 无条形边框
				    axisnames = FALSE,   # 不加坐标上的label
					horiz = FALSE, density = NULL, xlab = "", ylab = "Relative abundance (%)", axes = TRUE, cex.main = 1, cex.lab = 0.9, cex.axis = 0.8, ylim=c(0, ymax))
		text(x, y = -ymax/40, labels = names(data_for_plot), adj = 1, srt = 35, xpd = TRUE, cex = 0.6)
		dev.off()

	}

}


