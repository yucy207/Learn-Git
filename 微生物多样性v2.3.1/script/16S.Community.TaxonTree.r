# 物种分类树状图
rm(list = ls())

args = commandArgs(trailingOnly=T) 
if(length(args) != 2){

	cat("Usage: 
	Rscript 16S.Community.TaxonTree.r otustatfile outputpath
		parameters -> 
			otustatfile: [file -> always otu.tax.0.03.stat.xls derived from 16S.OTU.Modify.pl];
		         outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
}
 
# https://cran.r-project.org/src/contrib/mapplots_1.5.tar.gz install by R CMD
library(mapplots)

otustat = normalizePath(args[1])
outputpath = normalizePath(args[2])

# otustat = normalizePath("/home/yangkj/work/metagenomics/16B1102C/report/OTU/otu.tax.0.03.stat.xls")
# outputpath = normalizePath("/home/yangkj/work/metagenomics/16B1102C/report/Community/TaxonTree")

# otustat = normalizePath("/home/yangkj/work/metagenomics/16B1102C/report/OTU/otu.tax.0.03.stat.xls")
# outputpath = normalizePath("/home/zhengxq/workdir/meta_genomics/16C1025A/ITS2/report/Community/TaxonTree")

# =============================== Data prepare ===============================
otu_stat = read.table(otustat, header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "", quote = "")     # 970*25

rownames(otu_stat) = gsub("\\(", "", rownames(otu_stat))
rownames(otu_stat) = gsub("\\)", "", rownames(otu_stat))
#rownames(otu_stat) = gsub("\\[", "", rownames(otu_stat))
#rownames(otu_stat) = gsub("\\]", "", rownames(otu_stat))
rownames(otu_stat) = gsub("\\+", "", rownames(otu_stat))
rownames(otu_stat) = gsub("\\^", "", rownames(otu_stat))
		
loc = grep("Size|Abundance", colnames(otu_stat))

if (length(loc) > 0){

	data = as.matrix(otu_stat[,-loc])
	colnames(data) = colnames(otu_stat)[-loc]
	rownames(data) = rownames(otu_stat)
	
}else{

	data = as.matrix(otu_stat)

}

classi = strsplit(rownames(data), split="\\|")    # 转义，\也需要转义
classi_len = as.matrix(lapply(classi, length))    # 分类注释层数  7

# 识别分类及层数
ind = which(classi_len == max(as.numeric(classi_len)))[1]
tax = unlist(lapply(strsplit(classi[[ind]], split="_"),function(x) x[1]))
taxon = grep(paste(toupper(tax), collapse="|"), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), value=T)

# data中每层分类所在的位置
indexs = lapply(taxon, function(t) which(classi_len == grep(t,taxon)))
names(indexs) = taxon

# 相对丰度
proportions = lapply(indexs, function(i) {

	A = as.matrix(data[i,])
	
	if (nrow(A) != 1){
	
		if(rownames(A)[1] == colnames(data)[1]){		# 只有一个物种，被转置
	
			A = t(A)
		
		}
		
	}
	
	rownames(A) = rownames(data)[i]
	apply(A, 1, function(x) sum(x)/sum(data[indexs$Kingdom,]))
	
})

# 相对丰度大于1%的分类
indexs_for_plot = lapply(proportions, function(p) which(p > 0.01))

# 相对丰度大于1%的分类所对应的原始丰度值
data_for_plot = lapply(indexs_for_plot, function(i) {

	A = as.matrix(data[names(i),])
	
	if (nrow(A) != 1){
	
		if(rownames(A)[1] == colnames(data)[1]){		# 只有一个物种，被转置
	
			A = t(A)
		
		}
		
	}
	
	rownames(A) = names(i)
	return(A)
	
})

# 每层用于绘图的分类
taxons = lapply(indexs_for_plot, function(i)  unique(unlist(lapply(strsplit(names(i), split="__|\\|"), function(t) t[length(t)]))))

# 每层用于绘图分类的上层分类
taxons_Up_taxons = lapply(indexs_for_plot[-1], function(i)  unique(unlist(lapply(strsplit(names(i), split="__|\\|"), function(t) t[length(t)-2]))))
names(taxons_Up_taxons) = taxon[-length(taxon)]

# 添加无子分类点的分类
taxonsForPlot = rownames(as.matrix(data_for_plot[[length(taxon)]]))

for (i in 1:length(taxons_Up_taxons)){
  
	diff = setdiff(taxons[[i]], taxons_Up_taxons[[i]])
	if (length(diff) !=0 ){
		AddPlot = grep(paste(diff,collapse="|"), rownames(as.matrix(data_for_plot[[i]])), value=T)
	}else{
		AddPlot = NULL
	}
	taxonsForPlot = c(taxonsForPlot, AddPlot)
  
}

# 用于绘图的分类树
taxonsForPlot = sort(taxonsForPlot)

# =============================== Plot ===============================

mycol = c(552,81,142,53,460,566,107,476,258,551,576,30,426,430,564,561,619,84,521,644,10,46,51,53,128,101,371,373,386,435,429,528,512,375)
mycol = colors()[rep(mycol,20)]
color = mycol[1:ncol(data)] 

# 根据最终树的行数定义图像高度
height = 2.2*length(taxonsForPlot)    
if (height < 10)  height = 20 

# 根据最长字符串宽度及注释层数设定图像宽度
charlen = max(nchar(unlist(strsplit(taxonsForPlot, split = "\\|")))) 
classilen = as.numeric(lapply(strsplit(taxonsForPlot, split = "\\|"), length))
width = max(charlen)*max(classilen)/8 
if (width < 10)   width=10

# 用于记录所有注释点的位置
point_position = vector("list", length(taxon))    # 定义一个特定长度的空列表
names(point_position) = taxon

# 未注释到最底层的分类
naloc = which(classilen < max(classilen))
# 最底层注释所在的横坐标位置
xloc = (2*max(classilen)+1)/(2*max(classilen)+2)*width        # 每个饼图所在位置符合(2n+1)/(2n+2)
for (p in 1:length(taxon)) point_position[[p]][naloc] = paste(xloc,2*naloc,sep=",")

# 绘图
pdf(paste(outputpath, "TaxonTree.pdf", sep = "/"),width = width,height = height)

# 空画布
plot(seq(0, width, length = 10), seq(0, height, length = 10), axes = F, # without axis
	type = "n",   # plot nothing
    xlab = "", ylab = "")   # without lab 

if(length(colnames(data)) > 1)  legend("left", legend = colnames(data), bty = "n", fill = color, cex = 1.5)

terms = lapply(lapply(taxonsForPlot, function(t)  unlist(strsplit(t, split = "\\|\\w__"))), function(x) {x[1] = gsub("\\w__","",x[1]); return(x)})

ids = vector("list", length(taxon))    # 用于记录分类的index
taxonnames = vector("list", length(taxon))      # 用于记录分类的名称
position = vector("list", length(taxon))    # 用于记录相同分类的index，方便后续取均值

# 最底层注释
bottom = length(taxon)
ids[[bottom]] = which(length(taxon) <= classilen)
taxonnames[[bottom]] = unlist(lapply(ids[[bottom]], function(x) terms[[x]][bottom]))

for (b in ids[[bottom]]){
  
	x = (2*bottom+1)/(2*bottom+2)*width
	y = 2*b
  
	segments(2*bottom/(2*bottom+2)*width, y, x, y, lwd=1)
	points(x, y, pch = 16, cex = 0.7)
	
	size = as.numeric(data_for_plot[[bottom]][taxonsForPlot[b],])
	
	if(ncol(data_for_plot[[bottom]]) != 1){
		
		add.pie(size, x = x, y = y+0.7, radius = 0.6, labels = "", border = NA, col = color)
		
	}else{
	
		add.pie(size, x = x, y = y+0.08+size^(1/5)/20, radius = size^(1/5)/20, labels = "", clockwise = F, border = "cornflowerblue", col = "cornflowerblue")
		
	}
	
	text(x = x, y = y-0.45, labels = sum(size), cex = 1.3)
	text(x = x, y = y-0.18, labels = unlist(strsplit(taxonsForPlot[b],split = "\\|\\w__"))[bottom], cex = 1.3)
	point_position[[bottom]][b]  =  paste(x, y, sep = ",")
  
}

# 上层注释
for (i in (length(taxon)-1):1) {
  
	x = (2*i+1)/(2*bottom+2)*width    # 横坐标定位
  
	# 位置和分类名称
	ids[[i]] = which(i <= classilen)
	taxonnames[[i]] = unlist(lapply(ids[[i]], function(x) terms[[x]][i]))
  
	# 对于相同分类项
	for ( Repeat in names(table(taxonnames[[i]])[which(table(taxonnames[[i]]) > 1)]) ){
		
		index = grep(Repeat, taxonsForPlot)     # 相同分类项所在的index
    
		# 如果相同类别的父节点注释也相同，则可以合并
		if ( length(unique(unlist(lapply(strsplit(taxonsForPlot[index], split = Repeat), function(x) x[1])))) == 1 ){
      
			position[[i]] = c(position[[i]], index)

			# 所有子节点位置的均值
			miny = as.numeric(unlist(strsplit(point_position[[i+1]][min(index)], ","))[2])
			maxy = as.numeric(unlist(strsplit(point_position[[i+1]][max(index)], ","))[2])
			y = mean( c(miny, maxy) )
      
			segments((2*i+2)/(2*bottom+2)*width, miny, (2*i+2)/(2*bottom+2)*width, maxy, lwd = 1)
      
			# 如果到root节点，则左边不再多加半条线
			if(i == 1){
				segments((2*i+1)/(2*bottom+2)*width, y, (2*i+2)/(2*bottom+2)*width, y, lwd = 1)
				tex = paste(Repeat, "$", sep = "")
			}else{
				segments((2*i)/(2*bottom+2)*width, y, (2*i+2)/(2*bottom+2)*width, y, lwd = 1)
				tex = paste("_", Repeat, "$", sep = "")
			}
      
			points(x, y, pch = 16, cex = 0.7)
			size = as.numeric(data_for_plot[[i]][grep(tex, rownames(data_for_plot[[i]]), value = T),])

			if(ncol(data_for_plot[[bottom]]) != 1){
			
				add.pie(size, x = x, y = y+0.7, radius = 0.6, labels = "", border = NA, col = color)
      
			}else{
				
				radius = sum(size)^(1/5)/20
				add.pie(sum(size), x = x, y = y+0.08+radius, radius = radius, labels = "", clockwise = F, border = "cornflowerblue", col = "cornflowerblue")
			
			}
			
			text(x = x, y = y-0.45, labels = sum(size), cex = 1.3)	
			text(x = x, y = y-0.18, labels = Repeat, cex = 1.3)
      
			point_position[[i]][index] = paste(x, y, sep = ",")
     
		}

	}

	for ( p in setdiff(ids[[i]],position[[i]]) ){
  
		# 与子节点纵坐标相同
		y  =  as.numeric(unlist(strsplit(point_position[[i+1]][p],","))[2])
    
		# 如果到root节点，则左边不再多加半条线
		if(i == 1){
      
			if (length(terms[[p]]) == i){      # 如果不存在子节点
				segments((2*i+1)/(2*bottom+2)*width, y, (2*i+1)/(2*bottom+2)*width, y, lwd = 1)
			}else{
				segments((2*i+1)/(2*bottom+2)*width, y, (2*i+2)/(2*bottom+2)*width, y, lwd = 1)
			}
      
		}else{
      
			if (length(terms[[p]]) == i){      # 如果不存在子节点
				segments((2*i)/(2*bottom+2)*width, y, (2*i+1)/(2*bottom+2)*width, y, lwd = 1)
			}else{
				segments((2*i)/(2*bottom+2)*width, y, (2*i+2)/(2*bottom+2)*width, y, lwd = 1)
			}
      
		}

		tex2 = paste(paste(terms[[p]][1:i],collapse = ".*"),"$",sep = "")
		
		size = as.numeric(data_for_plot[[i]][grep(tex2, rownames(data_for_plot[[i]]), value = T),])
		
		if(ncol(data_for_plot[[bottom]]) != 1){
		
			add.pie(size, x = x, y = y+0.7, radius = 0.6, labels = "", border = NA, col = color)

		}else{

			radius = sum(size)^(1/5)/20
			add.pie(sum(size), x = x, y = y+0.08+radius, radius = radius, labels = "", clockwise = F, border = "cornflowerblue", col = "cornflowerblue")
		
		}
		
		text(x = x, y = y-0.45, labels = sum(size), cex = 1.3)
		text(x = x, y = y-0.18, labels = terms[[p]][i], cex = 1.3)

		point_position[[i]][p] = paste(x, y, sep = ",")
  
	}
  
}

dev.off()
