rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	"============== do Picrust different_analysis ===========
	Usage: 
	Rscript 16S.PI_diff_analysis.r cogfile keggfile groupfile outputpath
		parameters ->
			  cogfile: [file -> always cog_for_stamp.txt];
			  keggfile: [file -> always kegg_for_stamp.txt];
			  groupfile: [file -> always sample.groups];
		      outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

cogfile = normalizePath(args[1])
keggfile = normalizePath(args[2])
groupfile = normalizePath(args[3])
outputpath = normalizePath(args[4])

#library(ggplot2)
setwd(outputpath)
dir1 = paste(outputpath, "different_analysis",sep = "/")
system(paste("mkdir -p ", dir1, sep = ""))

SampleInfo = read.table(groupfile, header = F, sep = "\t",quote="")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Groups = unique(as.character(SampleInfo [,2]))


################################## 多组样本间差异ANOVA ##################################################
anova_analysis <-function(inputfile,SampleInfo){
	inputdata = read.table(inputfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote = "", fill = T)
	name = unlist(strsplit(unlist(strsplit(inputfile,"/"))[length(unlist(strsplit(inputfile,"/")))], "\\."))[1]
	# 样本信息
	#SampleInfo = read.table(groupfile, header = F, sep="\t")
	#rownames(SampleInfo) = SampleInfo[,1]
	#colnames(SampleInfo) = c("Sample", "Group")
	# prepare data
	count = inputdata[,rownames(SampleInfo)]
	rownames(count) = gsub("\"", "", rownames(count))   # 去除名字中的双引号,全局替换
	colnames(count) = gsub("\"", "", colnames(count))
	rownames(count) = as.character(strsplit(rownames(count), split = "{.*}", perl = T))

	count = count[rowSums(count) != 0 , ]

	loc = which(rownames(count) == "All")       # 去除 All 行
	if(length(loc) > 0) count = count[-loc, ]

	# 相对丰度
	re.count = apply(count, 2, function(d) d/sum(as.numeric(d)))
	#index = which(rowSums(re.count) > 0.001)
	#re.count = re.count[index, ]
	#count = count[index,]
	
	# compare
	DataForDEA = re.count[rownames(count),]
	
	res = sapply(1:nrow(DataForDEA), function(g) {
		mod = lm(DataForDEA[g,] ~ SampleInfo$Group)
		am = anova(mod)
		F1 = am[4][1, 1]   # F value
		P1 = am[5][1, 1]   # P value
		FP = cbind(F1, P1)
	})

	res = t(res)
	rownames(res) = rownames(DataForDEA)
	colnames(res) = c("Fvalue", "Pvalue")	

	fdr = p.adjust(res[,2], method = "BH")

	means = sapply(Groups, function(g) apply(DataForDEA[,as.character(SampleInfo[as.character(SampleInfo$Group) %in% g, 1])], 1, mean))
	colnames(means) = paste("Mean.In.", Groups, sep = "")
	
	result = data.frame(Taxon = rownames(DataForDEA), means, res, FDR = fdr)
	result = result[order(result$Pvalue),]
	write.table(result,paste(name, ".ANOVA_test_result.xls", sep = ""), sep = "\t", row.names = F, quote = F)

}


if(length(Groups) > 2){		# 多组样本间差异 ANOVA
	anova_analysis(cogfile,SampleInfo)
	anova_analysis(keggfile,SampleInfo)
}
#########################################################################################################


t_test <- function(data.x, data.y){

	res <- lapply(1:(ncol(data.x)-1), function(t){
  	data.1 <- unlist(data.x[,t])
  	res <- t.test(data.1 ~ data.y, var.equal = TRUE)
   	c(colnames(data.x)[t],res$estimate[1],res$estimate[2], res$conf.int[1],res$conf.int[2],mean(res$conf.int), res$p.value)
 })

 result <- data.frame(do.call(rbind, res))
 return(result)
}


diff_analysis <-function(inputfile,SampleInfo1,type,m,n){

	data = read.table(inputfile, header = T, sep = "\t", row.names = 1,check.name = F, comment.char = "", quote = "", fill = T)
	data1 = data [,as.character(SampleInfo1$Sample)]

	#去除和为0的行
	SUM = as.numeric(apply(data1,1,sum))
	data1 = data1[SUM!=0,]

	#获得相对丰度 
	Re = sapply(1:ncol(data1), function(x) 100*data1[,x]/sum(data1[,x]))
	rownames(Re) = rownames(data1)
	colnames(Re) = colnames(data1)

	#t.test
	new = as.data.frame(t(Re[,as.character(SampleInfo1$Sample)]))
	new$group = SampleInfo1$Group
	result = t_test(new,new$group)

	colnames(result) = c("Category",paste("mean",Groups[m], sep="_"),paste("mean",Groups[n], sep="_"),
						 "95.0% lower CI","95.0% upper CI","Difference between means","P.value")

	write.table(result, paste(type, "diff.all.xls",sep = "_"), row.names = F, sep = "\t", quote = F)

	#根据p值过滤
	p = apply(result[7], 1, as.numeric)
	df = result[( p<0.05 & p>0), ]

	if (nrow(df)>0){

		write.table(df, paste(type, "diff.P0.05.xls",sep = "_"), row.names = F, sep = "\t", quote = F)

		#################plot################
		#prepare data for plot
		#添加p值*号
		p_value = as.numeric(as.character(df$P.value))
		y <- numeric(length(p_value))
		y[p_value >= 0.01] = "*"
		y[p_value < 0.01 & p_value >= 0.001] = "**"
		y[p_value < 0.001 & p_value >= 0.0001] = "***"
		y[p_value < 0.0001] = "****"
		name = paste(df$Category,y,sep = "")


		mean1 = apply(df[2], 1, as.numeric)
		mean2 = apply(df[3], 1, as.numeric)
		mean =  as.data.frame(rbind(mean1,mean2))
		colnames(mean) = name
		data = as.matrix(mean)

		bar_cols = c("SkyBlue3","SandyBrown")
		height_barplot = 0.15 * nrow(df) * 2 + 0.15 * (nrow(df) - 1) + 1 + 0.3
		if(nrow(df) == 1){
			height =  height_barplot + 1
		}else{
			height =  height_barplot + 0.6
		}
		

		#layout
		if(type =="COG"){
			width = 10 
			pdf("COG.pdf",width = width ,height = height)

			layout(matrix(c(1,2), nrow = 2, byrow = T),heights = c(1,height_barplot))
			par(mai = c(0.1,3,0.2,3), xaxs = "i", yaxs = "i")
			plot(1:10, type = "n", yaxt = "n", xlab = "", xaxt = "n", ylab = "", bty = "n")

			rect(xleft = 2, ybottom = 4, xright = 3, ytop = 6, col = bar_cols[1],lwd = 1.5)
			rect(xleft = 5, ybottom = 4, xright = 6, ytop = 6, col = bar_cols[2],lwd = 1.5)
			text(x = 3*1.02, y = 5, adj = 0, labels = Groups[m], cex = 1.2)
			text(x = 6*1.02, y = 5, adj = 0, labels = Groups[n], cex = 1.2)

			par(las = 1, mai = c(1, 5.6, 0.3, 0.5),  xaxs = "i", yaxs = "i", mgp = c(3, 2, 1), lend = 2,lwd = 1.5)

		}
		if(type =="KEGG"){
			width = 12 
			pdf("KEGG.pdf",width = width ,height = height)

			layout(matrix(c(1,2), nrow = 2, byrow = T),heights = c(1,height_barplot))
			par(mai = c(0.1,3,0.2,3), xaxs = "i", yaxs = "i")
			plot(1:10, type = "n", yaxt = "n", xlab = "", xaxt = "n", ylab = "", bty = "n")

			rect(xleft = 2, ybottom = 4, xright = 3, ytop = 6, col = bar_cols[1],lwd = 1.5)
			rect(xleft = 5, ybottom = 4, xright = 6, ytop = 6, col = bar_cols[2],lwd = 1.5)
			text(x = 3*1.02, y = 5, adj = 0, labels = Groups[m], cex = 1.5)
			text(x = 6*1.02, y = 5, adj = 0, labels = Groups[n], cex = 1.5)


			par(las = 1, mai = c(1, 5.5, 0.3, 0.5),  xaxs = "i", yaxs = "i", mgp = c(3, 2, 1), lend = 2,lwd = 1.5)
		}


		#barplot
		wd = 0.8
		if(max(data) >= 1){

			x_max <- ceiling(max(data))
			a = barplot(data, beside = T, horiz = T, width = wd ,col = bar_cols, axes = FALSE,xlim = c(0,x_max))
			cols = 1:ncol(a)
			#添加每个柱子的灰度背景
			rect(xleft = 0, ybottom = a[1, cols %% 2 == 1] - wd / 2, xright = par("usr")[2], ytop = a[2, cols %% 2 == 1] + wd / 2  , col = "gray91", border = NA)
			barplot(data, beside = T, horiz = T, width = wd, col = bar_cols , add = TRUE, axes = FALSE,xlim = c(0,x_max))
			axis(side = 1, at = c(0, x_max), labels = c(0, x_max))

		}else{

			x_max<-format(max(data),scientific=TRUE,digit=2)
			a = barplot(data, beside = T, horiz = T, width = wd ,col = bar_cols, axes = FALSE,xlim = c(0,max(data)))
			cols = 1:ncol(a)
			#添加每个柱子的灰度背景
			rect(xleft = 0, ybottom = a[1, cols %% 2 == 1] - wd / 2, xright = par("usr")[2], ytop = a[2, cols %% 2 == 1] + wd / 2  , col = "gray91", border = NA)
			barplot(data, beside = T, horiz = T, width = wd, col = bar_cols , add = TRUE, axes = FALSE,xlim = c(0,max(data)))
			axis(side = 1, at = c(0, max(data)), labels = c(0, x_max))
		}
		title(xlab = "Mean proportion (%)")
		#添加标准误线条
		# x_coord = t(a)
		# mean = c(mean1,mean2)
		# stderr = c(df[,3],df[,5])

		# min_target_length = (max(data)) * 0.008
		# max_target_length = (max(data)) * 0.08
		# length = (stderr - min(stderr)) / (max(stderr) - min(stderr)) * (max_target_length - min_target_length) + min_target_length

		# arrows(y0 = x_coord, x0 = mean, y1 = x_coord, x1 = mean + length,   angle = 90, length = 0.05)


		dev.off()
	}

}


if(length(Groups) > 1){

    for(m in 1:(length(Groups)-1)){
                
        for (n in (m+1):length(Groups)){
        	#####创建两两差异分析的文件夹###########
			dir = paste(dir1, paste(Groups[m], "_vs_", Groups[n],sep = ""),sep = "/")
			system(paste("mkdir -p ", dir, sep = ""))
			setwd(dir)

			SampleInfo1 = rbind(SampleInfo[SampleInfo$Group == Groups[m],], SampleInfo[SampleInfo$Group == Groups[n],]) 
			diff_analysis(cogfile,SampleInfo1,"COG",m,n)
			diff_analysis(keggfile,SampleInfo1,"KEGG",m,n)
			
		}
    }
}








