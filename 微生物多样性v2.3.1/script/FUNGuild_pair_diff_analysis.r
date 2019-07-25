rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== do Picrust different_analysis ===========
	Usage: 
	Rscript FUNGuild_pair_diff_analysis.r inputfile groupfile outputpath
		parameters ->
			  inputfile: [file -> always guilds.xls];
			  groupfile: [file -> always sample.groups];
		      outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

inputfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

setwd(outputpath)
dir1 = paste(outputpath, "pair_group_different_analysis",sep = "/")
system(paste("mkdir -p ", dir1, sep = ""))

SampleInfo = read.table(groupfile, header = F, sep = "\t",quote="")
SampleInfo[,1] = gsub("-","_",SampleInfo[,1])   # 全局替换
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]

Groups = levels(as.factor(as.character(SampleInfo[,2])))

t_test <- function(data.x, data.y){

	res <- lapply(1:(ncol(data.x)-1), function(t){
		data.1 <- unlist(data.x[,t])
		res <- t.test(data.1 ~ data.y, var.equal = FALSE)
		c(colnames(data.x)[t],res$estimate[1],res$estimate[2], res$conf.int[1],res$conf.int[2],mean(res$conf.int), res$p.value)
	})

	result <- data.frame(do.call(rbind, res))
	return(result)
}


# set color by group
cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
mycol= colors()[rep(cols,20)][1:length(Groups)]
group2col = data.frame(as.character(Groups),mycol,stringsAsFactors=FALSE)
colnames(group2col) = c("gro","col")

diff_analysis <-function(inputfile,SampleInfo1,m,n){

	data = read.table(inputfile, header = T, sep = "\t", row.names = 1,check.name = F, comment.char = "", quote = "", fill = T)
	colnames(data) = gsub("-","_",colnames(data))     # 全局替换
	
	data1 = data[,as.character(SampleInfo1$Sample)]
	
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

	colnames(result) = c("Guild",paste("mean",Groups[m], sep="_"),paste("mean",Groups[n], sep="_"),
						 "95.0% lower CI","95.0% upper CI","Difference between means","P.value")

	write.table(result,"diff.all.xls", row.names = F, sep = "\t", quote = F)

	#根据p值过滤
	p = apply(result[7], 1, as.numeric)
	df = result[( p<0.05 & p>0), ]

	if (nrow(df)>0){
		df <- df[order(as.numeric(as.character(df$P.value)), decreasing = TRUE), ]
		write.table(df,"diff.P0.05.xls", row.names = F, sep = "\t", quote = F)

		#################plot################
		#prepare data for plot
		#添加p值*号
		p_value = as.numeric(as.character(df$P.value))
		y <- numeric(length(p_value))
		y[p_value >= 0.01] = "*"
		y[p_value < 0.01 & p_value >= 0.001] = "**"
		y[p_value < 0.001 & p_value >= 0.0001] = "***"
		y[p_value < 0.0001] = "****"
		name = paste(df$Guild,y,sep = "")

		### mean 
		mean1 = apply(df[2], 1, as.numeric)
		mean2 = apply(df[3], 1, as.numeric)
		mean =  as.data.frame(rbind(mean1,mean2))
		data = as.matrix(mean)
		colnames(data) = name
		
		#### CI
		CI1 = apply(df[4], 1, as.numeric)
		CI2 = apply(df[5], 1, as.numeric)
		pvalue = apply(df[7], 1, as.numeric)
		#pvalue = format(as.numeric(pvalue), scientific = TRUE, digits = 3)
		for (i in 1:length(pvalue)){
			if (as.numeric(pvalue[i]) > 0.01){
				pvalue[i] = round(as.numeric(pvalue[i]), digits = 3)
			}else{
				pvalue[i] = format(as.numeric(pvalue[i]), scientific = TRUE, digits = 3)
			}
		}
		CI  = as.data.frame(rbind(CI1,CI2))
		data2 = as.matrix(CI)
		colnames(data2) = pvalue

		#bar_cols = c("SkyBlue3","SandyBrown")
		bar_cols = group2col$col[group2col$gro %in% c(Groups[m],Groups[n])]
		
		b = 1.5
		l = max(nchar(name))*0.08 # 最大标签字符数占用边距
		t = 0.35
		r = 1.5
		c = 0.2
		w = 2
		h = ncol(data2)*0.3
		
		height = 1+b+h+t
		width  = l + w + c + c + w + r
		
		pdf("diff.P0.05.pdf",width = width ,height = height)
		
		#layout
		layout(matrix(c(1,3,2,4),2,2,byrow=T),width = c(l + w + c,c + w + r), heights = c(1,height-1))

		###
		par(mai = c(0,0.5,0,0), xpd = TRUE)
		plot.new()
		legend("bottomleft", legend = c(Groups[m],Groups[n]), ncol = 2, fill = bar_cols, bty = "n")

		###
		#坐标轴单位换算
		x_unite = max(data)/w
		y_unite = ncol(data)/h
		par(mai = c(b, l, t, c), xaxs = "i", yaxs = "i", mgp = c(3,2,1), xpd = TRUE)
		wd = 0.8
		a = barplot(data, 
					beside = T,     # TRUE并列条形图,FALSE堆积条形图
					horiz = T,      # TRUE条形图垂直Y轴，FALSE条形图垂直X轴
					width = wd,     # 条形图宽度
					col = "white",  # 条形图颜色
					border = NA,    # 条形图边框
					axes = F,       # TRUE显示坐标轴
					axisnames = F   # TRUE显示默认条形图标签
					)
		#添加条形图标签
		text(x = par("usr")[1]-0.1*x_unite, y = (a[1,] + a[2,])/2, labels = colnames(data), adj = 1, srt = 0, xpd = TRUE )
		#添加每个柱子的灰度背景
		cols = 1:ncol(data)
		rect(xleft = 0, ybottom = a[1, cols %% 2 == 1] - wd / 2, xright = par("usr")[2]+0.05*x_unite, ytop = a[2, cols %% 2 == 1] + wd / 2  , col = "gray91", border = NA)
		#重新画堆积条形图
		barplot(data, beside = T, horiz = T, width = wd, col = bar_cols , axes = TRUE, axisnames = F, add = TRUE)
		#添加标题、坐标轴标签
		title(xlab = "Mean proportion (%)")

		###
		par(mai = c(0,0,0,0), xpd = TRUE)
		plot.new()
		legend("bottomleft", legend = c("95.0% confidence intervals"), border = NA, fill = "white", bty = "n")

		###
		#坐标轴单位换算
		x_unite = (max(data2)-min(data2))/w 
		y_unite = ncol(data2)/h
		par(mai = c(b, c, t, r), xaxs = "i", yaxs = "i", mgp = c(3,2,1), xpd = TRUE)
		wd = 0.8
		a = barplot(data2, beside = T, horiz = T, width = wd, col = "white", border = NA, axes = TRUE, axisnames = F)
		#添加每个柱子的标签
		text(x = par("usr")[2]+0.1*x_unite, y = (a[1,] + a[2,])/2, labels = colnames(data2), adj = 0, srt = 0, xpd = TRUE )
		#添加奇数柱子的灰度背景
		cols = 1:ncol(data2)
		rect(xleft = par("usr")[1]-0.05*x_unite, ybottom = a[1, cols %% 2 == 1] - wd / 2, xright = par("usr")[2]+0.05*x_unite, ytop = a[2, cols %% 2 == 1] + wd / 2, col = "gray91", border = NA)
		#添加置信水平
		abline(v = 0, lty = 5, lwd = 1, xpd = FALSE)  #在x=0处有一条垂直虚线
		#添加线
		arrows(x0 = data2[1,], y0 = (a[1,] + a[2,])/2, x1 = data2[2,], y1 = (a[1,] + a[2,])/2, code = 3, angle = 90, length = 0.05, lty = 1, lwd = 2, xpd = TRUE)
		#添加点
		
		
		for (i in 1:length(a[1,])){

			if (data2[1,][i] < 0){
				if (data2[1,][i] < data2[2,][i]){
					legend(x = data2[2,][i] + (data2[1,][i] - data2[2,][i])/2 - 0.125*x_unite, y = ((a[1,] + a[2,])/2)[i], y.intersp = -1, legend = NA, pch = 19, pt.cex = 1.5, col =  bar_cols[2], bty = "n")
				}else{
					legend(x = data2[1,][i] + (data2[2,][i] - data2[1,][i])/2 - 0.125*x_unite, y = ((a[1,] + a[2,])/2)[i], y.intersp = -1, legend = NA, pch = 19, pt.cex = 1.5, col =  bar_cols[2], bty = "n")
				}
			}else{
				if (data2[1,][i] > data2[2,][i]){	
					legend(x = data2[2,][i] + (data2[1,][i] - data2[2,][i])/2 - 0.125*x_unite, y = ((a[1,] + a[2,])/2)[i], y.intersp = -1, legend = NA, pch = 19, pt.cex = 1.5, col =  bar_cols[1], bty = "n")
				}else{	
					legend(x = data2[1,][i] + (data2[2,][i] - data2[1,][i])/2 - 0.125*x_unite, y = ((a[1,] + a[2,])/2)[i], y.intersp = -1, legend = NA, pch = 19, pt.cex = 1.5, col =  bar_cols[1], bty = "n")
				}
			}
		}
		#添加标题、坐标轴标签
		text(x = par("usr")[2]+1*x_unite, y = par("usr")[4]/2, labels = "p-value", adj = 0, srt = 90, xpd = TRUE)
		title(xlab = "Difference in mean proportion(%)")
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
			
			diff_analysis(inputfile,SampleInfo1,m,n)
		}
    }
}








