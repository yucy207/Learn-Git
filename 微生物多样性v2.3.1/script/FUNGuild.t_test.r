rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== do FUNGuild different_analysis ===========
	Usage: 
	Rscript FUNGuild.t_test.r guildfile groupfile outputpath
		parameters ->
			  guildfile: [file -> always guilds.xls];
			  groupfile: [file -> always sample.groups];
		      outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

guildfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

setwd(outputpath)
dir1 = paste(outputpath, "pair_group_different_analysis",sep = "/")
system(paste("mkdir -p ", dir1, sep = ""))

SampleInfo = read.table(groupfile, header = F, sep = "\t",quote="")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Groups = unique(as.character(SampleInfo [,2]))

t_test <- function(data.x, data.y){
	res <- lapply(1:(ncol(data.x)-1), function(t){
  		data.1 <- unlist(data.x[,t])
  		res <- t.test(data.1 ~ data.y, var.equal = TRUE)
   		c(colnames(data.x)[t],res$estimate[1],res$estimate[2], res$conf.int[1],res$conf.int[2],mean(res$conf.int), res$p.value)
 	})

 	result <- data.frame(do.call(rbind, res))
 	return(result)
}


diff_analysis <-function(inputfile,SampleInfo1,m,n){

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

	write.table(result, "diff.all.xls", row.names = F, sep = "\t", quote = F)

	#根据p值过滤
	p = apply(result[7], 1, as.numeric)
	df = result[( p<0.05 & p>0), ]

	if (nrow(df)>0){

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
		name = paste(df$Category,y,sep = "")

		mean1 = apply(df[2], 1, as.numeric)
		mean2 = apply(df[3], 1, as.numeric)
		mean =  as.data.frame(rbind(mean1,mean2))
		colnames(mean) = name
		data = as.matrix(mean)

		bar_cols = c("SkyBlue3","SandyBrown")
		height_barplot = 0.15 * nrow(df) * 2 + 0.15 * (nrow(df) - 1) + 1 + 0.3
		height =  height_barplot + 1

		#layout
		width = 12 
		pdf("diff.P0.05.pdf",width = width ,height = height)

		layout(matrix(c(1,2), nrow = 2, byrow = T),heights = c(1,height_barplot))
		par(mai = c(0.1,3,0.2,3), xaxs = "i", yaxs = "i")
		plot(1:10, type = "n", yaxt = "n", xlab = "", xaxt = "n", ylab = "", bty = "n")

		rect(xleft = 2, ybottom = 4, xright = 2.8, ytop = 6, col = bar_cols[1],lwd = 1.5)
		rect(xleft = 5, ybottom = 4, xright = 5.8, ytop = 6, col = bar_cols[2],lwd = 1.5)
		text(x = 3*1.02, y = 5, adj = 0, labels = Groups[m], cex = 1.3)
		text(x = 6*1.02, y = 5, adj = 0, labels = Groups[n], cex = 1.3)
		
		x = 0.1* max(nchar(colnames(data)))
		if (x > 4){
			xx = x
		}else{
			xx = 4
		}		
		par(las = 1, mai = c(1, xx, 0.3, 0.5),  xaxs = "i", yaxs = "i", mgp = c(3, 2, 1), lend = 2,lwd = 1.5)

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
			diff_analysis(guildfile,SampleInfo1,m,n)
		}
    }
}








