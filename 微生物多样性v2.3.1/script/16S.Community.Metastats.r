#clean up 
rm(list = ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== multi_Metastats ===========
	Usage: 
	Rscript 16S.Community.Metastats.r taxonfile groupfile outputpath
		parameters ->
			taxonfile: [file -> always Community_structure/*.taxon.Abundance.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

taxonfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

setwd(outputpath)

######Data input###############
taxondata = read.table(taxonfile, header = T, sep = "\t", row.names = 1, check.names = F, quote = "", fill = T, comment.char = "")
sampleinfo = read.table(groupfile, header = F, sep = "\t")
colnames(sampleinfo) = c("sample", "Group")

######去除丰度表中非样本列
taxondata = taxondata[,as.character(sampleinfo$sample)]
loc = which(rownames(taxondata) == "All") 
if(length(loc) > 0) taxondata = taxondata[-loc,]

Groups = levels(as.factor(as.character(sampleinfo[,2])))
names = strsplit(strsplit(taxonfile,"/")[[1]][length(strsplit(taxonfile,"/")[[1]])],".",fixed=T)[[1]][1]

###两组metastats#####
if(length(Groups) > 1){

    for(m in 1:(length(Groups)-1)){
                
        for (n in (m+1):length(Groups)){
                
			#####创建两两差异分析的文件夹###########
			dir <- paste(outputpath, paste("Group", Groups[m], "_vs_Group", Groups[n], sep = ""), names, sep = "/")
			system(paste("mkdir -p ", dir, sep = ""))
			setwd(dir)
               
			#####创建新的分组信息###################           
			SampleInfo1 <- rbind(sampleinfo[sampleinfo$Group == Groups[m],], sampleinfo[sampleinfo$Group == Groups[n],])
            write.table(SampleInfo1, "sample.group", col.names = F, row.names = F, sep = "\t", quote = F)
					
			####构建差异分析的丰度表###############		
			count = taxondata[,as.character(SampleInfo1$sample)]
			count = count[rowSums(count) != 0,]
			countout = cbind(rownames(count), count)
			colnames(countout)[1] = names
			write.table(countout, paste(names, "_abundance.xls", sep = ""), row.names = F, sep = "\t", quote = F)
            
            # relative abundance
			if(nrow(count) == 1){
				data_relative_abundance = t(as.matrix(sapply(1:ncol(count),function(x) count[,x]/sum(count[,x]))))
			}else{
				data_relative_abundance = sapply(1:ncol(count),function(x) count[,x]/sum(count[,x]))
			}
			rownames(data_relative_abundance) = rownames(count)
			colnames(data_relative_abundance) = colnames(count)
			
            ####metastats analysis##################	
            system(paste("Rscript /home/zhengy/bin/modules/script/otu2shared.r", paste(names, "_abundance.xls", sep = ""), "sample.group", "0.03", ".", sep = " "))
            system("/home/panrf/Softwares/mothur/mothur '#metastats(shared = sample.shared, design = sample.group, processors = 2)'")

            res = read.table(list.files(pattern = "metastats"), skip = 5, header = T, sep = "\t")
            res$OTU = rownames(count)
            res$fdr = p.adjust(res$p.value, method = "BH")

			####set table index####################
			colnames(res) = c("Taxon", paste(Groups[n], "mean", sep="_"), paste(Groups[n], "variance", sep="_"), paste(Groups[n], "stderr", sep="_"), 
							paste(Groups[m], "mean", sep="_"), paste(Groups[m], "variance", sep="_"), paste(Groups[m], "stderr", sep="_"), "p_value", "q_value")
			
			#res = res[order(res$p_value),]

			###输出差异分析的所有结果############
            write.table(res, paste(names, "diff.all.xls", sep = "."), sep = "\t", row.names = F) 
            system("rm sample* mothur* *abundance*")

			##输出P<0.05差异分析结果##############
			res  = res[as.numeric(res[,2])>0.001 & as.numeric(res[,5])>0.001,]
            Sigres = subset(res, res[,8] < 0.05)
            Sigres =Sigres[which(Sigres[,1]!="No_Rank"),]
            Sigres =Sigres[which(Sigres[,1]!="Unassigned"),]
			# Boxplot
			if (nrow(Sigres) > 0){
			
				write.table(Sigres, paste(names, "diff.P0.05.xls", sep = "."), sep = "\t", row.names = F,quote = F)
				SigresForPlot = Sigres[which(Sigres[,1]!="Unassigned"),]

				DataForBoxplot = as.matrix(data_relative_abundance[rownames(data_relative_abundance) %in% SigresForPlot[,1], ] * 100)
				if(ncol(DataForBoxplot) == 1) DataForBoxplot = t(DataForBoxplot)				
				if(nrow(SigresForPlot) == 1) rownames(DataForBoxplot) = SigresForPlot[,1]

				# plot
				DataForBoxplot1 = data.frame(Taxon = rep(rownames(DataForBoxplot), ncol(DataForBoxplot)), Group = rep(as.character(SampleInfo1[,2]), each = nrow(DataForBoxplot)), Relative.Abundance = as.vector(as.matrix(DataForBoxplot)))

				if(length(levels(DataForBoxplot1$Taxon)) == 1){

					pdf(paste(names, ".taxon.DA.Boxplot.pdf", sep = ""),height = 6,width = 6)

				}else if(length(levels(DataForBoxplot1$Taxon)) == 2){

					pdf(paste(names, ".taxon.DA.Boxplot.pdf", sep = ""),height = 6,width = 12)
					par(mar = c(8,6,6,2), mfrow = c(1, 2),cex.axis = 1,cex.main = 1.2,cex.lab= 1.2)

				}else{

					n = ceiling(length(levels(DataForBoxplot1$Taxon)) / 3)
					width = 15
					height = 5*n
					pdf(paste(names, ".taxon.DA.Boxplot.pdf", sep = ""),height = height,width = width)
					par(oma = c(0,0,0,0),mar = c(8,6,6,2), mfrow = c(n, 3),cex.axis = 1.5,cex.main = 2,cex.lab= 2)

				}

				for (i in 1:length(SigresForPlot$Taxon)){

					tax = SigresForPlot$Taxon[i]
					new = DataForBoxplot1[which(DataForBoxplot1$Taxon %in% tax), ] 
					p_value = (SigresForPlot[which(SigresForPlot$Taxon %in% tax), ])$p_value

					if(p_value >= 0.01){p = "*"}
					if(p_value < 0.01 & p_value >= 0.001){p = "**"}
					if(p_value < 0.001 & p_value >= 0.0001){p = "***"}
					if(p_value < 0.0001){p = "****"} 

					boxplot(new$Relative.Abundance ~ factor(new$Group), new, 
						 	lwd = 2,
						 	xaxt = "n", 
						 	main = paste(tax,p), 
						 	ylab = "Relative.Abundance (%)",
						 	col = c("IndianRed1","DarkTurquoise"))
					
					text( x = c(1:2), y = (par("usr")[3]-(par("usr")[4] - par("usr")[3])*0.02), adj = c(1,1), srt = 45, cex = 1.2, labels = levels(factor(new$Group)), xpd = TRUE)
				}
				dev.off()
				
			}
			
        }
		
    }

}else{

    ########组数小于2组无法进行比较#####    
    stop("ERROR: *** Groups number less than 2 ! Please Cheack It! ***\n")

}
 
 
 
 
