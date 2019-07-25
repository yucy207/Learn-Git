rm(list = ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
    "============== do betaNTI ===========
    Usage: 
        Rscript NTI.r groupfile tmppath outputpath
            parameters ->
                       
                     groupfile: [file -> always sample.groups];
                       tmppath: [path -> path for tmp];
                    outputpath: [path -> path for output]; \n")
    options("show.error.messages" = F) 
    stop()

}

groupfile = normalizePath(args[1])
tmppath = normalizePath(args[2])
outputpath = normalizePath(args[3])

library(ggplot2)

SampleInfo = read.table(groupfile, header = F, sep = "\t",quote="")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Groups = unique(SampleInfo$Group)

dataforplot = data.frame()
Statistics = data.frame()

for(i in 1:length(Groups)){
	newSample = as.character(SampleInfo$Sample[which(SampleInfo$Group %in% Groups[i])])
	NTI = paste(tmppath,paste(Groups[i],"NTI.xls",sep ="."), sep ="/")
	NTIdata = read.table(NTI, header = T, check.name = F, stringsAsFactors = F, quote = "")
	betaNTI = NTIdata[(nrow(NTIdata)-length(newSample) + 1):nrow(NTIdata),]
	rownames(betaNTI) = betaNTI[,1]
	betaNTI = betaNTI[,-1]
	#提取Beta_NTI指数，并且格式化成列表
	low = betaNTI[lower.tri(betaNTI)]
	result = data.frame(sample_i = rep(colnames(betaNTI)[-ncol(betaNTI)],times = c((ncol(betaNTI)-1):1)), 
						sample_j = unlist(lapply(2:nrow(betaNTI),function(t){rownames(betaNTI)[t:nrow(betaNTI)]})), 
						Beta_NTI = as.numeric(low),stringsAsFactors = F)

	result$Group = Groups[i]
	result$Median = median(result$Beta_NTI)

	num = length(result$Beta_NTI[result$Beta_NTI<2 & result$Beta_NTI>-2])
	result$Deterministic = 100 - num/nrow(result) *100
	result$Stochastic    = num/nrow(result) *100
	
	dataforplot = rbind(dataforplot,result)

	#统计
	result1 = data.frame(a = Groups[i], b = nrow(result) - num, c = num, d = median(result$Beta_NTI),
						 e = 100 - num/nrow(result) *100, f = num/nrow(result) *100)

	Statistics = rbind(Statistics,result1)

}
write.table(dataforplot, paste(outputpath,"betaNTI.xls",sep = "/"), sep = "\t", row.names = F, quote=F)
colnames(Statistics) = c("Group","|BetaNTI|>2","|BetaNTI|<2","Median","Deterministic","Stochastic")
write.table(Statistics, paste(outputpath,"betaNTI_statistics.xls",sep = "/"), sep = "\t", row.names = F, quote=F)

#plot
pdf(paste(outputpath,"betaNTI.pdf",sep = "/"), height = 6, width = 8)
p = ggplot(dataforplot,aes(x = dataforplot$Group,y = as.numeric(dataforplot$Beta_NTI))) + geom_point() +
		xlab("scale") + ylab("beta NTI")  +  scale_y_continuous(breaks = c(-2,0,2)) +
		theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
		theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face="bold")) +
		geom_errorbar(width=0.5, aes(y = dataforplot$Median, ymax = dataforplot$Median, ymin = dataforplot$Median), colour="blue") +
		geom_hline(aes(yintercept=-2), colour="black", linetype="dashed") + 
		geom_hline(aes(yintercept=2), colour="black", linetype="dashed") 

print(p)
dev.off()
