# Differential Analysis Using ANOVA
rm(list=ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
    "Usage: 
    Rscript Anova.r countfile groupfile Threshold outputpath
        parameters ->
            countfile: [file -> always transcript_count_matrix.fmt.csv with row gene and column samples];
            groupfile: [file -> always sample.groups without header and only two column];
            outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()
	
}

countfile  = normalizePath(args[1])
groupfile  = normalizePath(args[2])
outputpath = normalizePath(args[3])

library(DESeq2)

countdata = read.table(countfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote = "", fill = T,stringsAsFactors = F)

SampleInfo = read.table(groupfile, header = F, sep="\t")
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("Sample", "Group")

# sort by groups
SampleInfo = SampleInfo[order(SampleInfo$Group),]
Groups = levels(SampleInfo$Group)

# prepare data
count = countdata[,rownames(SampleInfo)]
count = count[rowSums(count) != 0 , ]
#re.count = apply(count, 2, function(d) d/sum(as.numeric(d)))
#筛选相对丰度大于0.001的gene
# index = which(rowSums(re.count) > 0.001)
# DataForDEA = re.count[index, ]

#标准化
dds <- DESeqDataSetFromMatrix(countData = count, colData = SampleInfo, design = ~ Group)
dds <- DESeq(dds,betaPrior=FALSE)
cnt <- as.data.frame(counts(dds, normalized=TRUE))

DataForDEA = as.matrix(cnt)
	
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
	
result = data.frame(Transcript_id = rownames(DataForDEA), 
					Gene = countdata$gene_name[which(rownames(countdata) %in% rownames(DataForDEA))],
					DataForDEA, means, res, FDR = fdr ,Type = 0, check.names = F)

result = result[order(result$Pvalue),]
p_index = which(result$Pvalue < 0.05)
result[1:length(p_index),]$Type = "DEG";
result[(length(p_index) + 1):nrow(result),]$Type = "Not DEG";

anova_file = paste(outputpath, "/transcript_ANOVA_Test_Result.xls", sep="")
write.table(result, anova_file, row.names = F, quote = F,sep = "\t")




