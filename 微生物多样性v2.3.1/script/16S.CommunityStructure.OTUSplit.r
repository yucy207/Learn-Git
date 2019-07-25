# OTU.tax.0.03 ����ͬ����ˮƽ�����Ÿ�Ŀ�����֣����в��

rm(list=ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
	"Usage: 
	Rscript 16S.CommunityStructure.OTUSplit.r otufile outputpath
		parameters ->
			 otufile: [file -> always otu.tax.0.03.xls];
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

otufile = normalizePath(args[1])
outputpath = normalizePath(args[2])

# otufile = normalizePath("/home/zhengxq/workdir/meta_genomics/16B0930A/report/Group1/subsample_otu.tax.0.03.xls")
# outputpath = normalizePath("/home/zhengxq/workdir/meta_genomics/16B0930A/report/Group1/Community/Community_Structure/test")

Data = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
Data[,(ncol(Data)-6):ncol(Data)] = sapply( (ncol(Data)-6):ncol(Data), function(x) gsub("\"","",Data[,x]) )

Samples = colnames(Data)[1:(ncol(Data)-9)]
N = length(Samples)

classidata = data.frame(as.matrix(Data), check.names=FALSE, stringsAsFactors = FALSE)
classidata[classidata == ""] = "No_Rank"
colnames(classidata) = colnames(Data)

# ����ע��ˮƽ
tax = colnames(Data)[(ncol(Data)-6):ncol(Data)]    

# ���㵥��ˮƽ��,�������������еķ��
taxcount = function(t){

	# OTU�͸�ˮƽ������Ϣ
	A = classidata[,c(Samples,t)]     
	classi = unique(A[, ncol(A)])
	
	# ����ͬ����OTU���
	result = as.matrix(sapply(classi, function(c) apply(as.matrix(A[ A[, ncol(A)] %in% c, Samples]), 2, function(x) sum(as.numeric(x)))))
	
	if (nrow(result) != 1){
	
		if(rownames(result)[1] == colnames(Data)[1]){		# ֻ��һ�����֣���ת��
	
			result = t(result)
		
		}
		
	}
	
	colnames(result) = Samples
	
	# ��Ӹ�ˮƽ���ϲ�ע��
	Anno = sapply(rownames(result), function(r) { r = gsub("\\(", "\\\\(", r);r = gsub("\\)", "\\\\)", r);r = gsub("\\[", "\\\\[", r);r = gsub("\\-", "\\\\-", r);r = gsub("\\]", "\\\\]", r);r = gsub("\\+", "\\\\+", r);r = gsub("\\^", "\\\\^", r);classidata[grep(paste("^", r, "$", sep = ""), A[, ncol(A)])[1], tax[1:grep(t, tax)]]})
	
	if(is.matrix(Anno)){
		da = matrix(unlist(Anno), nrow = nrow(result), byrow = TRUE)
		rownames(da) = colnames(Anno)
		colnames(da) = rownames(Anno)
		da[which(rownames(da) == "No_Rank"), -ncol(da)] = "-"
		return(cbind(result, da))
	}else{
		return(cbind(result, as.matrix(Anno)))
	}
	
}

# ��������ˮƽ�£��������������еķ��
taxdata = lapply(tax, taxcount)

# output OTU size and relative abundance by tax
for (i in 2:length(tax)){

	out = as.matrix(taxdata[[i]][, Samples])
	if(ncol(out) == 1) next
	
	temp1 = cbind(out, apply(out, 1, function(o) sum(as.numeric(o))))      # Add sum column
	#temp2 = cbind(temp1, taxdata[[i]][,(N+1):ncol(taxdata[[i]])])      #  add annotation
	if(nrow(out) == 1){
		taxon = t(as.matrix(taxdata[[i]][,(N+1):ncol(taxdata[[i]])]))
	}else{
		taxon = as.matrix(taxdata[[i]][,(N+1):ncol(taxdata[[i]])])
	}
	temp2 = cbind(temp1, taxon)
	temp3 = temp2[order(as.numeric(temp2[,N+1]), decreasing = T), ]       # �� Abundance��������
	if(is.null(nrow(temp3))){
		temp3 = t(as.matrix(temp3))
		rownames(temp3) = rownames(taxdata[[i]])
	}
	colnames(temp3)[N+1] = "Abundance"
	
	# OTU raw size	
	rawdata = temp3[,c(ncol(temp3),1:(ncol(temp3)-1))]        # ���һ�з��ڵ�һ��
	if(is.null(nrow(rawdata))) rawdata = t(as.matrix(rawdata))
	write.table(rawdata, paste(outputpath, paste(tax[i], ".taxon", ".Abundance.xls", sep = ""),sep = "/"), row.names = F, sep = "\t", quote = F)
	rowsums = c("All", apply(temp1, 2, function(t) sum(as.numeric(t))))        # ���һ�У��Ӻ���Ϣ
	write.table(t(as.matrix(rowsums)), paste(outputpath, paste(tax[i], ".taxon", ".Abundance.xls", sep = ""),sep = "/"), row.names = F, col.names = F, sep = "\t", quote = F, append = T)
	
	# relative abundance
	ra = sapply(Samples, function(x) as.numeric(temp3[,x])/as.numeric(rowsums[x]))
	radata = cbind(ra, temp3[,(N+2):ncol(temp3)])        # add annotation
	radata = radata[,c(ncol(radata),1:(ncol(radata)-1))]      # ���һ�з��ڵ�һ��
	write.table(radata, paste(outputpath, paste(tax[i], ".taxon", ".RelativeAbundance.xls", sep = ""),sep = "/"), row.names = F, sep = "\t", quote = F)
	
}


