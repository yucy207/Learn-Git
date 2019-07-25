# Find the best subset of environmental variables, so that the Euclidean distances of scaled environmental variables have the maximum (rank) correlation with community dissimilarities.
# �õ���Ⱥ����������Ļ����������

rm(list=ls())

args = commandArgs(trailingOnly=T)
if(length(args) != 3){
	
	cat(
	" ============== Find the best subset of environmental variables =========
	Usage: 
	Rscript 16S.EnvironmentFactors.BioEnv.r envfile otufile outputpath
		parameters ->
			envfile: [file -> always env.xls with row env.xls and colum environmental factors];
			otufile: [file -> always otu.tax.0.03.xls];
		     outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

library(vegan)

envfile = normalizePath(args[1])
otufile = normalizePath(args[2])
outputpath = normalizePath(args[3])

envdata = read.table(envfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
otudata = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
otudata = otudata[, rownames(envdata)]

otudata = otudata[rowSums(otudata) != 0,]

bioenv = bioenv(t(otudata), envdata, method = "pearson")    # ���öԻ����������ݽ��б�׼���Ĳ���

a = summary(bioenv)  
res =  data.frame(EnvFactors = a$variables, Size = a$size, Correlation = a$correlation)
output = paste(outputpath, "SubEnvsCor.xls", sep = "/")
write.table(res, output, sep = "\t", row.names = F, quote = F)
#capture.output(summary(bioenv), file = paste(outputpath, "SubEnvsCor.csv", sep = "/"))    # print to file

# ================================ BioEnv vs. Mantel Test ==========================
# If you want to study the 'significance' of bioenv results, you can use function mantel or mantel.partial which use the same definition of correlation. 
# However, bioenv standardizes environmental variables to unit standard deviation using function scale and you must do the same in mantel for comparable results. 
# ����ʹ����ͬ������Է�����������������BioEnv���öԻ������������ݽ����˱�׼������Mantel��������Ҫʹ��scale()

