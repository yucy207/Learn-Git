rm(list = ls())
args = commandArgs(trailingOnly=T)
if(length(args) != 2){

	cat(
	"============== Use a variety of distance methods to do PCoA ===========
	Usage: 
	Rscript Enterotype.r taxofile outputpath
        parameters ->
            taxonfile: [file -> always *.taxon.Abundance.xls with row taxa and column samples];
            outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
}

taxonfile = normalizePath(args[1])
outputpath = normalizePath(args[2])
setwd(outputpath)

inputdata = read.table(taxonfile, header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T)    # 注意 quote="" 
loc = grep("OTUsize|Abundance", colnames(inputdata)) # 如果存在分类注释信息,则去除
if(length(loc) > 0) taxon = as.matrix(inputdata[,1:(loc-1)])
if( ncol(taxon) == 1){
	stop("Sorry: You must have at least 2 rows and 2 columns to plot heatmap!\n")
}else{
	taxon = taxon[rowSums(taxon) != 0, ]
}

index = which(rownames(taxon) == "All")       # 去除 All 行
if(length(index) > 0){
	taxon = as.matrix(taxon[-index, ])
	colnames(taxon) = colnames(inputdata)[1:(loc-1)]
}

if(nrow(taxon) > 1){
	data = sapply(1:ncol(taxon), function(t) taxon[,t]/sum(taxon[,t]))
	rownames(data) = rownames(taxon)
	colnames(data) = colnames(taxon)
}

# data=read.table("MetaHIT_SangerSamples.genus.txt", header=T, row.names=1, dec=".", sep="\t")
# data=data[-1,]

#############################
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        
	inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
}

data.dist=dist.JSD(data)

###################################
# library(cluster)  # pam() function
# pam(as.dist(data), k, diss=TRUE) # x is a distance matrix and k the number of clusters
require(cluster) # pam() function
pam.clustering = function(x,k) {
	cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
	return(cluster)
}

# data.cluster = pam.clustering(data.dist, k)

################################## evaluate the CH index for every number of clusters k: example
# require(clusterSim) # function index.G1()
# nclusters=NULL
	# for (k in 1:(ncol(data)-1)) { 
		# if (k==1) {
			# nclusters[k]=NA 
		# } else {
			# data.cluster_temp=pam.clustering(data.dist, k)
			# nclusters[k]=index.G1(t(data), data.cluster_temp, d = data.dist, centrotypes = "medoids")
		# }
	# }
## plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
# k = which.max(nclusters)
################################# This has shown that the optimal number of clusters for this particular dataset is k (k=3)

data.cluster=pam.clustering(data.dist, k = 3)

### Cluster validation
# obs.silhouette = mean(silhouette(data.cluster, data.dist)[,3]) # function silhouette() from the cluster

### Here is a function to remove the noise (removal low abundant genera): 
noise.removal <- function(dataframe, percent){
	Matrix <- dataframe
	bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
	Matrix_1 <- Matrix[bigones,]
	#print(percent)
	return(Matrix_1)
}
data.denoized = noise.removal(data, percent=0.00001)

library(ade4)
pdf("Enterotype.pdf")
par(mfrow=c(2,2))
################################### barplot
barplot(table(data.cluster), xlab="biotypes", ylab="Nb samples", col=as.numeric(levels(as.factor(as.numeric(data.cluster))))+1) #+1 to avoid black
box()
################################### PCA
obs.pca = dudi.pca(data.frame(t(data)), scannf=F, nf=3)
plot(obs.pca$li[,1], obs.pca$li[,2], main="PCA", pch=16, col=as.numeric(data.cluster)+1, xlab=paste("PC1"), ylab=paste("PC2"))
box()
################################### dPCoA
obs.dpcoa = dpcoa(data.frame(data.denoized), data.dist, scannf=F, nf=3)
plot(obs.dpcoa$dls[1:2], type="n", xlab=paste("PC1"), ylab=paste("PC2"), main="dPCoA")
s.class(obs.dpcoa$dls, fac=as.factor(data.cluster), cellipse=1,add.plot=TRUE)
points(obs.dpcoa$dls[1:2], col=as.numeric(data.cluster)+1, cex=1, pch=16)
################################### Between-class analysis (BCA)
obs.pca = dudi.pca(data.frame(t(data)), scannf=F, nf=3)
obs.bet = bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=3) #nf=k-1
plot(obs.bet$ls[,1],obs.bet$ls[,2], type="n", xlab=paste("PC1"), ylab=paste("PC2"), main="between class")
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cellipse=1, add.plot=TRUE) #plot the result using the s.class() function  sub="Between-class analysis",
points(obs.bet$ls[,1],obs.bet$ls[,2], col=as.numeric(data.cluster)+1, cex=1, pch=16)
dev.off()

# bet = cbind(obs.bet$ls,data.cluster)
# write.table(bet,"bet.xls",sep = "\t",col.names=NA)

#################################### PRINCIPAL COORDINATE ANALYSIS (PCoA)
# library(ggplot2)
# library(ggplotly)
# PC1 = obs.bet$ls[,1]
# PC2 = obs.bet$ls[,2]
# class = as.factor(data.cluster)

# ggplot(obs.bet$ls, aes(x = PC1, y = PC2, fill = class)) +
# geom_point() +
# stat_ellipse(geom = "polygon",alpha = 1/2)

# qplot(data = obs.bet$ls, x = PC1, y = PC2, colour = class)+
# stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = class))

# pdf("Enterotype_PCoA.pdf")
# obs.pcoa = dudi.pco(data.dist, scannf=F, nf=3)
# pcoa = cbind(obs.pcoa$li,data.cluster)
# write.table(pcoa,"pcoa.xls",sep = "\t",col.names=NA)
# plot(obs.pcoa$li[,1],obs.pcoa$li[,2], type="n", xlab=paste("Axis.1"), ylab=paste("Axis.2"), main="Principal coordiante")
# s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cellipse=1, add.plot=TRUE, col=c(2,3,4)) #plot the result using the s.class() function
### s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(2,3,4), sub="Principal coordiante analysis")
# points(obs.pcoa$li[,1],obs.pcoa$li[,2], col=as.numeric(data.cluster)+1, cex=1, pch=16)
# text(obs.pcoa$li[,1], obs.pcoa$li[,2], label=rownames(obs.pcoa$li), cex= 0.7, pos=3, xpd=TRUE)
# dev.off()

pdf("Enterotype_PCoA.pdf")
obs.pcoa = dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cellipse=1, sub="Principal coordiante analysis", possub ="topleft", col=c(2,3,4))
text(obs.pcoa$li[,1], obs.pcoa$li[,2], label=rownames(obs.pcoa$li), cex= 0.7, pos=3, xpd=TRUE)
dev.off()

cluster = data.frame(rownames(obs.pcoa$li),data.cluster)
colnames(cluster) = c("sample","cluster")
write.table(cluster,"cluster.xls",sep = "\t",col.names=NA)
