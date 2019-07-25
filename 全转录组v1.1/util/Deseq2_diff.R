#!/usr/bin/env Rscript

library(docopt)
"Usage: Deseq2_diff.R [options] INPUT CASE_GROUP CONTROL_GROUP CASE CONTROL LOG2FC OUTPUT
Options:
   -w --width=width    the width of viewport  [default: 12]
   -h --height=width   the height of viewport [default: 12]
Arguments:
  INPUT          the input file name
  CASE_GROUP     the case group name
  CONTROL_GROUP  the control group name
  CASE           the case group name
  CONTROL        the control group name
  LOG2FC         the log2 fold change
  OUTPUT         the output directory name" -> doc


opts          <- docopt(doc)
input         <- opts$INPUT
case_group    <- opts$CASE_GROUP
control_group <- opts$CONTROL_GROUP
case          <- opts$CASE
control       <- opts$CONTROL
log2fc        <- as.numeric(opts$LOG2FC)
output        <- opts$OUTPUT



library(DESeq2)

x <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)

control_count           <- sapply(unlist(strsplit(control, split=",")), function(t){x[[t]] })
colnames(control_count) <- unlist(strsplit(control, split=","))

case_count              <- sapply(unlist(strsplit(case, split=",")), function(t){ x[[t]] })
colnames(case_count)    <- unlist(strsplit(case, split=","))

count                   <- cbind(data.frame(control_count, check.names=F), data.frame(case_count, check.names=F))
rownames(count)         <- rownames(x)


countData <- count[apply(count, 1, sum) > 0 , ]



colData   <- data.frame(row.names = colnames(countData), 
	                    condition = rep(
	                    	         	c(control_group, case_group),
	                    	         	times  = c(ncol(control_count), ncol(case_count)) ,
	                    	         	levels = c(control_group, case_group)
	                    	         ) 	        
	                    )
colData$condition <- relevel(colData$condition, ref = control_group)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds <- DESeq(dds,betaPrior=FALSE)
cnt <- as.data.frame(counts(dds, normalized=TRUE))
sample_num <- ncol(cnt)
cnt$baseMeanA <- apply( sapply(unlist(strsplit(control, split=",")), function(t){cnt[[t]]}), 1, mean )
cnt$baseMeanB <- apply( sapply(unlist(strsplit(case, split=",")), function(t){cnt[[t]]}), 1, mean )
res <- as.data.frame(results(dds))
res <- cbind(cnt, res)
res$type <- "Not DEG"
res$type[res$pvalue < 0.05 & res$log2FoldChange >= log2fc ] <- "Up"
res$type[res$pvalue < 0.05 & res$log2FoldChange <= -(log2fc)] <- "Down"
res$type <- factor(res$type, levels = c("Up", "Down", "Not DEG"))
diff_file <- paste(output, "/diff.xls", sep="")
write.table(res, diff_file, sep="\t", quote=F,  col.names = NA)


# correlation plot
library(ggplot2)

correlation_file <- paste(output, "/correlation.pdf", sep="")
pdf(correlation_file)
ggplot(res, aes(x=log10(baseMeanA+ 0.00000000001),y= log10(baseMeanB+0.00000000001))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color="red") + 
  scale_x_continuous(limits=c(-2.5, 5)) + 
  scale_y_continuous(limits=c(-2.5, 5)) + 
  labs(x="log10(FPKM1)", y="log10(FPKM2)")
dev.off()

# MA plot
ma_file <- paste(output, "/MA.pdf", sep="")
pdf(ma_file)
ggplot(res, aes(x = baseMean, y = log2FoldChange , colour = type)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 2e+04)) + 
  scale_y_continuous(limits = c(-20, 20)) + 
  theme(legend.title = element_blank()) + 
  labs(x = "A", y="log2(FC)", tilte="MA plot")
dev.off()

# valcano plot
valcano_file <- paste(output, "/valcano.pdf", sep="")
pdf(valcano_file)
ggplot(res, aes(x = log2FoldChange, y = -log10(res$pvalue), color =type)) + 
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = bquote(paste(log[2],"(fold change)",sep="")), y = bquote(paste(-log[10],"(p value)",sep="")), title = "genes: CASE/CONTROL") +
  scale_x_continuous(limits=c(-10,10)) 
dev.off()


# heatmap plot

library(gplots)
myheatcol   = colorpanel(75, 'green','black','red')
data <- data.matrix(res[res$type != "Not DEG", 1:sample_num])
#res <- res[res$type != 'Not DEG', ]
#if( nrow(res) < 50 ){
#		data <- data.matrix(res[ , 1:sample_num])	   # 小于50,全画
#	}else{
#		res <- res[order(res$pvalue), ]
#		data <- data.matrix(res[1:50, 1:sample_num])           #  相对丰度最高的50个
#	}
heatmap_file <- paste(output, "/heatmap.pdf", sep="")
pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )
heatmap.2(
	  data,
	dendrogram = 'row', 
      Colv         = FALSE,
      col          = myheatcol, 
      scale="row", 
      margins      = c(8,10),
      density.info = "none",
      trace="none",
      key=TRUE,
      keysize=0.6,
      cexRow=0.01,
	#cexRow=0.8,
      cexCol=0.8,
     srtCol=90
)

# Top50 heatmap plot

library(gplots)
myheatcol   = colorpanel(75, 'green','black','red')
#data <- data.matrix(res[res$type != "Not DEG", 1:sample_num])
res <- res[res$type != 'Not DEG', ]
if( nrow(res) < 50 ){
               data <- data.matrix(res[ , 1:sample_num])          # 小于50,全画
	}else{
               res <- res[order(res$pvalue), ]
               data <- data.matrix(res[1:50, 1:sample_num])           #  相对丰度最高的50个
	}
heatmap_file <- paste(output, "/Top50.heatmap.pdf", sep="")
pdf(heatmap_file , width=as.numeric(opts$w), height = as.numeric(opts$h) )
heatmap.2(
          data,
        dendrogram = 'both',
      #Colv         = FALSE,
      col          = myheatcol,
      scale="row",
      margins	   = c(8,10),
      density.info = "none",
      trace="none",
      key=TRUE,
      keysize=0.6,
	cexRow=0.8,
      cexCol=0.8,
     srtCol=90
)
