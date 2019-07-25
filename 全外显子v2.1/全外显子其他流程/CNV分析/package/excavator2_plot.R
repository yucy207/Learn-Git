args           <- commandArgs(TRUE)
FileSeg        <- args[1]
FileCall       <- args[2]
PathOut        <- args[3]

y <- read.table(FileSeg,sep="\t",quote="\"",fill=T,header=T) # # HSLMResults_

log2RSeq <- as.numeric(as.character(y[,5]))
chrSeq <- as.character(y[,1])
SegSeq <- as.numeric(as.character(y[,6]))
PositionSeq <- as.numeric(as.character(y[,2]))

z <- read.table(FileCall,sep="\t",quote="\"",fill=T,header=T) # # FastCallResults_

chrCall <- as.character(z[,1])
StartCall <- as.numeric(as.character(z[,2]))
EndCall <- as.numeric(as.character(z[,3]))
Call <- as.numeric(as.character(z[,5]))

UniqueChr <- unique(chrCall) # # call的unique染色体

FilePlot <- file.path(PathOut, "excavator2_PlotResultsRegion.pdf")
pdf(FilePlot,height=5,width=15)
par(mfrow=c(1,1))

for (i in 1:length(UniqueChr))
{		
	indCall <- which(chrCall==UniqueChr[i]) # # 含有CNV的染色体
	if (length(indCall)!=0)
	{
		StartCallC <- StartCall[indCall]
		EndCallC <- EndCall[indCall]
		spareGap <- EndCallC - StartCallC	# # 两侧扩展长度为区域等长
		for (j in 1:(length(indCall)))
		{   
		    minPos <- StartCallC[j] - spareGap[j]
			maxPos <- EndCallC[j] + spareGap[j]
			if(minPos < 0) {minPos <- 0}
		    # # HSLMResults_ 中某一染色体的区域
		    indSeqRegion <- which((chrSeq == UniqueChr[i]) & (PositionSeq >= minPos) & (PositionSeq <= maxPos)) # # 两侧扩展
			
			PositionSeqC <- PositionSeq[indSeqRegion]   # # Position
	        log2RSeqC    <- log2RSeq[indSeqRegion]      # # log2
	        SegSeqC      <- SegSeq[indSeqRegion]        # # seg
			plotTitle    <- paste("chr",UniqueChr[i],":",StartCallC[j],"-",EndCallC[j])
			
	        plot(PositionSeqC,log2RSeqC,ylim=c(-3,3),xlim=c(minPos,maxPos),main=plotTitle,pch=16,cex=0.5,xlab="Position",ylab="log2ratio",col="Blue")
	        lines(PositionSeqC, SegSeqC,lwd=2,col="red")
	        abline(h=0,lty=2,lwd=1,col="black")
		}
	}
}
dev.off()