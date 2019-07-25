# correlation 
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
    "============== do Sample cor ===========
    Usage: 
    Rscript cor.r countfile outputpath
        parameters ->
            countfile: [file -> always count.norm.xls];
            outputpath: [path -> path for output]; \n")
    options("show.error.messages" = F) 
    stop()

}

countfile  <- normalizePath(args[1])
outputpath <- normalizePath(args[2])
setwd(outputpath)

library(Hmisc)

Count <- read.table(countfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote = "", fill = T)

# =================================== COR ===================================
# 相关性分析
result <- rcorr(as.matrix(Count),type = "spearman")
Cor    <- result$r
corres <- Cor[lower.tri(Cor) == T]

# 热图调色
colfunc <- colorRampPalette(c("Gold","orange","Red3"))
color   <- colfunc(length(corres))

count <- 0
scater_panel <- function(x, y, ...){
	
	ll <- par("usr")
	rect(ll[1], ll[3], ll[2], ll[4])
	points(x, y, cex = 0.5, pch = 16, col = "red")
	abline(lm(x ~ y), col = "black")
	
}

hist_panel <- function(x, ...){

	usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(usr[1] ,usr[2]  , 0, 1.5))  
   	tax <- table(x)
    if(length(tax) < 11){

		breaks <- as.numeric(names(tax))
		y <- tax/max(tax)
		interbreak <- min(diff(breaks))*(length(tax)-1)/41
		rect(breaks-interbreak,0,breaks + interbreak,y)

    }else {
   
		h <- hist(x, plot = FALSE)
		breaks <- h$breaks; nB <- length(breaks)
		y <- h$counts; y <- y/max(y)
		rect(breaks[-nB], 0, breaks[-1], y)

    }
	tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
	if(class(tryd) != "try-error"){

		d$y <- d$y/max(d$y)
		lines(d)

   	}
}


heatmap_panel <- function(x, y, ...){

	count <<- count + 1
	ll <- par("usr")
	bg <- color[which(order(corres) == count)]
	rect(ll[1], ll[3], ll[2], ll[4], col = bg)
	text((ll[2] - ll[1]) / 2 + ll[1] , ll[3] + (ll[4] - ll[3]) / 2, label = sprintf("%0.4f", corres[count]))

}

png("Samples.pearson.corrletaion.png")
pairs(Count, 
	lower.panel = scater_panel,
	diag.panel = hist_panel,
	upper.panel = heatmap_panel,
	pch = 1, cex = 1,  
	cex.labels = 1,
	cex.axis = 1,
	gap = 0, 
	font.labels = NULL,
	col="Black")

dev.off()



