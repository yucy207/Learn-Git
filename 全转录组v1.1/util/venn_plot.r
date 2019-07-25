args = commandArgs(trailingOnly=T)


len = length(args)

data <-  lapply(args[1:(len - 1)], function(t) { scan(t,  what = character())  })
names <- lapply(args[1:(len - 1)], function(t) { b <- unlist(strsplit(t, "/")); b[length(b)-2] })

names(data) <- names

common_gene <- Reduce(intersect, data)

output <- args[len]

gene_file <- paste(output, "/common_gene.tsv", sep = "")
write.table(common_gene,  gene_file, sep = "\t", quote = F, row.names = FALSE, col.names = FALSE )



library(VennDiagram)

mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,84,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)
mycol <- colors()[rep(mycol,20)]



venn.diagram.new = function (x, filename, height = 10, width = 10, resolution = 500, fill = mycol, # col = mycol 
    imagetype = "pdf", units = "in", compression = "lzw",
    main = NULL, sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
    main.fontfamily = "serif", main.col = "black", main.cex = 3, 
    main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
    sub.fontfamily = "serif", sub.col = "black", sub.cex = 3, 
    sub.just = c(0.5, 1), category.names = names(x),
    hyper.test = FALSE, total.population = NULL, lower.tail = TRUE) 
{

    for (i in 1:length(x)) {
        x[[i]] <- unique(x[[i]])
		x[[i]] <- x[[i]][!is.na(x[[i]])]
    }
        
    if (0 == length(x) | length(x) > 5) {
        stop("Incorrect number of groups.", call. = FALSE)
    }
    
	if (1 == length(x)) {
        list.names <- category.names
        if (is.null(list.names)) {
            list.names <- ""
        }
        grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), fill = fill[1],
            category = list.names, ind = FALSE, cat.dist)
    }

    if (2 == length(x)) {
        grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), area2 = length(x[[2]]), 
			cross.area = length(intersect(x[[1]], x[[2]])), fill = fill[1:2], category = category.names, ind = FALSE, cat.cex = 2, cex = 2)
    }
    
	if (3 == length(x)) {   # A B C
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        nab <- intersect(A, B)
        nbc <- intersect(B, C)
        nac <- intersect(A, C)
        nabc <- intersect(nab, C)
        grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), area2 = length(B), area3 = length(C), n12 = length(nab), 
                n23 = length(nbc), n13 = length(nac), n123 = length(nabc), fill = fill[1:3], category = category.names, ind = FALSE,cat.cex = 2, cex = 1.5)  # col = col[1:3]
    }
    
	if (4 == length(x)) {    # A C D B
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        D <- x[[4]]
        n12 <- intersect(A, B)
        n13 <- intersect(A, C)
        n14 <- intersect(A, D)
        n23 <- intersect(B, C)
        n24 <- intersect(B, D)
        n34 <- intersect(C, D)
        n123 <- intersect(n12, C)
        n124 <- intersect(n12, D)
        n134 <- intersect(n13, D)
        n234 <- intersect(n23, D)
        n1234 <- intersect(n123, D)
        grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), area2 = length(B), area3 = length(C), area4 = length(D), 
                n12 = length(n12), n13 = length(n13), n14 = length(n14), n23 = length(n23), n24 = length(n24), n34 = length(n34), 
                n123 = length(n123), n124 = length(n124), n134 = length(n134), n234 = length(n234), n1234 = length(n1234), 
				fill = fill[c(1,3,4,2)], category = category.names, ind = FALSE, cat.cex = 2, cex = 1.5)  # col = col[1:4]
    }
		
	if (5 == length(x)) {     # A C E D B
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        D <- x[[4]]
        E <- x[[5]]
        n12 <- intersect(A, B)
        n13 <- intersect(A, C)
        n14 <- intersect(A, D)
        n15 <- intersect(A, E)
        n23 <- intersect(B, C)
        n24 <- intersect(B, D)
        n25 <- intersect(B, E)
        n34 <- intersect(C, D)
        n35 <- intersect(C, E)
        n45 <- intersect(D, E)
        n123 <- intersect(n12, C)
		n124 <- intersect(n12, D)
        n125 <- intersect(n12, E)
        n134 <- intersect(n13, D)
        n135 <- intersect(n13, E)
        n145 <- intersect(n14, E)
        n234 <- intersect(n23, D)
        n235 <- intersect(n23, E)
        n245 <- intersect(n24, E)
        n345 <- intersect(n34, E)
		n1234 <- intersect(n123, D)
        n1235 <- intersect(n123, E)
        n1245 <- intersect(n124, E)
        n1345 <- intersect(n134, E)
        n2345 <- intersect(n234, E)
        n12345 <- intersect(n1234, E)
        grob.list <- VennDiagram::draw.quintuple.venn(area1 = length(A), area2 = length(B), area3 = length(C), area4 = length(D), 
                area5 = length(E), n12 = length(n12), n13 = length(n13), n14 = length(n14), n15 = length(n15), n23 = length(n23), 
                n24 = length(n24), n25 = length(n25), n34 = length(n34), n35 = length(n35), n45 = length(n45), n123 = length(n123), 
                n124 = length(n124), n125 = length(n125), n134 = length(n134), n135 = length(n135), n145 = length(n145), n234 = length(n234), 
                n235 = length(n235), n245 = length(n245), n345 = length(n345), n1234 = length(n1234), n1235 = length(n1235), 
                n1245 = length(n1245), n1345 = length(n1345), n2345 = length(n2345), n12345 = length(n12345), fill = fill[c(1,3,5,4,2)],  # col = col[1:5]
                category = category.names, ind = FALSE,cat.cex = 2, cex = 1.5) 
    }

    if (length(x) == 2 & !is.null(total.population) & hyper.test) {
        val.p = calculate.overlap.and.pvalue(x[[1]], x[[2]], 
            total.population, lower.tail = lower.tail)
        if (is.null(sub)) {
            sub = paste0("p = ", signif(val.p[3], digits = 2))
        }
        else {
            sub = paste0(sub, ", p = ", signif(val.p[3], digits = 2))
        }
    }
    if (!is.null(sub)) {
        grob.list <- add.title(gList = grob.list, x = sub, pos = sub.pos, fontface = sub.fontface, fontfamily = sub.fontfamily, col = sub.col, cex = sub.cex)
    }
    if (!is.null(main)) {
        grob.list <- add.title(gList = grob.list, x = main, pos = main.pos, fontface = main.fontface, fontfamily = main.fontfamily, col = main.col, cex = main.cex)
    }
	
    if (!is.null(filename)) {
        current.type <- getOption("bitmapType")
        if (length(grep("Darwin", Sys.info()["sysname"]))) {
            options(bitmapType = "quartz")
        }else {
            options(bitmapType = "cairo")
        }
	}
    
	if ("tiff" == imagetype) {
        tiff(filename = filename, height = height, width = width, units = units, res = resolution, compression = compression)
    }
    else if ("png" == imagetype) {
        png(filename = filename, height = height, width = width, units = units, res = resolution)
    }
    else if ("svg" == imagetype) {
        svg(filename = filename, height = height, width = width)
    }
	else if ("pdf" == imagetype) {
        pdf(file = filename, height = height, width = width)
    }
    else {
        stop("You have misspelled your 'imagetype', please try again")
    }
	
	# plot
    grid.draw(grob.list)
	
	# prepare text
	# all = paste("Unique objects: All = ", length(unique(unlist(DataForPlot))), sep = "")
	# other = sapply(1:length(DataForPlot), function(d) {len = length(unlist(DataForPlot[[d]])); paste(names(DataForPlot)[d], len, sep = " = ")})
	# str = paste(all, paste0(other, collapse = ";"), sep= ";")
	
	# add text
	#grid.text(label = str, vjust = 37.5)
    dev.off()
	
    options(bitmapType = current.type)
    return(1)

}

pdf_file <- paste(output, "/venn.pdf", sep = "")
venn.diagram.new(data, filename = pdf_file)





