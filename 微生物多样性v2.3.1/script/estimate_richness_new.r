estimate_richness_new = function (OTU, measures = NULL) {

    renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
    names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", "simpson", "invsimpson", "fisher")
    if (is.null(measures)) {
        measures = as.character(renamevec)
    }
    if (any(measures %in% names(renamevec))) {
        measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% measures]
    }
    if (!any(measures %in% renamevec)) {
        stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
    }
    outlist = vector("list")
    estimRmeas = c("Chao1", "Observed", "ACE")
    if (any(estimRmeas %in% measures)) {
        outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
    }
    if ("Shannon" %in% measures) {
        outlist <- c(outlist, list(shannon = diversity(OTU, index = "shannon")))
    }
    if ("Simpson" %in% measures) {
        outlist <- c(outlist, list(simpson = 1-diversity(OTU, index = "simpson")))
    }
    if ("InvSimpson" %in% measures) {
        outlist <- c(outlist, list(invsimpson = diversity(OTU, index = "simpson")))
    }
    if ("Fisher" %in% measures) {
        fisher = tryCatch(fisher.alpha(OTU, se = TRUE), warning = function(w) {
            warning("phyloseq::estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
            suppressWarnings(fisher.alpha(OTU, se = TRUE)[, c("alpha", "se")])})
        if (!is.null(dim(fisher))) {
            colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
            outlist <- c(outlist, list(fisher))
        }
        else {
            outlist <- c(outlist, Fisher = list(fisher))
        }
    }
	
    out = do.call("cbind", outlist)
    namechange = intersect(colnames(out), names(renamevec))
    colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
    colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), ignore.case = TRUE)
    out = out[, sort(unique(unlist(colkeep))), drop = FALSE]
    out <- as.data.frame(out)
    return(out)
	
}
