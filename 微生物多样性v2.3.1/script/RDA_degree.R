args <- commandArgs(trailingOnly=T)

if(length(args) != 3){
        
        cat(
        " ============== Calculate RDA sample and env factors degree =========
        Usage: 
        Rscript RDA_degree.R position_one position_two outputfile
                parameters ->
                     position_one:  [file -> example /home/lhj/project/metagenome/17B0122A-2/EnvironmentFactors/RDA_CCA/genus/sample.position.csv];
                     position_two:  [file -> example /home/lhj/project/metagenome/17B0122A-2/EnvironmentFactors/RDA_CCA/genus/environmentfactors.position.csv];
                     outputfile: [file -> output file example ]; \n")
        options("show.error.messages" = F) 
        stop()
        
}


sample_positon <- args[1]
env_position   <- args[2]
outputfile     <- args[3]

sample <- read.table(sample_positon, header = T, sep = ",", row.names = 1)

env    <- read.table(env_position, header = T, sep = ",", row.names = 1)


a <- lapply(1:nrow(sample), function(t){
	id <- rownames(sample)[t]
	tmp <- t(sapply(1:nrow(env), function(m){
		 env_id <- rownames(env)[m]
		 data <- matrix(c(0, 0, unlist(sample[t, ]), unlist(env[m, ])), nrow = 3, byrow = T)

		 len  <- as.vector(dist(data))
		 val  <- (len[1] ^ 2 + len[2] ^ 2 - len[3] ^ 2 )  / (2 * len[1] * len[2])
		 deg  <- acos(val) * 180 / pi
		 c(id, env_id, deg)
	}))
	data.frame(tmp, stringsAsFactors = F)

})

res <- do.call(rbind, a)
res[[3]] <- as.numeric(res[[3]])
res$type <- "positive"
res$type[res[[3]] > 90] <- "negative"

colnames(res) <- c("ID", "Env", "Degree", "Type")

write.table(res, outputfile, row.names = F, quote = F, sep = "\t")

