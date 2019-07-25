rm(list = ls())
args = commandArgs(trailingOnly=T)
if(length(args) != 4){

    cat(
    "============== do OTU cor ===========
    Usage: 
    Rscript stringdb.r diff_file organ tmpfile_path outfile_path
        parameters ->
                   diff_file   : [file -> always diff.xls];
                   organ       : [name -> for example Homo_sapiens/Mus_musculus/Rattus_norvegicus];
                   tmpfile_path：[path -> path for tmp];
                   outfile_path: [path -> path for output]; \n")
    options("show.error.messages" = F) 
    stop()

}

library(STRINGdb)

inputfile  = normalizePath(args[1])
organ      = args[2]
tmpputpath = normalizePath(args[3])
outputpath = normalizePath(args[4])

input = read.table(inputfile, header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "", quote = "", fill = T)
data  = input[input$type != "Not DEG" , ]
data$gene = rownames(data)

#获取数据库中物种的species_id
name2id    = get_STRING_species(version = "10", species_name = NULL)
organ      = gsub("_", " ", organ)
species_id = name2id$species_id[name2id$official_name == organ]

#加载数据库，将gene名和数据库中的蛋白标识符mapping
string_db   = STRINGdb$new(version = "10", species = species_id, score_threshold = 0, input_directory = "")
data_mapped = string_db$map(data,"gene",removeUnmappedRows = TRUE)
hits        = data_mapped$STRING_id



#ppi分析

cat("##### ppi start #####\n")

if(length(hits) < 400) {
    pdf(paste(outputpath, "ppi_string.pdf",sep = "/"))
    string_db$plot_network(hits) # ppi分析，绘图，该步骤会从官网下载数据库，大约35MB
    dev.off()
}



inter    = string_db$get_interactions(hits)
links    = data.frame(from  = unlist(lapply(inter$from, function(t){ data_mapped[data_mapped$STRING_id == t, 1][1] })),
                      to    = unlist(lapply(inter$to, function(t){ data_mapped[data_mapped$STRING_id == t, 1][1] })),
                      score = inter[,ncol(inter)], stringsAsFactors=FALSE )

#去除重复的边
res1 = unlist(lapply(1:nrow(links), function(t) {

  link_sort = sort(links[t,1:2]) 
  res = paste(link_sort[1],link_sort[2],sep = "_")

}))

links$id = res1
link_out = links[!duplicated(links[,4]),]
LINK = link_out[,1:3]
colnames(LINK) = c("gene1", "gene2", "score")

write.table(LINK, paste(tmpputpath, "interactions.xls", sep = "/"), sep = "\t", quote = F, row.names = F)


#igraph构建网络，聚类
library(igraph)
net = graph.data.frame(LINK, directed = F)
cfg = cluster_fast_greedy(net)

#输出结果
tmp = as.data.frame(as.matrix(membership(cfg)))
tmp$gene = rownames(tmp)
colnames(tmp) = c("group", "gene")
res = tmp[order(tmp[,1]),]
cluster_out = paste(outputpath, "cluster.xls", sep = "/")
write.table(res,cluster_out, sep = "\t", quote = F, row.names = F)

#绘图
pdf(paste(outputpath,"ppi.pdf",sep = "/"), height = 6, width = 8)
plot(cfg, net, vertex.size = 10, vertex.label.cex = 0.3 )
dev.off()




