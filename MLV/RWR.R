library("diffusr")


args <- commandArgs(TRUE)

graph_file <- args[1]
adjacency_matrix_file <- args[2]
result_file <- args[3]

table_file <- read.table(graph_file, header=TRUE, sep="\t")
temp_table <- table_file[,2:ncol(table_file)]
rownames(temp_table) <- table_file[,1]
graph <- as.matrix(temp_table)

adjacency_file <- read.table(adjacency_matrix_file, header=FALSE, sep="\t")
p0 <- as.vector(adjacency_file[,2])

pt <- random.walk(p0, graph)
pt <- pt$p.inf
rownames(pt) <- table_file[,1]

sink(result_file)
print (pt)
sink()

