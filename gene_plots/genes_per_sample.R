# http://rich-iannone.github.io/DiagrammeR/ndfs_edfs.html
# install.packages("DiagrammeR")
library("DiagrammeR")

data <- read.delim(file = 'gene_list.tsv', header = TRUE, sep = '\t', stringsAsFactors = FALSE)

uniquenodes <- unique(c(data[["genes"]], data[["sample"]]))

colors <- c( rep("blue", length(which(uniquenodes %in% data[["genes"]]))), 
             rep("red", length(which(uniquenodes %in% data[["sample"]]))))

nodes <- create_node_df(n = length(uniquenodes), 
                        type = "number", 
                        label = uniquenodes, 
                        color = colors, fillcolor = "white")



edges <- create_edge_df(from = match(as.character(data[["sample"]]), uniquenodes), 
                        to = match(as.character(data[["genes"]]), uniquenodes), 
                        rel = "related")
g <- create_graph(nodes_df = nodes, 
                  edges_df = edges)
render_graph(graph = g, title = "Genes per Sample")

# devtools::install_github('rich-iannone/DiagrammeRsvg')
# install.packages("rsvg")
export_graph(graph = g, file_name = "samples_genes.png", title = "Genes per Sample")

