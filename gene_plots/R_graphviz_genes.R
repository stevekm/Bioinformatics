#!/usr/bin/env Rscript

# visualize a character matrix
# output:
# ![genes_diagram](https://user-images.githubusercontent.com/10505524/34390683-508e08a8-eb0f-11e7-92f4-7b0ee9a4f3d5.png)

# install.packages("DiagrammeR")
# http://rich-iannone.github.io/DiagrammeR/ndfs_edfs.html
# https://github.com/rich-iannone/DiagrammeR#using-data-frames-to-define-graphviz-graphs
library("DiagrammeR")


data <- structure(list(gene = structure(c(8L, 18L, 10L, 3L, 4L, 11L, 
17L, 12L, 9L, 13L, 14L, 15L, 7L, 16L, 6L, 1L, 2L, 5L), 
.Label = c("ABHD17B", "ANKMY1_Average", "AURKAIP1", 
"BCL7C", "F11R_Average", "GEMIN8", "LMTK2", 
"MED20_Average", "PHC3_Average", "PISD", 
"PTCRA_Average", "PTPN9", "SCCPDH_Average", 
"SOCS2_Average", "SP2_Average", "TIMM10B", 
"VPS53_Average", "YIPF2_Average"), 
class = "factor"), 
sample = structure(c(2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L,4L),
.Label = c("Anhydroicaritin", "Glycyrrhizic_acid_rep_1", 
"Hydroxysafflor_yellow_A", "Hyperoside"), 
class = "factor")), 
.Names = c("gene","sample"),
class = "data.frame", row.names = c(NA, -18L))

head(data)
# gene                  sample
# 1 MED20_Average Glycyrrhizic_acid_rep_1
# 2 YIPF2_Average Glycyrrhizic_acid_rep_1
# 3          PISD Glycyrrhizic_acid_rep_1
# 4      AURKAIP1 Glycyrrhizic_acid_rep_1
# 5         BCL7C Glycyrrhizic_acid_rep_1
# 6 PTCRA_Average Hydroxysafflor_yellow_A


uniquenodes <- unique(c(as.character(data[["gene"]]), as.character(data[["sample"]])))
nodes <- create_node_df(n = length(uniquenodes), 
                        type = "number", 
                        label = uniquenodes)
edges <- create_edge_df(from = match(as.character(data[["sample"]]), uniquenodes), 
                        to = match(as.character(data[["gene"]]), uniquenodes), 
                        rel = "related")
g <- create_graph(nodes_df=nodes, 
                  edges_df=edges)
render_graph(g)
# devtools::install_github('rich-iannone/DiagrammeRsvg')
# install.packages("rsvg")
export_graph(g, "genes_diagram.png")
