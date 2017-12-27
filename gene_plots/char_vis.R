#!/usr/bin/env Rscript

# visualize a character matrix
# example output: 
# ![genes](https://user-images.githubusercontent.com/10505524/34389094-9ef7b890-eb05-11e7-962b-321e03b2c5f4.png)

library("ggplot2")
# install.packages("ggrepel")
library("ggrepel") # for spreading text labels on the plot

lines <- "
MED20_Average   Glycyrrhizic_acid_rep_1

YIPF2_Average   Glycyrrhizic_acid_rep_1

PISD    Glycyrrhizic_acid_rep_1

AURKAIP1    Glycyrrhizic_acid_rep_1

BCL7C   Glycyrrhizic_acid_rep_1

PTCRA_Average   Hydroxysafflor_yellow_A

VPS53_Average   Hydroxysafflor_yellow_A

PTPN9   Hydroxysafflor_yellow_A

PHC3_Average    Anhydroicaritin

SCCPDH_Average  Anhydroicaritin

SOCS2_Average   Anhydroicaritin

SP2_Average Anhydroicaritin

LMTK2   Anhydroicaritin

TIMM10B Anhydroicaritin

GEMIN8  Anhydroicaritin

ABHD17B Anhydroicaritin

ANKMY1_Average  Hyperoside

F11R_Average    Hyperoside

"

con <- textConnection(lines)
data <- read.delim(con, header = FALSE, sep = "")
close(con)

colnames(data) <- c("gene", "sample")

head(data)
# gene                  sample
# 1 MED20_Average Glycyrrhizic_acid_rep_1
# 2 YIPF2_Average Glycyrrhizic_acid_rep_1
# 3          PISD Glycyrrhizic_acid_rep_1
# 4      AURKAIP1 Glycyrrhizic_acid_rep_1
# 5         BCL7C Glycyrrhizic_acid_rep_1
# 6 PTCRA_Average Hydroxysafflor_yellow_A

str(data)
# 'data.frame':	18 obs. of  2 variables:
#     $ gene  : Factor w/ 18 levels "ABHD17B","ANKMY1_Average",..: 8 18 10 3 4 11 17 12 9 13 ...
# $ sample: Factor w/ 4 levels "Anhydroicaritin",..: 2 2 2 2 2 3 3 3 1 1 ...

ggplot(data = data, aes(y = as.numeric(gene), x = as.numeric(gene), label = gene)) + 
    geom_dotplot(alpha = 0) + 
    facet_grid(sample~.) + 
    coord_cartesian(ylim = c(max(as.numeric(data[["gene"]])) + 1, 0 )) +
    theme(
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          ) +
    geom_text_repel(aes(y = 0, label = gene), show.legend = FALSE, segment.alpha = 0, force = 5) + 
    ggtitle("Genes per Sample")
