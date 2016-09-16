


# print(Diff_plot_workflow(sample_file_list))

# pdf(file = paste0(diff_cl4_ctrl_file,"_volcano_plot.pdf"),height = 8,width = 8)
# mycat(paste0("## ", basename(diff_cl4_ctrl_file), ' {.tabset}\n'))
# diff_df <- read.delim(file = diff_cl4_ctrl_file,header = TRUE,sep = ',')
# DiffBind_volcano_plotly_top_protein_coding_promoters(diff_df = diff_df)
# DiffBind_volcano_plot_top_protein_coding_promoters(diff_df = diff_df)
# Diff_stats(diff_df = diff_df)
# # my_plot <- DiffBind_volcano_plotly_top_protein_coding_promoters(diff_df = diff_df)
# # htmltools::tagList(list(as.widget(my_plot)))
# # dev.off()
# 
# # pdf(file = paste0(diff_cl4_cl5_file,"_volcano_plot.pdf"),height = 8,width = 8)
# mycat(paste0("## ", basename(diff_cl4_cl5_file), ' {.tabset}\n'))
# diff_df <- read.delim(file = diff_cl4_cl5_file,header = TRUE,sep = ',')
# # DiffBind_volcano_plot_top_protein_coding_promoters(diff_df = diff_df)
# DiffBind_volcano_plotly_top_protein_coding_promoters(diff_df = diff_df)
# Diff_stats(diff_df = diff_df)
# # my_plot <- DiffBind_volcano_plotly_top_protein_coding_promoters(diff_df = diff_df)
# # htmltools::tagList(list(as.widget(my_plot)))
# 
# # dev.off()
# 
# # pdf(file = paste0(diff_cl5_ctrl_file,"_volcano_plot.pdf"),height = 8,width = 8)
# mycat(paste0("## ", basename(diff_cl5_ctrl_file), ' {.tabset}\n'))
# diff_df <- read.delim(file = diff_cl5_ctrl_file,header = TRUE,sep = ',')
# # DiffBind_volcano_plot_top_protein_coding_promoters(diff_df = diff_df)
# DiffBind_volcano_plotly_top_protein_coding_promoters(diff_df = diff_df)
# Diff_stats(diff_df = diff_df)
# # my_plot <- DiffBind_volcano_plotly_top_protein_coding_promoters(diff_df = diff_df)
# # htmltools::tagList(list(as.widget(my_plot)))
# 
# # dev.off()

```{r, eval=FALSE}
mycat('# Plotly Test\n\n')
# https://plot.ly/r/knitr/
# install.packages("plotly")
suppressPackageStartupMessages(library("plotly"))
# p <- plot_ly(economics, x = date, y = unemploy / pop)
# p
p <- plot_ly(economics, x = date, y = uempmed, name = "unemployment")
p <- p %>% add_trace(y = fitted(loess(uempmed ~ as.numeric(date))), x = date)

p2 <- plot_ly(data = iris, x = Sepal.Length, y = Petal.Length, mode = "markers", color = Species)

htmltools::tagList(list(as.widget(p), as.widget(p2)))

# test the plotly..

# base R 
# # add base volcano plot; all points
#    with(diff_df_min, plot(Fold, -log10(p.value), pch=20, xlim=c(min(Fold)-1,max(Fold)+1),col=plot_colors[1],xlab = "log2(Fold Change)"))
#    
#    # Add colored points for data subsets
#    with(subset(diff_df_min, p.value<signif_p_value ), points(Fold, -log10(p.value), pch=20, col=plot_colors[2]))
#    with(subset(diff_df_min, abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[3]))
#    with(subset(diff_df_min, p.value<signif_p_value & abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[4]))

# # melt it into long format
#     if (melt_df) diff_log_ratio <- reshape2::melt(diff_log_ratio,id.vars="gene",variable.name="sample",value.name="diff_peak_log_ratio")
# 
# diff_df_min_sub <-diff_df_min[c("external_gene_name", "gene_biotype", "Fold", "p.value")]
# diff_df_min_sub["group"] <- "NotSignif"
# diff_df_min_sub[which(diff_df_min_sub['p.value'] < signif_p_value & abs(diff_df_min_sub['Fold']) < fold_change_value ),"group"] <- "signif"
# diff_df_min_sub[which(diff_df_min_sub['p.value'] > signif_p_value & abs(diff_df_min_sub['Fold']) > fold_change_value ),"group"] <- "FC"
# diff_df_min_sub[which(diff_df_min_sub['p.value'] < signif_p_value & abs(diff_df_min_sub['Fold']) > fold_change_value ),"group"] <- "signif_FC"
# 
# head(diff_df_min_sub, 100)
# 
# plot_ly(data = diff_df_min_sub, x = Fold, y = -log10(p.value), mode = "markers", color = group)
# <!-- # References -->
# citation_package: natbib
# bibliography: references.bib
# biblio-style: apsr
# <!-- [@DiffBind] -->
```


# to test:
# diff_df_min[which(diff_df_min['external_gene_name'] == "CROCC"),]
#        external_gene_name   gene_biotype Fold  p.value     group
# 12              CROCC protein_coding 5.46 3.44e-14 signif_fc


# which(diff_df_min)
# diff_df_min[which(diff_df_min['gene_biotype'] == "protein_coding"),]
# 
# head(within(subset(diff_df_min, gene_biotype == "protein_coding"), order(-Fold, p.value))) 
# 
# diff_df_min[with(diff_df_min[which(diff_df_min[['gene_biotype']] == "protein_coding", arr.ind = TRUE),], order(-Fold, p.value) ),"gene_biotype"][1:top_gene_number]
# 
# diff_df_min[with(diff_df_min[which(diff_df_min['gene_biotype'] == "protein_coding"),], order(-Fold, p.value) ),"gene_biotype"][1:top_gene_number] <- "top_signif_fc"
# diff_df_min[with(diff_df_min[which(diff_df_min['gene_biotype'] == "protein_coding"),], order(-Fold, p.value) ),"group"][1:top_gene_number] <- "top_signif_fc"
# # diff_df_min[with(diff_df_min, order(Fold, p.value)),"group"][1:top_gene_number] <- "top_signif_fc" 



# require(data.table)
# d <- data.table(mtcars, key="cyl")
# d[, head(.SD, 4), by=cyl]
# 
# setDT(diff_df_min,key = "external_gene_name")
# diff_df_min <- diff_df_min[gene_biotype == "protein_coding"][order(-Fold, p.value), head(.SD, 10)][,group := "top_signif_fc"]
# diff_df_min[gene_biotype == "protein_coding"][order(Fold, p.value), head(.SD, top_gene_number)][,group := "top_signif_fc"]

# diff_df_min[gene_biotype == "protein_coding" & order(Fold, p.value), head(.SD, top_gene_number)][,group := "top_signif_fc"]

# diff_df_min[gene_biotype == "protein_coding" & order(Fold, p.value)]

# mydf %>%
# group_by(yearmonth) %>%
# arrange(desc(count)) %>%
# slice(1:2)

# head(diff_df_min)
# dim(diff_df_min)

# diff_df_min %>%
#     filter(gene_biotype == "protein_coding") %>% # subset for protein coding genes
#     arrange(-Fold, p.value) %>% # Sort by Fold change, then by p value
#     slice(1:10) %>% # take the top 10 entries... 
#     mutate(group = "top_signif_fc") # ... and change the "group" column value to "top_signif_fc"
# 

# diff_df_minNew <- diff_df_min %>%
#                   mutate(ind = gene_biotype  != "protein_coding") %>% 
#                   arrange(ind, -Fold, p.value) %>% 
#                   mutate(group = ifelse(row_number() < 11, "top_signif_fc", group)) %>%
#                   select(-ind)

# p.value<signif_p_value & abs(Fold)>fold_change_value

```{r, eval=FALSE}
# mycat('# Plotly Test\n\n')
# https://plot.ly/r/knitr/
# install.packages("plotly")
suppressPackageStartupMessages(library("plotly"))
# p <- plot_ly(economics, x = date, y = unemploy / pop)
# p
p <- plot_ly(economics, x = date, y = uempmed, name = "unemployment")
p <- p %>% add_trace(y = fitted(loess(uempmed ~ as.numeric(date))), x = date)

p2 <- plot_ly(data = iris, x = Sepal.Length, y = Petal.Length, mode = "markers", color = Species)

htmltools::tagList(list(as.widget(p), as.widget(p2)))

```