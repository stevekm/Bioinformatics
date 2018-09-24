#!/usr/bin/env Rscript
library("ggplot2")
library("data.table")

args = commandArgs(TRUE)
aa_change_file <- args[1]
agg_counts_file <- args[2]
aa_pdf <- args[3]
agg_pdf <- args[4]

# load files
aa_df <- read.delim(aa_change_file, sep = '\t', header = FALSE)
agg_df <- read.delim(agg_counts_file, sep = '\t', header = FALSE)

# set column names
names(aa_df) <- c("AAChange", "Count", "Study", "Project", "Caller")
names(agg_df) <- c("Sample", "BAM_UUID", "Count", "Study", "Project", "Caller", "TMB")

# sort the Project levels based on median value
agg_medians <- aggregate( TMB ~ Project, data = agg_df, FUN = median)
agg_median_order <- as.character(agg_medians[order(agg_medians[["TMB"]]), ][["Project"]])
agg_df[["Project"]] <- factor(agg_df[["Project"]], levels = agg_median_order)

# get the number of cases per project
agg_df2 <- agg_df
agg_df2[["n"]] <- 1
agg_n <- aggregate( n ~ Project, data = agg_df2, FUN = sum)
agg_n[["Project"]] <- factor(agg_n[["Project"]], levels = agg_median_order)

# re order the AA change df
aa_df <- aa_df[order(aa_df[["Count"]], decreasing = TRUE), ]
aa_df[["Project"]] <- factor(aa_df[["Project"]], levels = agg_median_order)

# subset for the top 10 types for each Project
# aa_dt <- data.table(aa_df)
# aa_dt <- aa_dt[,.SD[order(-Count)[1:10]], by = .(Project, Caller)]
# setDF(aa_dt)
# aa_dt[["Project"]] <- factor(aa_dt[["Project"]], levels = agg_median_order)



# get the variant callers in project
callers <- unique(c(levels(agg_df[["Caller"]]), levels(aa_df[["Caller"]])))

for(caller in callers){
    # base plot
    agg_plot <- ggplot(data = agg_df[ which(agg_df[["Caller"]] == caller ), ], aes(x = Project, y = Count)) + 
        geom_violin() + 
        geom_boxplot(width=.1, outlier.size=0) +
        scale_y_log10() + 
        theme_bw() +
        theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        ggtitle(sprintf("Tumor Mutation Burden per sample across TCGA projects: %s\n\n", caller))
    
    aa_plot <- ggplot(data = aa_df[which(aa_df[["Caller"]] == caller), ], aes(x = Project, y = Count, fill = AAChange)) + 
        geom_bar(stat = 'identity', position = "fill") + 
        scale_y_continuous(labels = scales::percent) +
        theme(legend.position="bottom") + 
        theme_bw() + 
        theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        ylab("Proportion") + 
        ggtitle(sprintf("Amino Acid Change Proportions across TCGA projects: %s", caller))
    
    pdf(file = sprintf('%s.%s.pdf', agg_pdf, caller), width = 12)
    print(agg_plot)
    dev.off()
    
    pdf(file = sprintf('%s.%s.pdf', aa_pdf, caller), width = 22)
    print(aa_plot)
    dev.off()
}

save.image()

# plot, just for one variant caller
# caller <- "varscan"

# base plot
# agg_plot <- ggplot(data = agg_df[ which(agg_df[["Caller"]] == caller ), ], aes(x = Project, y = Count)) + 
#     geom_violin() + 
#     geom_boxplot(width=.1, outlier.size=0) +
#     scale_y_log10() + 
#     theme_bw() +
#     theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
#     ggtitle(sprintf("Tumor Mutation Burden per sample across TCGA projects: %s\n\n", caller))
# 
# # plot with n top axis
# agg_plot_n <- ggplot(data = agg_df[ which(agg_df[["Caller"]] == caller ), ], aes(x = Project, y = Count)) + 
#     geom_violin() + 
#     geom_boxplot(width=.1, outlier.size=0) +
#     scale_y_log10() + 
#     theme_bw() +
#     theme(panel.grid.minor = element_blank()) + 
#     ggtitle(sprintf("Tumor Mutation Burden per sample across TCGA projects: %s\n\n", caller)) + 
#     scale_x_discrete(name = "number of samples", position = "top", labels = agg_n[["n"]])
# 
# # https://stackoverflow.com/questions/36487283/ggplot2-2-1-0-broke-my-code-secondary-transformed-axis-now-appears-incorrectly/36761846#36761846
# # extract gtable
# # library("gtable")
# # library("grid")
# # https://stackoverflow.com/questions/21026598/ggplot2-adding-secondary-transformed-x-axis-on-top-of-plot
# # agg_plot_g <- ggplot_gtable(ggplot_build(agg_plot))
# # agg_plot_n_g <- ggplot_gtable(ggplot_build(agg_plot_n))
# # # overlap the panel of the 2nd plot on that of the 1st plot
# # pp <- c(subset(agg_plot_g$layout, name=="panel", se=t:r))
# # g <- gtable_add_grob(agg_plot_g, agg_plot_n_g$grobs[[which(agg_plot_n_g$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
# # 
# # g <- gtable_add_grob(agg_plot_g, agg_plot_g$grobs[[which(agg_plot_g$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
# # grid.draw(g)
# # 
# # g <- gtable_add_grob(g1, g1$grobs[[which(g1$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
# ## extract gtable
# # g1 <- ggplot_gtable(ggplot_build(p1))
# # g2 <- ggplot_gtable(ggplot_build(p2))
# # 
# # ## overlap the panel of the 2nd plot on that of the 1st plot
# # pp <- c(subset(g1$layout, name=="panel", se=t:r))
# # 
# # g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b,
# #                      pp$l)
# 
# 
# aa_plot <- ggplot(data = aa_df[which(aa_df[["Caller"]] == caller), ], aes(x = Project, y = Count, fill = AAChange)) + 
#     geom_bar(stat = 'identity', position = "fill") + 
#     scale_y_continuous(labels = scales::percent) +
#     theme(legend.position="bottom") + 
#     theme_bw() + 
#     theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#     ylab("Proportion") + 
#     ggtitle(sprintf("Amino Acid Change Proportions across TCGA projects: %s", caller))
# 
# pdf(file = agg_pdf, width = 12)
# print(agg_plot)
# dev.off()
# 
# pdf(file = aa_pdf, width = 22)
# print(aa_plot)
# dev.off()


