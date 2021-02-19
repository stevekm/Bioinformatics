#!/usr/bin/env Rscript
# Script for plotting the Nextflow trace.txt file to get some aggregate metrics such as execution time per sample
# depends on the 'tag' for each process having ':' delim values
# USAGE:
# $ module load R/R-3.6.3
# $ ./plot_trace.R results_aggregate/trace.txt
library("ggplot2")
args <- commandArgs(T)
input_file <- args[1]
trace <- read.delim(file = input_file, header = TRUE, sep = '\t', check.names = FALSE)
save.image(".loaded.Rdata")

# get all the tags and split them on ':'
tags <- as.data.frame(do.call(rbind, strsplit(x = as.character(trace[["tag"]]), split = ':')))
colnames(tags) <- c("num_threads", "num_samples", "total_samples")
trace <- cbind(trace, tags)

# fix some data for plotting & aggregating
trace[["num_threads"]] <- factor(x = trace[["num_threads"]], levels = sort(as.numeric(unique(levels(trace[["num_threads"]])))))
trace[["num_samples"]] <- factor(x = trace[["num_samples"]], levels = sort(as.numeric(unique(levels(trace[["num_samples"]])))))

trace[["realtime"]] <- as.numeric(as.character(trace[["realtime"]]))

trace[["time_per_sample"]] <- trace[["realtime"]] / as.numeric(as.character(trace[["num_samples"]]))

trace[["n"]] <- 1

# make plots
pdf(file = "total_time.pdf", width = 8, height = 8)
ggplot(data = trace, aes(x = num_threads,
                         y = ((realtime / 1000) / 60),
                         fill = num_samples)) +
    geom_bar(stat = 'identity', position = "dodge") +
    ggtitle("GetBaseCountsMultiSample Execution Time") +
    ylab("time (min)") +
    theme_bw()
dev.off()


pdf(file = "time_per_sample.pdf", width = 8, height = 8)
ggplot(data = trace, aes(x = num_threads,
                         y = ((time_per_sample / 1000) / 60) ,
                         fill = num_samples)) +
    geom_bar(stat = 'identity', position = "dodge") +
    ggtitle("GetBaseCountsMultiSample Execution Time") +
    ylab("time per sample (min)") +
    theme_bw()

dev.off()


# make aggregate metrics table
agg <- aggregate(time_per_sample ~ num_samples * num_threads, trace, mean)
sd_df <- aggregate(time_per_sample ~ num_samples * num_threads, trace, sd)
colnames(sd_df) <- c("num_samples", "num_threads", "sd")
agg <- merge(agg, sd_df, all = TRUE)

n_df <- aggregate(n ~ num_samples * num_threads, trace, sum)
agg <- merge(agg, n_df, all = TRUE)

agg <- agg[ order(agg[["time_per_sample"]]) , ]

write.table(x = agg, file = "aggregate.tsv", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

save.image(".final.Rdata")
