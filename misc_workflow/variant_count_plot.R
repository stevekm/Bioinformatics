#!/usr/bin/env Rscript
annot_files <- c("/ifs/data/molecpathlab/cancer1/sns_WES_copy-steve/VCF-GATK-HC-annot.all.txt",
                 "/ifs/data/molecpathlab/cancer1/Normal_WES/sns_out__bed/VCF-GATK-HC-annot.all.txt",
                 "/ifs/data/molecpathlab/braincancer-analysis/WES_analysis/results2/VCF-GATK-HC-annot.all.txt")
print(annot_files)


# load all the data

cancer1_df <- read.delim(file = "/ifs/data/molecpathlab/cancer1/sns_WES_copy-steve/VCF-GATK-HC-annot.all.txt", 
                          header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)

cancer1_df[["SET"]] <- "cancer1"

Normal_df <- read.delim(file = "/ifs/data/molecpathlab/cancer1/Normal_WES/sns_out__bed/VCF-GATK-HC-annot.all.txt", 
                          header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)

Normal_df[["SET"]] <- "Normal"

braincancer_df <- read.delim(file = "/ifs/data/molecpathlab/braincancer-analysis/WES_analysis/results2/VCF-GATK-HC-annot.all.txt", 
                       header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)

braincancer_df[["SET"]] <- "braincancer"

df <- rbind(cancer1_df, Normal_df)

df <- rbind(df, braincancer_df)

# export the combinded data
write.table(x = df, file = "VCF-GATK-HC-annot.all.merged.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


dput(colnames(df))
# c("#MUT", "SAMPLE", "CHR", "POS", "QUAL", "DEPTH", "FREQ", "Ref", 
#   "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", 
#   "ExonicFunc.refGene", "AAChange.refGene", "dbSNP_147", "gnomAD_exome_ALL", 
#   "gnomAD_genome_ALL", "Kaviar_AF", "cosmic80", "CADD13_PHRED", 
#   "FATHMM_noncoding", "FATHMM_coding", "SET")

dim(df)
# [1] 1615934      24


# add a column for dbSNP present/absence
df[["dbSNP"]] <- NA
df[which(df[["dbSNP_147"]] == '.'), "dbSNP"] <- "absent"
df[which(df[["dbSNP_147"]] != '.'), "dbSNP"] <- "present"

# subset for only chr1 entries
df_chr1 <- df[which(df[["CHR"]] == "chr1"), ]



dput(colnames(df_chr1))
# c("#MUT", "SAMPLE", "CHR", "POS", "QUAL", "DEPTH", "FREQ", "Ref", 
#   "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", 
#   "ExonicFunc.refGene", "AAChange.refGene", "dbSNP_147", "gnomAD_exome_ALL", 
#   "gnomAD_genome_ALL", "Kaviar_AF", "cosmic80", "CADD13_PHRED", 
#   "FATHMM_noncoding", "FATHMM_coding", "SET", "dbSNP")
dim(df_chr1)
# [1] 174515     24
# 




# make plot
library("ggplot2")
pdf("VCF-GATK-HC-dbSNP-chr1.pdf")
print(ggplot(data = df_chr1, aes(x = SAMPLE, fill = dbSNP)) + geom_bar() + coord_flip() + facet_grid(.~SET) + ggtitle("Chr1 Variants - WES GATK HaploType Caller"))
dev.off()



# ~~~~~~~ SUMMARY TABLES ~~~~~~ # 
# aggregate the analysis summary tables as well and plot those too...

cancer1_summary <- read.delim(file = "/ifs/data/molecpathlab/cancer1/sns_WES_copy-steve/summary-combined.wes.csv", 
                               header = TRUE, sep = ',', check.names = FALSE, stringsAsFactors = FALSE)


Normal_summary <- read.delim(file = "/ifs/data/molecpathlab/cancer1/Normal_WES/sns_out__bed/summary-combined.wes.csv",
                            header = TRUE, sep = ',', check.names = FALSE, stringsAsFactors = FALSE)


braincancer_summary <- read.delim(file = "/ifs/data/molecpathlab/braincancer-analysis/WES_analysis/results2/summary-combined.wes.csv",
                                 header = TRUE, sep = ',', check.names = FALSE, stringsAsFactors = FALSE)

Normal_summary[["SET"]] <- "Normal"
cancer1_summary[["SET"]] <- "cancer1"
braincancer_summary[["SET"]] <- "braincancer"


summary_df <- rbind(Normal_summary, cancer1_summary)
summary_df <- rbind(summary_df, braincancer_summary)

dput(colnames(summary_df))
# c("#SAMPLE", "R1 RAW READS", "R2 RAW READS", "TRIM INPUT READS", 
#   "TRIM SURVIVING READS", "TRIM SURVIVING READS %", "INPUT READS", 
#   "MAPPED READS (MQ10)", "MAPPED %", "CHIMERIC %", "MAPPED READS", 
#   "DEDUPLICATED READS", "DUPLICATES %", "ON-TARGET", "ON-TARGET 100BP PAD", 
#   "ON-TARGET 500BP PAD", "MEAN COVERAGE", "MEDIAN COVERAGE", "%_bases_above_10", 
#   "%_bases_above_50", "%_bases_above_100", "%_bases_above_500", 
#   "SET")
dim(summary_df)
# [1] 50 23


# combine with dbSNP values
# [1] "SAMPLE" "SET"    "dbSNP" 
dbSNP_present <- aggregate(dbSNP ~ SAMPLE + SET, data = df_chr1[which(df_chr1[["dbSNP"]] == "present"), ], FUN = length)
colnames(dbSNP_present) <- c("SAMPLE", "SET", "chr1_dbSNP_present")
dbSNP_absent <- aggregate(dbSNP ~ SAMPLE + SET, data = df_chr1[which(df_chr1[["dbSNP"]] == "absent"), ], FUN = length)
colnames(dbSNP_absent) <- c("SAMPLE", "SET", "chr1_dbSNP_absent")

summary_df <- merge(x = summary_df, y = dbSNP_present, by.x = c("#SAMPLE", "SET"), by.y = c("SAMPLE", "SET"), all.x = TRUE)
summary_df <- merge(x = summary_df, y = dbSNP_absent, by.x = c("#SAMPLE", "SET"), by.y = c("SAMPLE", "SET"), all.x = TRUE)

library('data.table')
setnames(x = summary_df, old = c('#SAMPLE', "DEDUPLICATED READS", "MEAN COVERAGE"), 
         new = c("SAMPLE", "DEDUPLICATED_READS", "MEAN_COVERAGE"))


dput(colnames(summary_df))
# c("SAMPLE", "SET", "R1 RAW READS", "R2 RAW READS", "TRIM INPUT READS", 
#   "TRIM SURVIVING READS", "TRIM SURVIVING READS %", "INPUT READS", 
#   "MAPPED READS (MQ10)", "MAPPED %", "CHIMERIC %", "MAPPED READS", 
#   "DEDUPLICATED_READS", "DUPLICATES %", "ON-TARGET", "ON-TARGET 100BP PAD", 
#   "ON-TARGET 500BP PAD", "MEAN_COVERAGE", "MEDIAN COVERAGE", "%_bases_above_10", 
#   "%_bases_above_50", "%_bases_above_100", "%_bases_above_500", 
#   "chr1_dbSNP_present", "chr1_dbSNP_absent")

# melt df for plot
library("reshape2")
summary_df_long <- reshape2::melt(data = summary_df[, c("SAMPLE", "SET", "DEDUPLICATED_READS", "MEAN_COVERAGE", "chr1_dbSNP_present")], 
               id.vars = c("SAMPLE", "SET"), variable.name = "type", value.name = 'value')

library("ggplot2")
summ_plot <- ggplot(data = summary_df_long, aes(x = SAMPLE, y = value, fill = type)) + geom_bar(stat = "identity")+ coord_flip() + facet_grid(SET~type, scales="free")

pdf("summary-combined.wes_merged.pdf", width = 10, height = 10)
print(summ_plot)
dev.off() 

write.table(x = summary_df, file = "summary-combined.wes_merged.csv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ',')
save.image("image.Rdata")
