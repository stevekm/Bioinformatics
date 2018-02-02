#!/usr/bin/env Rscript
annot_files <- c("/ifs/data/molecpathlab/cancer1/sns_WES_copy-steve/VCF-GATK-HC-annot.all.txt",
                 "/ifs/data/molecpathlab/cancer1/Normal_WES/sns_out_bed/VCF-GATK-HC-annot.all.txt",
                 "/ifs/data/molecpathlab/brain-analysis/WES_analysis/results2/VCF-GATK-HC-annot.all.txt")
print(annot_files)


# load all the data

cancer1_df <- read.delim(file = "/ifs/data/molecpathlab/cancer1/sns_WES_copy-steve/VCF-GATK-HC-annot.all.txt", 
                          header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)

cancer1_df[["SET"]] <- "cancer1"

Normal_df <- read.delim(file = "/ifs/data/molecpathlab/cancer1/Normal_WES/sns_out_bed/VCF-GATK-HC-annot.all.txt", 
                          header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)

Normal_df[["SET"]] <- "Normal"

brain_df <- read.delim(file = "/ifs/data/molecpathlab/brain-analysis/WES_analysis/results2/VCF-GATK-HC-annot.all.txt", 
                       header = TRUE, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)

brain_df[["SET"]] <- "brain"

df <- rbind(cancer1_df, Normal_df)

df <- rbind(df, brain_df)

# export the combinded data
write.table(x = df, file = "VCF-GATK-HC-annot.all.merged.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
save.image()

dput(colnames(df))
# c("#MUT", "SAMPLE", "CHR", "POS", "QUAL", "DEPTH", "FREQ", "Ref", 
#   "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", 
#   "ExonicFunc.refGene", "AAChange.refGene", "dbSNP_147", "gnomAD_exome_ALL", 
#   "gnomAD_genome_ALL", "Kaviar_AF", "cosmic80", "CADD13_PHRED", 
#   "FATHMM_noncoding", "FATHMM_coding", "SET")



# add a column for dbSNP present/absence
df[["dbSNP"]] <- NA
df[which(df[["dbSNP_147"]] == '.'), "dbSNP"] <- "absent"
df[which(df[["dbSNP_147"]] != '.'), "dbSNP"] <- "present"

# subset for only chr1 entries
df_chr1 <- df[which(df[["CHR"]] == "chr1"), ]

# make plot
library("ggplot2")
pdf("VCF-GATK-HC-dbSNP-chr1.pdf")
print(ggplot(data = df_chr1, aes(x = SAMPLE, fill = dbSNP)) + geom_bar() + coord_flip() + facet_grid(.~SET) + ggtitle("Chr1 Variants - WES GATK HaploType Caller"))
dev.off()







