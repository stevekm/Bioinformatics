# test the deconstructSigs package 
# https://github.com/raerose01/deconstructSigs
# deconstructSigs aims to determine the contribution of known mutational processes to a tumor sample. By using deconstructSigs, one can:
# 
# Determine the weights of each mutational signature contributing to an individual tumor sample
# 
# Plot the reconstructed mutational profile (using the calculated weights) and compare to the original input sample

# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
# install.packages("deconstructSigs")
library("deconstructSigs")


# ~~~~~ LOAD DATA ~~~~~ # 
# broad public RNASeq tumor dataset from here https://portals.broadinstitute.org/ccle/data
ccle_data_file <- "ccle2maf_081117.txt"
ccle <- read.delim(file = ccle_data_file, header = TRUE, sep = '\t', check.names = FALSE)

dim(ccle)
# [1] 1159663      32

dput(colnames(ccle))
# c("Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chromosome", 
#   "Start_position", "End_position", "Strand", "Variant_Classification", 
#   "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "dbSNP_RS", 
#   "dbSNP_Val_Status", "Genome_Change", "Annotation_Transcript", 
#   "Tumor_Sample_Barcode", "cDNA_Change", "Codon_Change", "Protein_Change", 
#   "isDeleterious", "isTCGAhotspot", "TCGAhsCnt", "isCOSMIChotspot", 
#   "COSMIChsCnt", "ExAC_AF", "WES_AC", "WGS_AC", "SangerWES_AC", 
#   "SangerRecalibWES_AC", "RNAseq_AC", "HC_AC", "RD_AC")

length(levels(ccle[["Tumor_Sample_Barcode"]]))
# [1] 1463
head(levels(ccle[["Tumor_Sample_Barcode"]]))
# [1] "201T_LUNG"                     "22RV1_PROSTATE"                "2313287_STOMACH"              
# [4] "253J_URINARY_TRACT"            "253JBV_URINARY_TRACT"          "42MGBA_CENTRAL_NERVOUS_SYSTEM"

levels(ccle[["Chromosome"]])
# [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "21" "22" "3"  "4"  "5"  "6"  "7"  "8" 
# [22] "9"  "M"  "X"  "Y" 

# check number of entries per sample
length(table(ccle[["Tumor_Sample_Barcode"]]))
# [1] 1463

# check the number of tumor types with less than 50 entries
length(which(data.frame(table(ccle[["Tumor_Sample_Barcode"]]))[["Freq"]] < 50))
# [1] 13




# ~~~~~ CONVERT DATA ~~~~~ # 
# fix the chrom names
ccle[["Chromosome"]] <- sprintf('chr%s', as.character(ccle[["Chromosome"]]))

# keep only entries with chroms in the reference data
ccle <- ccle[which(ccle[["Chromosome"]] %in% seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)), ]

# need at least 50 variants per sample, remove samples that dont have enough variants
# actually looks like you need at least 55 variants? 
tumor_samples <- data.frame(table(ccle[["Tumor_Sample_Barcode"]]))
remove_samples <- as.character(tumor_samples[which(tumor_samples[["Freq"]] < 55), "Var1"])
ccle <- ccle[which(! as.character(ccle[["Tumor_Sample_Barcode"]]) %in% remove_samples), ]
ccle <- droplevels(ccle)


# convert to signatures format
sigs.input <- mut.to.sigs.input(mut.ref = ccle, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele1")

# I dunno what this means:
# Warning message:
#     In mut.to.sigs.input(mut.ref = ccle, sample.id = "Tumor_Sample_Barcode",  :
#                              Check ref bases -- not all match context:
#                              22RV1_PROSTATE:chrM:9477:G:A, 22RV1_PROSTATE:chrM:14323:G:A, 253JBV_URINARY_TRACT:chrM:10589:G:A, 253JBV_URINARY_TRACT:chrM:13768:T:C, AML193_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:12720:A:G, AMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:6179:G:A, AMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:8684:C:T, AMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:14470:T:C, AMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:15148:G:A, AMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:15355:G:A, AMO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:15487:A:T, BICR31_UPPER_AERODIGESTIVE_TRACT:chrM:14793:A:G, BICR31_UPPER_AERODIGESTIVE_TRACT:chrM:15218:A:G, BICR56_UPPER_AERODIGESTIVE_TRACT:chrM:8697:G:A, BICR56_UPPER_AERODIGESTIVE_TRACT:chrM:13928:G:C, BICR56_UPPER_AERODIGESTIVE_TRACT:chrM:14905:G:A, BICR6_UPPER_AERODIGESTIVE_TRACT:chrM:6365:T:C, BL41_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:11251:A:G, BL41_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE:chrM:12612:A:G, BL41_HAEMATOPOIETIC_AND_LY [... truncated]




# ~~~~~ FIND SIGNATURES ~~~~~ # 
# Determine the signatures contributing an already normalized sample, just the first 5 tumor types
#  using signatures.nature2013
sample_signatures <- list()
pdf(file = "mut_signatures.pdf")
for(i in seq(5)){
    sampleID <- rownames(sigs.input)[i]
    signatures <- whichSignatures(tumor.ref = sigs.input, 
                                  signatures.ref = signatures.nature2013, 
                                  sample.id = sampleID,
                                  contexts.needed = TRUE,
                                  tri.counts.method = 'default')
    
    # make plots
    print(plotSignatures(signatures, sub = 'signatures.nature2013'))
    
    print(makePie(signatures, sub = 'signatures.nature2013'))
    
    # add it to the list
    sample_signatures[[sampleID]] <- signatures
}
dev.off()

