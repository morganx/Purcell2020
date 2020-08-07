# Libraries ####
library(phyloseq)
library(dada2)

# Set working directory
setwd("/SCRATCH2/DAN_DATA/purcellProject/")
# Set path to data files
path <- "/SCRATCH2/DAN_DATA/purcellProject/data/MGS00332_Rachel Purcell_Delivery/MGS00332_1/processed"


# dada2 Code ####
# Assign files that end with –R1.fq.gz to forward (fnFs), –R2.fq.gz to reverse (FnRs)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz"))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz"))

# What are the actual names of our samples? Let’s parse them out of the files.
sample.names <- sapply(strsplit(sub("_", ".", fnFs), ".", fixed = T), `[`, 2)
sample.names <- sapply(strsplit(sample.names, "_"), `[`, 1)

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Use this command to filter and trim. It will truncate reads to an #indicated length
# (informed by your RQC data), truncate poor-quality #bases,
# and set the maximum acceptable error rate 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncQ = 2, maxN = 0, maxEE = 1,
                     rm.phix = FALSE, multithread = TRUE)
out
print("Filter and Trim Complete")

# Learn the profile of sequencing errors in the data
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)


derepFs <- derepFastq(filtFs, verbose = FALSE)
derepRs <- derepFastq(filtRs, verbose = FALSE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = FALSE)


seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN),
               rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track
dd
silva <- ("/DATABASES/SILVA/silva_v138/silva_nr_v138_train_set.fa")
taxa <- assignTaxonomy(seqtab.nochim, silva, tryRC = TRUE, multithread = TRUE)
unname(head(taxa))
print("dada2 Complete")

# Handoff to phyloseq ####
# Make phyloseq files
samples.out <- rownames(seqtab.nochim)
subject <- samples.out
samdf <- data.frame(Subject = subject)
rownames(samdf) <- samples.out

# Construct phyloseq object (no tree)
ps_notree <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               sample_data(samdf),
               tax_table(taxa))

# Save phyloseq object
saveRDS(ps_notree, file = "ps_notree.rds")
print("Phyloseq Object Created")

# Session Info ####
sessionInfo()
