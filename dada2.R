library(BiocManager)
# To install dada2 bioconductor R package
# BiocManager::install('dada2')
library(dada2)

# path for the folder that contains your Fastq files, it is better to be 
# in the same folder with your dada2.R file
path<- 'data' # 


list.files(path)

# Forward and reverse fastq filenames have format: 
# SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[c(1,4)])
plotQualityProfile(fnRs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=F)
head(out)
errF <- learnErrors(filtFs, multithread=TRUE) 

errR <- learnErrors(filtRs, multithread=TRUE) 
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers) # dim(seqtab)  91 6617
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

sum(seqtab.nochim)/sum(seqtab)  # 0.9852159
###
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) # summarizing the results
# 


taxa <- assignTaxonomy(seqtab.nochim, "RefSeq-RDP16S_v2_May2018.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)




rownames(seqtab.nochim)  # [1] "0004L0" "0009L0"

unqs.mock <- seqtab.nochim["0004L0",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock

cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
# DADA2 inferred 323 sample sequences present in the Mock community.

##
unqs.mock <- seqtab.nochim["0009L0",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
##
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
## DADA2 inferred 109 sample sequences present in the Mock community.
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
#########
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
# phyloseq object

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
###remove mock samples
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
##
plot_richness(dietswap, x="group", measures=c("Shannon", "Simpson"), color="sex")



