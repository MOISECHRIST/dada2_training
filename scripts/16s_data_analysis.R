#!/usr/bin/R

## Import libraries
library(dada2)

## Load data
data_dir_path <- "data/MiSeq_SOP"
list.files(data_dir_path)
fnFs <- sort(list.files(data_dir_path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(data_dir_path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Quality check 
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

## Trimming and filtering 
result_path <- "results"
dir.create(file.path(result_path, "filtered"))
filtFs <- file.path(result_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(result_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(240, 160), maxN = 0, maxEE = c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
plot(reads.out ~ reads.in, out)
abline(a=0, b=1)

## Quality check 
plotQualityProfile(filtFs)
plotQualityProfile(filtRs)


## Lear Error rate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

## Plot errors 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


## Sample Inference
### NOTE : to benchmark various error correction models. 
### We will change the default errorEstimationFunction parameter by desired one.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, 
               errorEstimationFunction = loessErrfun)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, 
               errorEstimationFunction = loessErrfun)

### Inspect the sample result
dadaFs[[1]]
dadaRs[[1]]


## Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
View(seqtab)

### Distribution of sequence lengths
plot(table(nchar(getSequences(seqtab))), xlab = "Sequence length", ylab="Frequency")

## Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## Assign Taxonomy 
taxa <- assignTaxonomy(seqtab.nochim, "data/ref_db/SILVA-v138.2-16s/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
