## Best to create conda environment
## install.packages("BiocManager")
## BiocManager::install("dada2")
## BiocManager::install("decontam")
## BiocManager::install("phyloseq")
## BiocManager::install("tidyverse")
## BiocManager::install("microbiome")

library(dada2); packageVersion("dada2")
library(tidyverse)
path <- "raw"
list.files(path)

## Get forward and reverse files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

## Sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Setting for filter reads: Forward and reverse
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## Filter reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# Estimate error forward and reverse
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

### De-replication and sample inference
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

## Set name into the dereplicated sample
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## Merge pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Convert to sequence tab
seqtab <- makeSequenceTable(mergers)

## Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

## Backup
saveRDS(seqtab.nochim, "seqtab_nochim.RDS")

## Reads tracking
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.csv(track, "reads_tracking.csv")

## Asign taxonomy and species
taxa <- assignTaxonomy(seqtab.nochim, "./tax/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "./tax/silva_species_assignment_v138.fa.gz")

saveRDS(taxa, "taxa_silva.RDS")

### DECOMTAM pipe
## Need inport file here
## experiment file has sample name (nothing after _R*) and sample_or_control 
samdf <- read.csv("experiment.csv", header = TRUE, row.names=1)

## Convert to phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

## Remove Mock community
ps <- prune_samples(sample_names(ps) != "Mock", ps)
saveRDS(ps, "ps_no_ChimeMock.RDS")

## Decontam package
## Contamination removal
contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
table(contamdf.freq$contaminant)






## convert sample name to dataframe
df <- as.data.frame(sample_data(ps))
## library size = nonchim values?
df$LibrarySize <- sample_sums(ps)
## order based on library size
df <- df[order(df$LibrarySize),]
## number of rows and start the count from 1 to nrow.
df$Index <- seq(nrow(df))

## Draw figures
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave("librarysize.png", width = 9, height = 6)

## Specific the control
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
## decontam prevalence method
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
## Show contaminations
table(contamdf.prev$contaminant)

## More strict method
#contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
#table(contamdf.prev05$contaminant)

## Figure drawing
## transform
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
## get negative samples
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
## get positive samples
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
## new dataframe
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("contamination.png", width = 9, height = 6)

## Remove contamination
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)

# extract otu table and coerce to data.frame
otu = microbiome::abundances(pSeq)
counts = tibble::as_tibble(otu)

# taxonomic classifications
taxonomy = as_tibble(pSeq@tax_table@.Data) %>% bind_cols(., ID = rownames(otu))

# compute total abundance, total prevalence, meanAbundanceWhenDetected
numbers = tibble::tibble(
  totalCounts = rowSums(otu), 
  prevalence = apply(otu, 1, function(ab){ sum(ab > 0) }),
  meanAbundanceWhenDetected = totalCounts / prevalence,
  minAbundanceWhenDetected = apply(otu, 1, function(ab){ abv = ab[ab > 0]; ifelse(length(abv) > 0, min(abv), NA) }),
  maxAbundance = apply(otu, 1, max)
)

# combine info
info = bind_cols(taxonomy, numbers)







# extract otu table and coerce to data.frame
otu = ps@otu_table@.Data
counts = tibble::as_tibble(otu)

taxonomy = as_tibble(ps@tax_table@.Data) %>% bind_cols(., ID = rownames(otu))

# compute total abundance, total prevalence, meanAbundanceWhenDetected
numbers = tibble::tibble(
  totalCounts = rowSums(otu), 
  prevalence = apply(otu, 1, function(ab){ sum(ab > 0) }),
  meanAbundanceWhenDetected = totalCounts / prevalence,
  minAbundanceWhenDetected = apply(otu, 1, function(ab){ abv = ab[ab > 0]; ifelse(length(abv) > 0, min(abv), NA) }),
  maxAbundance = apply(otu, 1, max)
)





dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
r <- plot_richness(ps, x="Duration", measures=c("Shannon", "Simpson"), color="When")
ggsave("figure5.png", width = 9, height = 6)





































