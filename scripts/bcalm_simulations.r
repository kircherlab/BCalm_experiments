#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(BCalm))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))

args = commandArgs(trailingOnly=TRUE)
input = "data/simulations/variants/simulated_HepG2_variants.tsv.gz"
nr_reps = as.integer(args[1])  
aggregate = as.logical(args[2])
if (aggregate == TRUE) {
	output = "results/simulations/mpralm/"
} else {
	print("BCalm  on barcodes")
	output = "results/simulations/BCalm/"
}


df <- read.table(file=gzfile(input), sep='\t', header=TRUE)
df$variant_id <- df$name

df_subset <- df %>% select(variant_id, allele, matches("count")) %>% as.data.frame()

if (aggregate == TRUE) {
	df_summed <- df_subset %>% 
		group_by(variant_id, allele) %>% 
		summarise(across(everything(), sum),.groups = 'drop') %>%
		as.data.frame()

	rna <- df_summed[, grep("variant_id|allele|rna", colnames(df_summed), ignore.case=TRUE)]
	dna <- df_summed[, grep("variant_id|allele|dna", colnames(df_summed), ignore.case=TRUE)]

	rna <- cbind(filter(rna,allele=="ref")[,-2], filter(rna,allele=="alt")[,3:(nr_reps+2)])
	dna <- cbind(filter(dna,allele=="ref")[,-2], filter(dna,allele=="alt")[,3:(nr_reps+2)])

	names(rna)[-1] <- paste0("sample", c(1:nr_reps), "_", rep(c("ref", "alt"), each = nr_reps))
	row.names(rna) <- rna$variant_id
	rna$variant_id <- NULL

	names(dna)[-1] <- paste0("sample", c(1:nr_reps), "_", rep(c("ref", "alt"), each = nr_reps))
	row.names(dna) <- dna$variant_id
	dna$variant_id <- NULL

	block_vector <- rep(1:nr_reps, 2)

	
} else {
	dna <- create_dna_df(df_subset)
	rna <- create_rna_df(df_subset)
	# number of barcodes is number of columns, divided by (nr of samples * number of alleles)
	bcs <- ncol(dna) / (nr_reps)
	# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector. 
	block_vector <- rep(1:nr_reps, each=bcs)
}


mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra)))

mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)

# Finding significant variants
toptab_allele <- topTreat(mpralm_fit, coef = 2, number = Inf)
toptab_allele$variant_id <- row.names(toptab_allele)

write_feather(toptab_allele, paste(output, ".feather", sep=""))




