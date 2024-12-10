#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(BCalm))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))

args = commandArgs(trailingOnly=TRUE)
input = "data/simulations/variants/simulated_HepG2_variants.tsv.gz"
fraction = as.numeric(args[1])  
iteration = as.integer(args[2])
aggregate = as.logical(args[3])
nr_reps = 6
if (aggregate == TRUE) {
	output = paste("results/simulations/mpralm/outliers/fraction",gsub("\\.","", fraction),"/iteration",iteration,".feather",sep="")
} else {
	output = paste("results/simulations/BCalm/outliers/fraction",gsub("\\.","", fraction),"/iteration",iteration,".feather",sep="")
}

if (!dir.exists(dirname(output))) {
  dir.create(dirname(output), recursive = TRUE)
}

df <- read.table(file=gzfile(input), sep='\t', header=TRUE)
df$variant_id <- df$name

num_outliers <- ceiling(nrow(df) * fraction)
outlier_indices <- sample(1:nrow(df), num_outliers)
num_selected_columns <- sample(1:nr_reps, 1, prob = c(0.85, 0.05, 0.04, 0.03, 0.02, 0.01))
selected_columns <- sample(grep("rna", colnames(df)), num_selected_columns)
df[outlier_indices,selected_columns] <- df[outlier_indices,selected_columns] * 25  


if (aggregate == TRUE) {
	df <- df %>% select(variant_id, allele, matches("count")) %>% as.data.frame()
	df_summed <- df %>% 
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
	dna <- create_dna_df(df)
	rna <- create_rna_df(df)
	# number of barcodes is number of columns, divided by (nr of samples * number of alleles)
	bcs <- ncol(dna) / (nr_reps)
	# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector. 
	block_vector <- rep(1:nr_reps, each=bcs)
}


mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra)))
mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)

# Finding significant variants
toptab_allele <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele$variant_id <- row.names(toptab_allele)
write_feather(toptab_allele, output)





