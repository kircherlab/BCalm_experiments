#!/usr/bin/env Rscript

library("MPRAnalyze")
library(dplyr)
suppressPackageStartupMessages(library(arrow))
library(devtools)
library(parallel)
suppressPackageStartupMessages(load_all('BCalm'))


# Processing input arguments
args = commandArgs(trailingOnly=TRUE)
input = args[1] #"data/simulations/variants/simulated_HepG2_first_100_variants.tsv.gz"
print(input)
nr_reps = 6  
aggregate = as.logical(args[2])
num_cores <- as.integer(args[3])
if (aggregate == TRUE) {
	output = "results/simulations/mpranalyze_aggr/simulated_6reps"
} else {
	print("mpranalyze on barcodes")
	output = "results/simulations/mpranalyze_bc/simulated_6reps"
}

df <- read.table(file=gzfile(input), sep='\t', header=TRUE)
df$variant_id <- df$name


# Code for running mpranalyze on a single combination
run_mpranalyze <- function() {
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

		annot <- data.frame(name = grep("sample", colnames(dna), value=TRUE))
		annot$batch <- sub("sample([1-6])_.*", "\\1", annot$name)
		annot$condition <- sub("sample.*_(ref|alt)", "\\1", annot$name)

		annot$batch <- as.factor(annot$batch)
		annot$condition <- as.factor(annot$condition)
		
	} else {
		dna <- create_dna_df(df_subset)
		rna <- create_rna_df(df_subset)

		annot <- data.frame(name = grep("count", colnames(dna), value=TRUE))
		annot$batch <- sub("sample_count_([1-6]+)_.*", "\\1", annot$name)
		annot$condition <- sub("sample_count_.*_(ref|alt)", "\\1", annot$name)
		annot$barcode_allelic <- sub("sample_count_.*_(bc.*)", "\\1", annot$name)
			
		# make sure mpranalyze doesn't treat these as numeric, to prevent the error: L-BFGS-B needs finite values of 'fn'
		annot$batch <- as.factor(annot$batch)
		annot$condition <- as.factor(annot$condition)
		annot$barcode_allelic <- as.factor(annot$barcode_allelic)
	}
	
	rna <- as.matrix(rna)
	dna <- as.matrix(dna)

	rownames(annot) <- annot$name
	annot$name <- NULL

	mpra <- MpraObject(dnaCounts = dna, rnaCounts = rna, dnaAnnot=annot, rnaAnnot=annot)

	# upper quartile normalization
	print("UQ normalization")
	mpra <- estimateDepthFactors(mpra, lib.factor = c("batch"),
								which.lib = "dna", 
								depth.estimator = "uq")
	mpra <- estimateDepthFactors(mpra, lib.factor = c("batch"),
								which.lib = "rna", 
								depth.estimator = "uq")

	print("Running MPRAnalyze")

	if (aggregate == TRUE) {
		fit <- analyzeComparative(obj = mpra, 
							dnaDesign = ~ batch + condition, 
							rnaDesign = ~ condition, 
							reducedDesign = ~ 1)
	} else {
		fit <- analyzeComparative(obj = mpra, 
							dnaDesign = ~ batch + condition + barcode_allelic, 
							rnaDesign = ~ condition, 
							reducedDesign = ~ 1)
	}




	res <- testLrt(fit)
	res$variant_id <- row.names(res)
	write_feather(res, paste(output, ".feather", sep=""))
}

run_mpranalyze()

print("Done!")
