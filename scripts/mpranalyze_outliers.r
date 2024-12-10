#!/usr/bin/env Rscript

library("MPRAnalyze")
library(dplyr)
suppressPackageStartupMessages(library(arrow))
library(devtools)
library(parallel)
suppressPackageStartupMessages(load_all('BCalm'))


# Processing input arguments
args = commandArgs(trailingOnly=TRUE)
input = "data/simulations/variants/simulated_HepG2_variants.tsv.gz"
fraction = as.numeric(args[1])  
iteration = as.integer(args[2])
aggregate = as.logical(args[3])
nr_reps = 6
if (aggregate == TRUE) {
	output = paste("results/simulations/mpranalyze_aggr/outliers/fraction",gsub("\\.","", fraction),"/iteration",iteration,".feather",sep="")
} else {
	print("mpranalyze on barcodes")
	output = paste("results/simulations/mpranalyze_bc/outliers/fraction",gsub("\\.","", fraction),"/iteration",iteration,".feather",sep="")
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

	annot <- data.frame(name = grep("sample", colnames(dna), value=TRUE))
	annot$batch <- sub("sample([1-6])_.*", "\\1", annot$name)
	annot$condition <- sub("sample.*_(ref|alt)", "\\1", annot$name)

	annot$batch <- as.factor(annot$batch)
	annot$condition <- as.factor(annot$condition)
	
} else {
	dna <- create_dna_df(df)
	rna <- create_rna_df(df)

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
write_feather(res, output)
