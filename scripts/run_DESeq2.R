#!/usr/bin/env Rscript

# Load required libraries
library(DESeq2)
library(apeglm)
library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript deseq2_analysis.R <count_matrix> <coldata_file> <output_dir>")
}

# Assign arguments
count_matrix_file <- args[1]
coldata_file <- args[2]
output_dir <- args[3]

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load and Prepare Count Matrix
message("Loading count matrix...")
merged_data <- read_tsv(count_matrix_file, col_names = TRUE)
cts <- as.matrix(merged_data[,-1])
rownames(cts) <- merged_data[[1]]

# 2. Load Coldata
message("Loading column data...")
coldata <- read_csv(coldata_file, col_names = FALSE)
colnames(coldata) <- c("sample", "condition", "type")
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# Ensure column names match
stopifnot(all(colnames(cts) %in% rownames(coldata)))
cts <- cts[, rownames(coldata)]

# 3. Create DESeqDataSet
message("Creating DESeqDataSet...")
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ condition
)

# 4. Filtering and Preprocessing
message("Filtering low-count genes...")
# Keep genes with at least 10 reads in at least 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Set reference level
dds$condition <- relevel(dds$condition, ref = "control")

# 5. Run DESeq
message("Running DESeq analysis...")
dds <- DESeq(dds)

# 6. Extract and Save Results
message("Extracting results...")
condition_levels <- levels(dds$condition)

# Loop through conditions (excluding reference)
for (i in seq(2, length(condition_levels))) {
    condition_name <- condition_levels[i]
    coef_name <- paste0("condition_", condition_name, "_vs_control")
    
    # Shrink log fold change
    resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    
    # Convert to data frame and add gene column
    res_df <- as.data.frame(resLFC) %>%
        rownames_to_column("gene")
    
    # Create output file path
    output_file <- file.path(output_dir, paste0(coef_name, ".tsv"))
    
    # Write results
    write_tsv(res_df, output_file)
    message(paste("Saved results for", condition_name, "to", output_file))
}

message("DESeq2 analysis complete!")