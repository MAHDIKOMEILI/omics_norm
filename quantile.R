library(dplyr)
library(ggplot2)

# Make a copy of merged_counts_final so we don't overwrite it
merged_counts_final_fullnorm <- merged_counts_final

non_gene_cols <- c("sampleID", "source")  # plus any other columns you want to exclude
all_genes <- setdiff(colnames(merged_counts_final_fullnorm), non_gene_cols)

metabric_data <- merged_counts_final_fullnorm[merged_counts_final_fullnorm$source == "metabric", ]

adjust_gene <- function(target_values, ref_values, nquantiles = 10) {
  probs <- seq(0, 1, length.out = nquantiles + 1)
  breaks_target <- quantile(target_values, probs = probs, na.rm = TRUE)
  breaks_ref    <- quantile(ref_values,    probs = probs, na.rm = TRUE)
  
  adjusted <- target_values
  
  for (i in 1:nquantiles) {
    if (i < nquantiles) {
      idx <- which(target_values >= breaks_target[i] & target_values < breaks_target[i + 1])
      idx_ref <- which(ref_values    >= breaks_ref[i]    & ref_values    < breaks_ref[i + 1])
    } else {
      idx <- which(target_values >= breaks_target[i] & target_values <= breaks_target[i + 1])
      idx_ref <- which(ref_values    >= breaks_ref[i]    & ref_values    <= breaks_ref[i + 1])
    }
    
    if (length(idx) > 0) {
      mean_target <- mean(target_values[idx], na.rm = TRUE)
      mean_ref <- if (length(idx_ref) == 0) {
        mean(ref_values, na.rm = TRUE)
      } else {
        mean(ref_values[idx_ref], na.rm = TRUE)
      }
      if (mean_target == 0) mean_target <- 1e-6
      adjusted[idx] <- target_values[idx] / mean_target * mean_ref
    }
  }
  return(adjusted)
}

# Identify the non-reference sources
non_metabric_sources <- setdiff(unique(merged_counts_final_fullnorm$source), "metabric")

for (src in non_metabric_sources) {
  cat("Adjusting dataset:", src, "\n")
  
  # Indices for samples in this source
  target_idx <- which(merged_counts_final_fullnorm$source == src)
  
  for (gene in all_genes) {
    # Get target (non-metabric) values
    target_vals <- as.numeric(merged_counts_final_fullnorm[target_idx, gene])
    # Get reference (metabric) values
    ref_vals <- as.numeric(metabric_data[[gene]])
    
    # If there are NAs, you may want to impute or skip
    if (any(is.na(target_vals)) || any(is.na(ref_vals))) {
      cat("Warning: Column", gene, "in dataset", src, "contains NA values.\n")
      # Optionally handle NAs (impute or skip)
    }
    
    # Adjust the target values
    adjusted_vals <- adjust_gene(target_vals, ref_vals, nquantiles = 10)
    
    # Store back in the new object
    merged_counts_final_fullnorm[target_idx, gene] <- adjusted_vals
  }
}

##################PCAPLOT unnormalised

# --- STEP 1: Define which columns represent gene expression
# Assume merged_counts_final has columns like sampleID, source, and gene expression columns.
non_gene_cols <- c("sampleID", "source")
gene_cols <- setdiff(colnames(merged_counts_final), non_gene_cols)

# --- STEP 2: Split the data into metabric and non-metabric subsets
metabric_data <- merged_counts_final %>%
  filter(source == "metabric") %>%
  select(all_of(gene_cols))

non_metabric_data <- merged_counts_final %>%
  filter(source != "metabric") %>%
  select(all_of(gene_cols))

# --- STEP 3: Identify gene columns that have complete (non-NA) data in BOTH groups
# For metabric samples, get the names of genes with no NA values:
complete_met_genes <- colnames(metabric_data)[apply(metabric_data, 2, function(x) all(!is.na(x)))]
# For non-metabric samples, do the same:
complete_non_met_genes <- colnames(non_metabric_data)[apply(non_metabric_data, 2, function(x) all(!is.na(x)))]
# Take the intersection of these gene sets:
common_expr_cols <- intersect(complete_met_genes, complete_non_met_genes)
cat("Number of gene columns common (complete) in both groups:", length(common_expr_cols), "\n")

# --- STEP 4: Create a data frame for PCA using the intersection of complete gene columns
# We keep sampleID and source columns from the original data.
df_for_pca <- merged_counts_final %>%
  select(sampleID, source, all_of(common_expr_cols))

# --- STEP 5: Convert the expression portion to a numeric matrix
# We assume that after removing sampleID and source, the rest are expression values.
expr_matrix <- as.data.frame(lapply(df_for_pca %>% select(-sampleID, -source), as.numeric))
rownames(expr_matrix) <- df_for_pca$sampleID

# --- STEP 6: Run PCA on the intersection expression matrix
pca_result <- prcomp(expr_matrix, center = TRUE, scale. = TRUE)

# --- STEP 7: Create a plotting data frame with the PCA results
pca_df <- data.frame(
  sampleID = df_for_pca$sampleID,
  source = df_for_pca$source,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# --- STEP 8: Plot the PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = source)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(
    title = "PCA Plot Using Intersection of Complete Gene Columns",
    x = "PC1",
    y = "PC2"
  )
##############################PCA PLOT normalised###################

# --- STEP 1: Identify Gene Columns in Normalized Data ---

# Specify columns that are not gene expression (identifiers, etc.)
non_gene_cols <- c("sampleID", "source")

# Gene columns: all columns except the non-gene ones
gene_cols <- setdiff(colnames(merged_counts_final_fullnorm), non_gene_cols)

# --- STEP 2: Split Normalized Data into Metabric and Non-Metabric Subsets ---

metabric_data_norm <- merged_counts_final_fullnorm %>%
  filter(source == "metabric") %>%
  select(all_of(gene_cols))

non_metabric_data_norm <- merged_counts_final_fullnorm %>%
  filter(source != "metabric") %>%
  select(all_of(gene_cols))

# --- STEP 3: Find the Intersection of Complete Gene Columns ---

# Identify gene columns with no NA values in metabric data
complete_met_genes_norm <- colnames(metabric_data_norm)[
  apply(metabric_data_norm, 2, function(x) all(!is.na(x)))
]

# Identify gene columns with no NA values in non-metabric data
complete_non_met_genes_norm <- colnames(non_metabric_data_norm)[
  apply(non_metabric_data_norm, 2, function(x) all(!is.na(x)))
]

# Take the intersection of these two sets:
common_expr_cols_norm <- intersect(complete_met_genes_norm, complete_non_met_genes_norm)
cat("Number of complete gene columns in normalized data:", length(common_expr_cols_norm), "\n")

# --- STEP 4: Prepare Data for PCA ---

# Create a data frame for PCA that includes sampleID, source, and the common gene columns.
df_for_pca_norm <- merged_counts_final_fullnorm %>%
  select(sampleID, source, all_of(common_expr_cols_norm))

# Convert the gene expression portion to a numeric data frame/matrix
expr_matrix_norm <- as.data.frame(lapply(df_for_pca_norm %>% select(-sampleID, -source), as.numeric))
rownames(expr_matrix_norm) <- df_for_pca_norm$sampleID

# --- STEP 5: Run PCA on the Normalized Data ---

pca_result_norm <- prcomp(expr_matrix_norm, center = TRUE, scale. = TRUE)

# --- STEP 6: Create a Data Frame for Plotting PCA Results ---

pca_df_norm <- data.frame(
  sampleID = df_for_pca_norm$sampleID,
  source = df_for_pca_norm$source,
  PC1 = pca_result_norm$x[, 1],
  PC2 = pca_result_norm$x[, 2]
)

# --- STEP 7: Plot the PCA for Normalized Data ---

ggplot(pca_df_norm, aes(x = PC1, y = PC2, color = source)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(
    title = "PCA Plot (Normalized Data)",
    x = "Principal Component 1",
    y = "Principal Component 2"
  )
# Suppose your PCA object is pca_result
pca_var <- pca_result_norm$sdev^2                 # Variances of the PCs
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 1)  # Percentage of total variance

ggplot(pca_df_norm, aes(x = PC1, y = PC2, color = source)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "PCA Plot (Normalized Data)",
    x = paste0("PC1 (", pca_var_percent[1], "%)"),
    y = paste0("PC2 (", pca_var_percent[2], "%)")
  ) +
  theme_bw()

