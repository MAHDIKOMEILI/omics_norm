#######################
######Libraries########
#######################
library(qmap)
library(dplyr)
library(tidyr)
library(ggplot2)
# Define metadata columns and your genes of interest
non_gene_cols <- c("sampleID", "source")
genes_of_interest <- c("DRC3", "NEMP1", "GTSE1", "PCLAF", "DEPDC1", "PLK4", "STAT5A", "UBE2C")

# Separate reference and target samples
metabric_data <- merged_counts_final[merged_counts_final$source == "metabric", ]
target_data   <- merged_counts_final[merged_counts_final$source != "metabric", ]

# Get intersection of gene columns
common_genes <- intersect(
  setdiff(colnames(metabric_data), non_gene_cols),
  setdiff(colnames(target_data), non_gene_cols)
)

# Extract metadata and expression matrices
metabric_meta <- metabric_data[, non_gene_cols]
target_meta   <- target_data[, non_gene_cols]

metabric_expr <- metabric_data[, common_genes]
target_expr   <- target_data[, common_genes]

##############

# Create copy to hold normalized values
target_expr_qmap <- target_expr

for (gene in common_genes) {
  ref_vector <- metabric_expr[[gene]]
  mod_vector <- target_expr[[gene]]
  
  # Skip genes with insufficient variability
  if (length(unique(na.omit(ref_vector))) < 2 || length(unique(na.omit(mod_vector))) < 2) {
    cat("Skipping gene due to insufficient unique values:", gene, "\n")
    next
  }
  
  # Fit and apply quantile mapping
  fit_obj <- fitQmapQUANT(obs = ref_vector, mod = mod_vector, qstep = 0.01, nboot = 1, wet.day = FALSE)
  target_expr_qmap[[gene]] <- doQmapQUANT(mod_vector, fobj = fit_obj, type = "linear")
}

################

# Combine metadata and gene expression
metabric_final <- cbind(metabric_meta, metabric_expr)
target_final_qmap <- cbind(target_meta, target_expr_qmap)

# Merge into one dataset
combined_qmap <- bind_rows(metabric_final, target_final_qmap)

###############

df_long <- combined_qmap %>%
  dplyr::select(sampleID, source, all_of(genes_of_interest)) %>%
  pivot_longer(
    cols = all_of(genes_of_interest),
    names_to = "Gene",
    values_to = "Expression"
  )



ggplot(df_long, aes(x = Expression, color = source)) +
  geom_density(na.rm = TRUE) +
  facet_wrap(~ Gene, scales = "free") +
  theme_bw() +
  labs(title = "Quantile-Normalized Expression Distributions",
       x = "Expression", y = "Density")

###################

# Subset and prepare expression matrix
expr_pca <- combined_qmap %>% 
  dplyr::select(all_of(genes_of_interest)) %>%
  as.matrix()

# Run PCA
pca_res <- prcomp(expr_pca, center = TRUE, scale. = TRUE)

# Get explained variance
pca_var <- round((pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100, 1)

# Construct PCA plot data
pca_df <- combined_qmap %>%
  dplyr::select(sampleID, source) %>%
  mutate(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = source)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw() +
  labs(
    title = "PCA of Quantile-Normalized Expression (Genes of Interest)",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  )
