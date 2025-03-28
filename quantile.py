import numpy as np

# Load the TSV file and treat "NA" as missing (converted to np.nan)
data = np.genfromtxt(
    "counts.tsv",
    delimiter="\t",
    names=True,
    dtype=None,
    encoding="utf-8",
    missing_values="NA",
    filling_values=np.nan
)

# Extract sample IDs and source labels
sample_ids = data["sampleID"]
sources = data["source"]

# Identify gene expression columns (all columns except "sampleID" and "source")
all_columns = data.dtype.names
gene_cols = [col for col in all_columns if col not in ("sampleID", "source")]

# Build a 2D NumPy array for gene expression data (samples x genes)
X = np.column_stack([data[col] for col in gene_cols])

# Optional: Print shapes to verify the import.
print("Number of samples:", sample_ids.shape[0])
print("Number of gene columns:", X.shape[1])


# Optional: Print shapes to verify the import.
print("Number of samples:", sample_ids.shape[0])
print("Number of gene columns:", X.shape[1])

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# ASSUMPTIONS:
# sample_ids: 1D array-like of sample IDs (strings)
# sources: 1D array-like of source labels (strings)
# X: 2D array (n_samples x n_genes) of gene expression values (floats)
#
# For example, you might have loaded your data from a CSV file and extracted:
#   sample_ids = np.array([...])
#   sources = np.array([...])
#   X = np.array([...])
#
# All computations below use only NumPy for data manipulation.
# =============================================================================

# For safety, make a copy of the expression matrix to apply normalization
X_fullnorm = X.copy()

# Identify indices corresponding to the metabric samples (reference)
metabric_idx = np.where(sources == "metabric")[0]

# -----------------------------------------------------------------------------
# Define a custom quantile normalization function similar to your R adjust_gene
# -----------------------------------------------------------------------------
def adjust_gene(target_values, ref_values, nquantiles=10):
    """
    Adjust target gene expression values by aligning their quantiles to those
    of the reference values.
    
    Parameters:
        target_values (np.array): 1D array for the target gene expression.
        ref_values (np.array): 1D array for the reference (metabric) gene expression.
        nquantiles (int): Number of quantile intervals.
        
    Returns:
        np.array: Adjusted target values.
    """
    probs = np.linspace(0, 1, nquantiles + 1)
    # Compute quantile breakpoints, using nanquantile in case of missing values
    breaks_target = np.nanquantile(target_values, probs)
    breaks_ref = np.nanquantile(ref_values, probs)
    
    adjusted = target_values.copy()
    
    for i in range(nquantiles):
        if i < nquantiles - 1:
            # For non-final quantiles: half-open interval [break_i, break_i+1)
            idx = np.where((target_values >= breaks_target[i]) & (target_values < breaks_target[i+1]))[0]
            idx_ref = np.where((ref_values >= breaks_ref[i]) & (ref_values < breaks_ref[i+1]))[0]
        else:
            # Last quantile: include the upper bound
            idx = np.where((target_values >= breaks_target[i]) & (target_values <= breaks_target[i+1]))[0]
            idx_ref = np.where((ref_values >= breaks_ref[i]) & (ref_values <= breaks_ref[i+1]))[0]
        
        if idx.size > 0:
            mean_target = np.nanmean(target_values[idx])
            # If no samples fall in the reference interval, use overall reference mean
            if idx_ref.size == 0:
                mean_ref = np.nanmean(ref_values)
            else:
                mean_ref = np.nanmean(ref_values[idx_ref])
            # Prevent division by zero
            if mean_target == 0:
                mean_target = 1e-6
            adjusted[idx] = target_values[idx] / mean_target * mean_ref
    return adjusted

# -----------------------------------------------------------------------------
# Apply quantile normalization for each non-metabric source and each gene
# -----------------------------------------------------------------------------
# Get unique sources and exclude "metabric"
unique_sources = np.unique(sources)
non_metabric_sources = unique_sources[unique_sources != "metabric"]

n_samples, n_genes = X.shape

# Apply quantile normalization for each non-metabric source and each gene
for src in non_metabric_sources:
    print("Adjusting dataset:", src)
    # Find indices of samples with the current source
    target_idx = np.where(sources == src)[0]
    for gene in range(n_genes):
        # Get target gene values and reference (metabric) values as floats
        target_vals = X_fullnorm[target_idx, gene].astype(float)
        ref_vals = X_fullnorm[metabric_idx, gene].astype(float)
        
        # If there are any NA values, skip processing this gene
        if np.any(np.isnan(target_vals)) or np.any(np.isnan(ref_vals)):
            print(f"Warning: Gene {gene} in dataset {src} contains NA values. Skipping normalization for this gene.")
            continue
        
        # Adjust the target values using the custom quantile normalization function
        adjusted_vals = adjust_gene(target_vals, ref_vals, nquantiles=10)
        X_fullnorm[target_idx, gene] = adjusted_vals

# =============================================================================
# PCA on Normalized Data using only NumPy (no pandas)
# =============================================================================

# Identify gene columns that have complete (non-NA) data in both groups.
# Create boolean masks for complete data in metabric and non-metabric groups.
complete_met = ~np.isnan(X_fullnorm[metabric_idx, :]).any(axis=0)
non_metabric_idx = np.where(sources != "metabric")[0]
complete_non_met = ~np.isnan(X_fullnorm[non_metabric_idx, :]).any(axis=0)

# Select columns (gene indices) that are complete in both groups.
common_expr_cols_norm = np.where(complete_met & complete_non_met)[0]
print("Number of complete gene columns in normalized data:", len(common_expr_cols_norm))

# Prepare the expression matrix for PCA (only the complete gene columns)
X_pca = X_fullnorm[:, common_expr_cols_norm]

# Standardize the data (center and scale) similar to prcomp(center=TRUE, scale.=TRUE)
mean_vals = np.mean(X_pca, axis=0)
std_vals = np.std(X_pca, axis=0, ddof=1)
# Avoid division by zero in case of zero standard deviation
std_vals[std_vals == 0] = 1
X_scaled = (X_pca - mean_vals) / std_vals

# Compute PCA via Singular Value Decomposition (SVD)
U, S, Vt = np.linalg.svd(X_scaled, full_matrices=False)
# The principal components are given by U * S.
PCs = U * S

# Compute percentage variance explained
explained_variance = (S ** 2) / (X_scaled.shape[0] - 1)
explained_variance_ratio = explained_variance / explained_variance.sum()
pc1_label = f"PC1 ({explained_variance_ratio[0]*100:.1f}%)"
pc2_label = f"PC2 ({explained_variance_ratio[1]*100:.1f}%)"

# =============================================================================
# Plot the PCA results using matplotlib
# =============================================================================
plt.figure(figsize=(8, 6))
# Iterate through each unique source and plot its samples
for src in np.unique(sources):
    idx = np.where(sources == src)[0]
    plt.scatter(PCs[idx, 0], PCs[idx, 1], label=src, alpha=0.7)
plt.xlabel(pc1_label)
plt.ylabel(pc2_label)
plt.title("PCA Plot (Normalized Data)")
plt.legend()
plt.tight_layout()
plt.show()
