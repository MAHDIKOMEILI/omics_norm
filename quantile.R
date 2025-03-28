#########################
#######Libraries#########
#########################

library(dplyr)
library(readr)


#############################
#######set directory#########
#############################

# Set your directory path
data_dir <- "~/er_positive/MAS5"

# List all CSV files in the directory
files <- list.files(data_dir, pattern = "\\.tsv$", full.names = TRUE)

###########################################
#########Separate clin and counts##########
###########################################

for (file in files) {
  # 1. Get a short name for new data frames
  var_name <- tools::file_path_sans_ext(basename(file))
  
  # 2. Read the TSV
  dataset <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 3. Identify all columns that start with "clin_"
  clin_cols <- grep("^clin_", colnames(dataset), value = TRUE)
  
  # 4. Create the clinical subset: these are the columns that start with "clin_"
  #    (this includes "clin_sample_name" automatically, if it exists)
  dataset_clin <- dataset[, clin_cols, drop = FALSE]
  
  # 5. Create the counts subset: everything that doesn't start with "clin_"
  #    BUT we also force "clin_sample_name" to stay in the counts so we can track sample IDs
  non_clin_cols <- setdiff(colnames(dataset), clin_cols)
  
  # Force "clin_sample_name" into the counts subset if it exists in the data
  if ("clin_sample_name" %in% colnames(dataset)) {
    non_clin_cols <- union("clin_sample_name", non_clin_cols)
  }
  
  rm(dataset,dataset_clin,dataset_counts,metabric_clin,metabric_counts)
  
  dataset_counts <- dataset[, non_clin_cols, drop = FALSE]
  
  # 6. Assign the subsets to variables named <fileName>_clin and <fileName>_counts
  assign(paste0(var_name, "_clin"), dataset_clin)
  assign(paste0(var_name, "_counts"), dataset_counts)
}

# Loop over objects in the global environment with names ending in _clin or _counts
for (obj_name in ls(pattern = "(_clin$|_counts$)")) {
  obj <- get(obj_name)
  # Check if the object is a data frame and has the target column
  if (is.data.frame(obj) && "clin_sample_name" %in% colnames(obj)) {
    colnames(obj)[colnames(obj) == "clin_sample_name"] <- "sampleID"
    assign(obj_name, obj)  # Update the object in the workspace
  }
}

###################################################
##############MAS MERGED COUNTS####################
###################################################

# Exclude metabric_counts from the counts data frames:
counts_names <- setdiff(ls(pattern = "_counts$"), "metabric_counts")
counts_list <- mget(counts_names)

# Coerce each object to a data frame (if it isn't already)
counts_list <- lapply(counts_list, function(x) {
  if (!is.data.frame(x)) as.data.frame(x) else x
})

# Identify common columns (features) across these counts data frames
common_features <- Reduce(intersect, lapply(counts_list, colnames))
print(common_features)  # Check which features are common

# Subset each counts data frame to include only the common columns
counts_list_subset <- lapply(counts_list, function(df) {
  df <- as.data.frame(df)  # ensure df is a data frame
  df[, common_features, drop = FALSE]
})

# Merge the counts data frames by row-binding them together and store in mas_merged_counts
mas_merged_counts <- do.call(rbind, counts_list_subset)

# Set the 'clin_sample_name' column as row names
rownames(mas_merged_counts) <- mas_merged_counts$sampleID

# (Optional) Verify the row names
head(rownames(mas_merged_counts))


####################################################
###################METABRIC#########################
####################################################

clinical_cols <- c("X.Patient.Identifier", "ER.status.measured.by.IHC", "HER2.status.measured.by.SNP6", "Sex", "Age.at.Diagnosis", "Overall.Survival..Months.", "Overall.Survival.Status", "Patient.s.Vital.Status", "DFS", "Event")

# Path to your metabric file
file_path <- "~/er_positive/MAS5/metabric.tsv"

# Read the data
metabric_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)


# 3) Subset the clinical data
metabric_clin <- metabric_data[, clinical_cols, drop = FALSE]

# 4) Subset the counts data: all columns not listed as clinical.
#    But if you want to keep the sampleID column in the counts data as well, force it in.
counts_cols <- setdiff(colnames(metabric_data), clinical_cols)
if ("X.Patient.Identifier" %in% colnames(metabric_data)) {
  counts_cols <- union("X.Patient.Identifier", counts_cols)
}
metabric_counts <- metabric_data[, counts_cols, drop = FALSE]

# 5) (Optional) Set the sampleID column as the row names for the counts data,
#    while still keeping the sampleID column in the data frame.
if ("X.Patient.Identifier" %in% colnames(metabric_counts)) {
  rownames(metabric_counts) <- metabric_counts$X.Patient.Identifier
  # If you later decide to remove the column, uncomment the next line:
  # metabric_counts$sampleID <- NULL
}

# Now you have:
#   - metabric_clin: containing the clinical information columns
#   - metabric_counts: containing the counts data plus the sampleID column

colnames(metabric_counts)

# Define the names of the columns you want to remove
remove_cols <- c("Lymph.nodes.examined.positive", "Nottingham.prognostic.index", "Cellularity", "Chemotherapy", "Cohort", "Hormone.Therapy", "Inferred.Menopausal.State", "Integrative.Cluster", "Pam50...Claudin.low.subtype", "X3.Gene.classifier.subtype", "Primary.Tumor.Laterality", "Radio.Therapy", "Tumor.Other.Histologic.Subtype", "Type.of.Breast.Surgery")  # Replace with your actual column names

# Subset metabric_counts to exclude these columns
metabric_counts <- metabric_counts[, !(colnames(metabric_counts) %in% remove_cols)]

colnames(metabric_counts)[colnames(metabric_counts) == "X.Patient.Identifier"] <- "sampleID"
colnames(metabric_clin)[colnames(metabric_clin) == "X.Patient.Identifier"] <- "sampleID"

##############################
######Merge counts data#######
##############################

# ---------- Step 1: Intersection Merge ----------

# Compute common (intersecting) columns between mas_merged_counts and metabric_counts

common_features_final <- intersect(colnames(mas_merged_counts), colnames(metabric_counts))

# Create an intersection merge: subset each object to the common columns, then row-bind

merged_counts_intersect <- rbind(
  mas_merged_counts[, common_features_final, drop = FALSE],
  metabric_counts[, common_features_final, drop = FALSE]
)

# ---------- Step 2: Identify Extra (Uncommon) Columns ----------

# Extra columns present in mas_merged_counts but not in the common set

extra_mas <- setdiff(colnames(mas_merged_counts), common_features_final)

# Extra columns present in metabric_counts but not in the common set

extra_met <- setdiff(colnames(metabric_counts), common_features_final)

# The union of extra columns from both datasets

extra_all <- union(extra_mas, extra_met)

"KIAA0101" %in% extra_all

# ---------- Step 3: Patch Extra Columns into the Merged Data ----------

# Add each extra column to the merged_counts_intersect object, initializing with NA

for (col in extra_all) {
  merged_counts_intersect[[col]] <- NA
}

# Define a helper function to patch extra column values from an original dataset

patch_extra_columns <- function(original_df, merged_df, extra_cols) {
  # Use the correct unique identifier "sampleID"
  sample_ids <- original_df$sampleID
  idx <- match(sample_ids, merged_df$sampleID)
  for (col in extra_cols) {
    if (col %in% colnames(original_df)) {
      merged_df[idx, col] <- original_df[[col]]
    }
  }
  return(merged_df)
}

# Patch values for extra columns from mas_merged_counts and metabric_counts respectively:

merged_counts_final <- merged_counts_intersect
merged_counts_final <- patch_extra_columns(mas_merged_counts, merged_counts_final, extra_mas)
merged_counts_final <- patch_extra_columns(metabric_counts, merged_counts_final, extra_met)

# ---------- Final Output ----------
# merged_counts_final now contains:
#  - All rows (samples) from both mas_merged_counts and metabric_counts.
#  - The common columns merged from both datasets.
#  - Extra columns from either dataset added with NA values where a sample didn't have that information.
# The unique sample identifier is in 'clin_sample_name'.
# Optionally, you can set row names based on 'clin_sample_name':

rownames(merged_counts_final) <- merged_counts_final$sampleID

#################################################
######Merging Alternative Gene name columns######
#################################################

anyNA(merged_counts_final$PCLAF)

merged_counts_final <- merged_counts_final %>%
  mutate(
    DRC3 = coalesce(DRC3, LRRC48),
    NEMP1 = coalesce(NEMP1, TMEM194A),
    PCLAF = coalesce(PCLAF, KIAA0101)
  ) %>%
  dplyr::select(-LRRC48, -TMEM194A, -KIAA0101)

###########################
##########Tagging##########
###########################

# Initialize a new "source" column with NA in the merged_counts_final data frame.

merged_counts_final$source <- NA

# Tag each sample by checking membership in the respective clinical dataset.
merged_counts_final$source[ merged_counts_final$sampleID %in% GSE25066_primary_clin$sampleID ] <- "GSE25066"
merged_counts_final$source[ merged_counts_final$sampleID %in% MAINZ_primary_clin$sampleID ]    <- "MAINZ"
merged_counts_final$source[ merged_counts_final$sampleID %in% STK_primary_clin$sampleID ]      <- "STK"
merged_counts_final$source[ merged_counts_final$sampleID %in% TRANSBIG_primary_clin$sampleID ] <- "TRANSBIG"
merged_counts_final$source[ merged_counts_final$sampleID %in% MSK_primary_clin$sampleID ]      <- "MSK"
merged_counts_final$source[ merged_counts_final$sampleID %in% UPP_primary_clin$sampleID ]      <- "UPP"
merged_counts_final$source[ merged_counts_final$sampleID %in% VDX_primary_clin$sampleID ]      <- "VDX"
merged_counts_final$source[ merged_counts_final$sampleID %in% metabric_clin$sampleID ]         <- "metabric"

# Check the assignment:

table(merged_counts_final$source)

##########################
######Optional Save#######
##########################

write.table(merged_counts_final, file = "counts_unnorm.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

###########################################
########Clinical data cleaning#############
###########################################

# Find all objects in the environment that end with _ER_Positive_clin
clin_datasets <- ls(pattern = ".*_ER_Positive_clin$")

# Loop over each dataset
for (dataset in clin_datasets) {
  # Retrieve the clinical dataset
  clin_data <- get(dataset)
  
  # Generate a new object name by replacing _ER_Positive_clin with _primary_clin
  new_name <- sub("_ER_Positive_clin$", "_primary_clin", dataset)
  
  # Rename the specified columns
  names(clin_data)[names(clin_data) == "clin_days_to_tumor_recurrence"] <- "dfs"
  names(clin_data)[names(clin_data) == "clin_recurrence_status"]        <- "dfs_status"
  names(clin_data)[names(clin_data) == "clin_days_to_death"]              <- "os"
  names(clin_data)[names(clin_data) == "clin_vital_status"]               <- "os_status"
  names(clin_data)[names(clin_data) == "clin_dmfs_days"]                  <- "dmfs"
  names(clin_data)[names(clin_data) == "clin_dmfs_status"]                <- "dmfs_status"
  names(clin_data)[names(clin_data) == "clin_er"]                <- "ER"
  names(clin_data)[names(clin_data) == "clin_her2"]                <- "her2"
  names(clin_data)[names(clin_data) == "clin_age_at_initial_pathologic_diagnosis"]                <- "age"
  
  # Select only the renamed columns
  desired_cols <- c("dfs", "dfs_status", "os", "os_status", "dmfs", "dmfs_status", "ER", "her2", "age", "sampleID")
  clin_data <- clin_data[, desired_cols, drop = FALSE]
  
  # Assign the modified dataset with the new name in the global environment
  assign(new_name, clin_data)
}

###############

# Rename the columns in metabric_clin as specified:
names(metabric_clin)[names(metabric_clin) == "X.Patient.Identifier"]       <- "sampleID"
names(metabric_clin)[names(metabric_clin) == "Overall.Survival..Months."]    <- "os"
names(metabric_clin)[names(metabric_clin) == "Overall.Survival.Status"]      <- "os_status"
names(metabric_clin)[names(metabric_clin) == "Patient.s.Vital.Status"]       <- "vital_status"
names(metabric_clin)[names(metabric_clin) == "DFS"]                          <- "dfs"
names(metabric_clin)[names(metabric_clin) == "DFS_event"]                    <- "dfs_status"
names(metabric_clin)[names(metabric_clin) == "ER.status.measured.by.IHC"]      <- "ER"
names(metabric_clin)[names(metabric_clin) == "HER2.status.measured.by.SNP6"]    <- "her2"
names(metabric_clin)[names(metabric_clin) == "Sex"]                          <- "sex"
names(metabric_clin)[names(metabric_clin) == "Age.at.Diagnosis"]             <- "age"

#################################################
###########Combine non metabric_clin#############
#################################################

clin_meta <- bind_rows(
  GSE25066_primary_clin,
  MAINZ_primary_clin,
  STK_primary_clin,
  TRANSBIG_primary_clin,
  MSK_primary_clin,
  UPP_primary_clin,
  VDX_primary_clin
)

clin_meta$dfs <- clin_meta$dfs / 30
clin_meta$dmfs <- clin_meta$dmfs / 30
clin_meta$os <- clin_meta$os / 30

################################################
#####full join clin_meta and metabric_clin######
###############################################

#--- Step 1: Add a source indicator to each data frame
clin_meta$source <- "clin_meta"
metabric_clin$source <- "metabric_clin"

#--- Step 2: Determine the common columns (including the key and source)
common_cols <- intersect(colnames(clin_meta), colnames(metabric_clin))
# (Make sure that your key column—for example, "sampleID"—is present in both.)

#--- Step 3: Merge by taking only the intersection (using rbind)
merged_clin <- rbind(clin_meta[, common_cols, drop = FALSE],
                     metabric_clin[, common_cols, drop = FALSE])

#--- Step 4: Identify extra columns in each file
extra_clin_meta <- setdiff(colnames(clin_meta), common_cols)
extra_metabric <- setdiff(colnames(metabric_clin), common_cols)

#--- Step 5: Add extra columns from clin_meta
for (col in extra_clin_meta) {
  # Create a new column in merged_clin and initialize with NA
  merged_clin[[col]] <- NA
  # For rows that came from clin_meta, fill in the values.
  # We match by the 'source' column.
  idx <- which(merged_clin$source == "clin_meta")
  # Use the original order from clin_meta (assuming rbind preserved the order)
  merged_clin[idx, col] <- clin_meta[[col]]
}

#--- Step 6: Add extra columns from metabric_clin
for (col in extra_metabric) {
  merged_clin[[col]] <- NA
  idx <- which(merged_clin$source == "metabric_clin")
  merged_clin[idx, col] <- metabric_clin[[col]]
}

merged_clin <- merged_clin %>%
  dplyr::select(sampleID, dfs, dfs_status, os, os_status, everything())

##########################
#####Optional Save########
##########################

write_tsv(clin_meta, "clin_meta.tsv")
