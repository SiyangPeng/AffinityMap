######### This file contains script used to process DIANN output data ##########

# Set global variables
working_directory <- r"(path/to/DIANN/processed/file)"
# Define output file name
output_file = 'data_preprocessed_separate.csv'
# String that is in all columns with signal intensity but no other columns
data_identifier <- "CW" # example
# Column names from DIANN output
gene_column <- "Genes"
protein_column <- "Protein.Ids"
pg_filename <- "report.pg_matrix.tsv"
pr_filename <- "report.pr_matrix.tsv"

# Load libraries
library("MSnbase")
library("MSstats")
library("limma")
library("DEP")
library("dplyr")

# Specify columns used for downstream processing in each treatment group
data_columns <- c(6:101)

# Separate by probe concentration
conc1_column_names <- c(6:29)
conc2_column_names <- c(30:53)
conc3_column_names <- c(54:77)
conc4_column_names <- c(78:101)

# Load DIANN-generated protein group file
data <- read.csv(sprintf("% s/% s", working_directory,pg_filename), sep = "\t")
# Generates unique identifiers for a proteomics dataset 
# gene_column = "name", protein_column = "id"
data_unique <- make_unique(data, gene_column, protein_column) 
# Rename index of rows to gene name
rownames(data_unique) <- data_unique[,"name"]

# Select metadata columns to be used 
names_column <- grep(data_identifier, colnames(data_unique))

# Load DIANN-generated peptide precursor file
# Extract number of peptides per protein from precursor matrix file
data_peptides <- read.csv(sprintf("% s/% s", working_directory,pr_filename), sep = "\t")
unique_peptides <- distinct(data_peptides)
unique_peptides <- unique_peptides[order(unique_peptides$Stripped.Sequence),]
data_peptides_counted <- add_count(unique_peptides, Protein.Ids)
data_peptides_counted <- distinct(data_peptides_counted, Genes, Protein.Ids, n)
data_peptides_counted_unique <- as.data.frame(make_unique(data_peptides_counted, "Genes", "Protein.Ids"))
colnames(data_peptides_counted_unique)[3] <- "Combined.Total.Peptides"
rownames(data_peptides_counted_unique) <- data_peptides_counted_unique[,"name"]
data_unique <- merge(data_unique, 
                     data_peptides_counted_unique[3],
                     by = 0)
data_unique <- data_unique[,-1]
rownames(data_unique) <- data_unique[,"name"]

# Format data for imputation and normalization
# Each probe concentration independently normalized and imputed
data_unique$NA_percent <- rowSums(apply(is.na(data_unique[,names_column]),2,as.numeric))/ncol(data_unique[,names_column])

# Calculate the number of missing values for each protein at each defined probe concentration
data_unique$NA_percent_conc1 <- rowSums(apply(is.na(data_unique[,conc1_column_names]),2,as.numeric))/ncol(data_unique[,conc1_column_names])
data_unique$NA_percent_conc2 <- rowSums(apply(is.na(data_unique[,conc2_column_names]),2,as.numeric))/ncol(data_unique[,conc2_column_names])
data_unique$NA_percent_conc3 <- rowSums(apply(is.na(data_unique[,conc3_column_names]),2,as.numeric))/ncol(data_unique[,conc3_column_names])
data_unique$NA_percent_conc4 <- rowSums(apply(is.na(data_unique[,conc4_column_names]),2,as.numeric))/ncol(data_unique[,conc4_column_names])
# Store in defined columns
NA_filter_columns = c("NA_percent_conc1", "NA_percent_conc2","NA_percent_conc3","NA_percent_conc4")

# Convert into treatment group-unique msnset objects for downstream processing
msnset_conc1 <- readMSnSet2(data_unique, conc1_column_names, "name")
msnset_conc2 <- readMSnSet2(data_unique, conc2_column_names, "name")
msnset_conc3 <- readMSnSet2(data_unique, conc3_column_names, "name")
msnset_conc4 <- readMSnSet2(data_unique, conc4_column_names, "name")

# Normalize each treatment group separately
msnset_conc1_norm <- normalise(msnset_conc1, "vsn")
msnset_conc2_norm <- normalise(msnset_conc2, "vsn")
msnset_conc3_norm <- normalise(msnset_conc3, "vsn")
msnset_conc4_norm <- normalise(msnset_conc4, "vsn")
  
# Impute each treatment group separately, generate matrix of imputed data
msnset_conc1_norm_imputed <- MSnbase::impute(msnset_conc1_norm, "MinProb")
msnset_conc2_norm_imputed <- MSnbase::impute(msnset_conc2_norm, "MinProb")
msnset_conc3_norm_imputed <- MSnbase::impute(msnset_conc3_norm, "MinProb")
msnset_conc4_norm_imputed <- MSnbase::impute(msnset_conc4_norm, "MinProb")

# Combine separately processed datasets into one msnset object
msnset_normed_imputed <- MSnbase::combine(msnset_conc1_norm_imputed,
                                          msnset_conc2_norm_imputed,
                                          msnset_conc3_norm_imputed,
                                          msnset_conc4_norm_imputed)

# Combine processed data and metadata of interest into dataframe to be used for further processing in python, output csv file
metadata_columns <- c('Protein.Ids','Protein.Names', 'Genes', 
                      'First.Protein.Description', "name", "ID", 
                      "Combined.Total.Peptides", "NA_percent")
data_combined_preprocessing_names <- msnset_normed_imputed@featureData@data[,metadata_columns]
data_combined_preprocessing_exprs <- 2^msnset_normed_imputed@assayData[["exprs"]]
data_combined_NAs <-msnset_normed_imputed@featureData@data[,NA_filter_columns]
data_preprocessed <- cbind(data_combined_preprocessing_names, data_combined_preprocessing_exprs,data_combined_NAs)
write.csv(data_preprocessed, file = sprintf("% s/% s", working_directory,output_file))
