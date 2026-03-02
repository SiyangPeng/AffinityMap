#Set working directory
working_directory_windows <- r"(path/to/DIANN/processed/file)"

filter_percent_cutoff <- .1
imputation_filter <- "strict"
peptide_count_cutoff <- 1
abundance_cutoff <- 0
data_identifier <- "Slot"
pos_identifier <- ""
neg_identifier <- ""
protein_of_interest_1 <- "name1"
protein_of_interest_2 <- "name2"
protein_of_interest_3 <- "name3"

#Define sample group columns
#1-3
group_a <- c(6,8,10)
#4-6   
group_b <- c(9,11,12)

#Define sample groups to be compared
pos_column_names <- group_a
neg_column_names <- group_b

#Volcano plot labeling cutoffs for plotting
PValue_cutoff = 1e-2
logFC_cutoff = -10

#Define data type. For msfragger, "DDA"; for diann, "DIA"
data_type <- "DIA"

# Define processing type: uniform (normalization and imputation done together) or separate (arms treated independently for normalization and imputation).

processing_type <- "separate_arms"

#Data structure definition
if(data_type == "DDA") {
  gene_column <- "Gene"
  protein_column <- "Protein.ID"
  metadata_columns <- 11
  filename <- "%scombined_protein.tsv"
}

if(data_type == "DIA") {
  gene_column <- "Genes"
  protein_column <- "Protein.Group"
  metadata_columns <- 5
  filename <- "%sreport.pg_matrix.tsv"
}

# Load libraries
library("MSnbase")
library("MSstats")
library("limma")
library("EnhancedVolcano")
library("DEP")
library("dplyr")

# Convert working directory to forward slash structure
working_directory <- gsub("[^A-Za-z0-9 _ -]", "/", working_directory_windows)
working_directory <- sub("/", ":", working_directory)
working_directory <- paste(working_directory, "/", sep = "")

# Load data, make unique entry column, assign NA as required (MSfragger reports 0)
data <- read.csv(sprintf(filename, working_directory), sep = "\t")
data_unique <- make_unique(data, gene_column, protein_column)
rownames(data_unique) <- data_unique[,"name"]
names_column <- grep(data_identifier, colnames(data_unique))
data_unique[,names_column][data_unique[,names_column] == 0] <- NA

#For DIA data, extract number of peptides per protein from precursor matrix file
if(data_type == "DIA") {
  data_peptides <- read.csv(sprintf("%sreport.pr_matrix.tsv", working_directory), sep = "\t")
  unique_peptides <- distinct(data_peptides)
  unique_peptides <- unique_peptides[order(unique_peptides$Stripped.Sequence),]
  data_peptides_counted <- add_count(unique_peptides, Protein.Group)
  data_peptides_counted <- distinct(data_peptides_counted, Genes, Protein.Group, n)
  data_peptides_counted_unique <- as.data.frame(make_unique(data_peptides_counted, "Genes", "Protein.Group"))
  colnames(data_peptides_counted_unique)[3] <- "Combined.Total.Peptides"
  rownames(data_peptides_counted_unique) <- data_peptides_counted_unique[,"name"]
  data_unique <- merge(data_unique, 
                       data_peptides_counted_unique[3],
                       by = 0)
  data_unique <- data_unique[,-1]
  rownames(data_unique) <- data_unique[,"name"]
}

# Rename sample group-containing columns 
colnames(data_unique)[pos_column_names] <- "high_a"
colnames(data_unique)[neg_column_names] <- "low_a"

# Define imputation cases for each protein, filter proteins based on missing value cutoff. Filtration of missing values is on a one-arm basis for separate arm processing,
# so that data is retained if one arm passes the cutoff and "all-negative" arms are retained. For uniform analysis missing values are filtered globally.
data_unique$NA_percent_neg <- rowSums(apply(is.na(data_unique[,neg_column_names]),2,as.numeric))/ncol(data_unique[,neg_column_names])
data_unique$NA_percent_pos <- rowSums(apply(is.na(data_unique[,pos_column_names]),2,as.numeric))/ncol(data_unique[,pos_column_names])
data_unique$fully_imputed <- as.numeric(data_unique[,"NA_percent_pos"] == 1 | data_unique[,"NA_percent_neg"] == 1)
data_unique$partially_imputed <- as.numeric(!data_unique[,"NA_percent_pos"] %in% c(0,1) | !data_unique[,"NA_percent_neg"] %in% c(0,1))
data_unique$not_imputed <- as.numeric(data_unique[,"NA_percent_pos"] == 0 & data_unique[,"NA_percent_neg"] == 0)
data_unique_filtered <- subset(data_unique, Combined.Total.Peptides > peptide_count_cutoff)

# Define imputation conditions: strict requires that both arms exceed the cutoff filter, lax dictates that only one needs to satisfy the cutoff filter
if(imputation_filter == "strict") {
  data_unique_filtered <- subset(data_unique_filtered, NA_percent_neg < filter_percent_cutoff & NA_percent_pos < filter_percent_cutoff)
}

if(imputation_filter == "lax") {
  data_unique_filtered <- subset(data_unique_filtered, NA_percent_neg < filter_percent_cutoff | NA_percent_pos < filter_percent_cutoff)
}


#For imputation and normalization of separate experimental arms:
if(processing_type == "separate_arms") {
  # Convert into treatment-group unique msnset objects for downstream processing
  msnset_neg <- readMSnSet2(data_unique_filtered, neg_column_names, "name")
  msnset_pos <- readMSnSet2(data_unique_filtered, pos_column_names, "name")
  
  # Normalize in each treatment group separately
  msnset_neg_norm <- normalise(msnset_neg, "vsn")
  msnset_pos_norm <- normalise(msnset_pos, "vsn")
  
  # Combine and generate matrix of unimputed data for optional inspection
  combined_preprocessing_postnorm <- MSnbase::combine(msnset_neg_norm, msnset_pos_norm)
  data_combined_preprocessing_postnorm_names <- combined_preprocessing_postnorm@featureData@data[1:metadata_columns]
  data_combined_preprocessing_postnorm_exprs <- combined_preprocessing_postnorm@assayData[["exprs"]]
  data_preprocessed_postnorm <- cbind(data_combined_preprocessing_postnorm_names, data_combined_preprocessing_postnorm_exprs)
  data_preprocessed_postnorm_unique <- make_unique(data_preprocessed_postnorm, gene_column, protein_column)
  write.csv(data_preprocessed_postnorm_unique, file = sprintf("%sreport_no_imputation.csv", working_directory))
  
  # Impute in each treatment group separately, combine into one msnset object, generate matrix of imputed data
  msnset_neg_norm_imputed <- MSnbase::impute(msnset_neg_norm, "MinProb")
  msnset_pos_norm_imputed <- MSnbase::impute(msnset_pos_norm, "MinProb")
  combined_preprocessing <- MSnbase::combine(msnset_neg_norm_imputed,
                                             msnset_pos_norm_imputed)
}

#For uniform imputation and normalization of experimental arms:
if(processing_type == "uniform") {
  msnset <- readMSnSet2(data_unique_filtered, names_column, "name")
  
  # Normalize globally
  msnset_norm <- normalise(msnset, "vsn")
  
  # Generate matrix of unimputed data
  data_combined_preprocessing_postnorm_names <- msnset_norm@featureData@data[1:metadata_columns]
  data_combined_preprocessing_postnorm_exprs <- msnset_norm@assayData[["exprs"]]
  data_preprocessed_postnorm <- cbind(data_combined_preprocessing_postnorm_names, data_combined_preprocessing_postnorm_exprs)
  data_preprocessed_postnorm_unique <- make_unique(data_preprocessed_postnorm, gene_column, protein_column)
  write.csv(data_preprocessed_postnorm_unique, file = sprintf("%sreport_no_imputation.csv", working_directory))
  
  # Impute missing values
  combined_preprocessing <- MSnbase::impute(msnset_norm, "MinProb")
}


#For processing without normalization and imputation:
if(processing_type == "none") {
  data_unique_filtered[,names_column] <- log(data_unique_filtered[names_column], 2)
  msnset <- readMSnSet2(data_unique_filtered, names_column, "name")
  
  # Pass dataset
  msnset_norm <- msnset
  data_combined_preprocessing_postnorm_names <- msnset_norm@featureData@data[1:metadata_columns]
  data_combined_preprocessing_postnorm_exprs <- msnset_norm@assayData[["exprs"]]
  data_preprocessed_postnorm <- cbind(data_combined_preprocessing_postnorm_names, data_combined_preprocessing_postnorm_exprs)
  data_preprocessed_postnorm_unique <- make_unique(data_preprocessed_postnorm, gene_column, protein_column)
  write.csv(data_preprocessed_postnorm_unique, file = sprintf("%sreport_no_imputation.csv", working_directory))
  
  # Pass dataset
  combined_preprocessing <- msnset_norm
}

# Combine and generate matrix of imputed data for optional inspection
data_combined_preprocessing_names <- combined_preprocessing@featureData@data[1:metadata_columns]
data_combined_preprocessing_exprs <- combined_preprocessing@assayData[["exprs"]]
data_preprocessed <- cbind(data_combined_preprocessing_names, data_combined_preprocessing_exprs)
data_preprocessed_unique <- make_unique(data_preprocessed, gene_column, protein_column)
write.csv(data_preprocessed_unique, file = sprintf("%sreport_imputed.csv", working_directory))

# Find data in imputed matrix, convert into suitable format for limma, create limma design matrix
data_limma <- data.matrix(data_preprocessed_unique[c(grep("low_a", colnames(data_preprocessed_unique)), 
                                                     grep("high_a", colnames(data_preprocessed_unique)))])
colnames <- 1:sum(length(pos_column_names), length(neg_column_names))
colnames(data_limma) <- 1:sum(length(pos_column_names), length(neg_column_names))
design <- cbind(grp1=1,grp2=c(rep(0, length(neg_column_names)),rep(1, length(pos_column_names))))

# Run limma 
limma_output <- eBayes(lmFit(data_limma, design))

# Extract all DEP data from limma output, combine with descriptive columns and generate export matrix.
# For DIA data incorporate a column for the total peptides used for each protein.
limma_output_toptable <- topTable(limma_output, coef=2, n = 10000)
limma_output_toptable$neg_log_p <- -log10(limma_output_toptable[,"P.Value"])
limma_output_toptable <- cbind(row.names(limma_output_toptable), limma_output_toptable)
colnames(limma_output_toptable)[1] <- "name"
full_report <- merge(limma_output_toptable, 
                     data_unique_filtered[grep("imputed", colnames(data_unique_filtered))],
                     by = 0)
full_report <- full_report[,-1]
rownames(full_report) <- full_report[,"name"]
if(data_type == "DIA") {
  full_report <- merge(full_report, 
                       data_unique_filtered[grep("Combined.Total.Peptides", colnames(data_unique_filtered))],
                       by = 0)
  full_report <- full_report[,-1]
  colnames(full_report)[1] <- "name"
  rownames(full_report) <- full_report[,"name"]
}
full_report <- merge(full_report, 
                     data_preprocessed_unique[c(grep("low_a", colnames(data_preprocessed_unique)),grep("high_a", colnames(data_preprocessed_unique)))],
                     by = 0)
full_report <- full_report[,-1]
colnames(full_report)[1] <- "name"
rownames(full_report) <- full_report[,"name"]
full_report <- merge(full_report, 
                     data_unique_filtered[,1:metadata_columns],
                     by = 0)
full_report <- full_report[,-1]
rownames(full_report) <- full_report[,"name"]
colnames(full_report)[1] <- "name"
full_report <- full_report[order(full_report$P.Value),]
full_report <- subset(full_report, AveExpr > abundance_cutoff)
write.csv(full_report, file = sprintf("%sDEA_report_imputed.csv", working_directory))

# Reorder report to put POIs on bottom, Define volcano plot custom colors, labels, and cutoffs
proteins_of_interest <- c(protein_of_interest_1, protein_of_interest_2, protein_of_interest_3)
poi_positions <- c(grep(protein_of_interest_1, full_report[,"name"]), 
                   grep(protein_of_interest_2, full_report[,"name"]), 
                   grep(protein_of_interest_3, full_report[,"name"]))

if(!isEmpty(poi_positions)) {
  proteins_of_interest <- c(protein_of_interest_1, protein_of_interest_2, protein_of_interest_3)
  rownames(full_report) <- c()
  full_report_no_poi <- full_report[-poi_positions,]
  full_report_no_poi$is_poi<- 0
  full_report_poi <- full_report[poi_positions,]
  full_report_poi$is_poi<- 1
  full_report <- rbind(full_report_no_poi, full_report_poi)
  poi_positions <- c(grep(protein_of_interest_1, full_report[,"name"]), 
                     grep(protein_of_interest_2, full_report[,"name"]), 
                     grep(protein_of_interest_3, full_report[,"name"]))
  keyvals <- ifelse(full_report$name %in% proteins_of_interest, "orange3",
                    ifelse(full_report$fully_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick1",
                           ifelse(full_report$partially_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick3",
                                  ifelse(full_report$not_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff,  "darkred",
                                         "black"))))
  
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == "orange3"] <- "proteins_of_interest"
  names(keyvals)[keyvals == 'firebrick1'] <- 'fully imputed'
  names(keyvals)[keyvals == 'firebrick3'] <- 'partially imputed'
  names(keyvals)[keyvals == 'darkred'] <- 'not imputed' 
}

if(isEmpty(poi_positions)) {
  full_report$is_poi<- 0
  keyvals <- ifelse(full_report$fully_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick1",
                    ifelse(full_report$partially_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick3",
                           ifelse(full_report$not_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff,  "darkred",
                                  "black")))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == "orange3"] <- "proteins_of_interest"
  names(keyvals)[keyvals == 'firebrick1'] <- 'fully imputed'
  names(keyvals)[keyvals == 'firebrick3'] <- 'partially imputed'
  names(keyvals)[keyvals == 'darkred'] <- 'not imputed'  
}

plot_max_logfc <- max(full_report$logFC) + 0.5
plot_max_logp <- max(full_report$neg_log_p) + 0.5
plot_min_logfc <- min(full_report$logFC) -0.5 

# Plot volcano plot
EnhancedVolcano(full_report, 
                lab = full_report$name,
                x = 'logFC',
                y = 'P.Value', 
                xlim = c(plot_min_logfc, plot_max_logfc),
                ylim = c(0,plot_max_logp),
                col=c('black', 'black', 'black', 'red3'),
                FCcutoff = logFC_cutoff,
                pCutoff = PValue_cutoff,
                labSize = 4,
                pointSize = c(ifelse(full_report$is_poi>0, 3, 3)),
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                arrowheads = FALSE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                subtitle = NULL,
                legendPosition = "none",
                title = NULL,
                cutoffLineType = "blank",
                colCustom = keyvals,
                selectLab = c((full_report$name)[which(names(keyvals) %in% c('fully imputed', 'partially imputed', 'not imputed', "proteins_of_interest"))]))

# Export volcano plot as .tif file
tiff(sprintf("%sRPlot.tif", working_directory), width = 6000, height = 6000, res = 600, compression = c("lzw"))
EnhancedVolcano(full_report, 
                lab = full_report$name,
                x = 'logFC',
                y = 'P.Value', 
                xlim = c(plot_min_logfc, plot_max_logfc),
                ylim = c(0,plot_max_logp),
                col=c('black', 'black', 'black', 'red3'),
                FCcutoff = logFC_cutoff,
                pCutoff = PValue_cutoff,
                labSize = 4,
                pointSize = c(ifelse(full_report$is_poi>0, 3, 3)),
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                arrowheads = FALSE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                subtitle = NULL,
                legendPosition = "none",
                title = NULL,
                cutoffLineType = "blank",
                colCustom = keyvals,
                selectLab = c((full_report$name)[which(names(keyvals) %in% c('fully imputed', 'partially imputed', 'not imputed', "proteins_of_interest"))]))
dev.off()


